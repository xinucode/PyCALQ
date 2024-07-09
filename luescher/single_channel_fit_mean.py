import logging
import general.data_reader as dr
import general.plotting_handler as ph
import fvspectrum.sigmond_util as sigmond_util

############################################################################################
##      Required packages
import cmath                    #math library (complex)
import csv                      # read csv files
import math                     # math library, mostly use it for math.pi, math.pow()
import matplotlib.pyplot as plt #plotting library            
import numpy as np              #basic functions, linear algebra, etc. 
import random
import scipy.special as sp      #special functions
import scipy.integrate as integrate #needed to do numerical integration, alt?
import yaml
import os
import luescher.tools.parametrizations as parametrizations 
import luescher.tools.kinematics as kinematics
from luescher.tools.zeta import Z
##
from matplotlib.lines import Line2D  # Import Line2D for custom legend handles
#from tqdm import tqdm  # Import the tqdm library
from scipy import optimize,integrate
from scipy.optimize import fsolve,minimize
from scipy.integrate import quad


## Z function that is fast


doc = '''
essential documentation

tasks param should have which channel to analyze
as well as how many levels for each data to use 
by default it is one level for each state
'''

class SingleChannelFitMean:

    @property
    def info(self):
        return doc
    

    def __init__( self, task_name, proj_handler, general_configs, task_params ):
        print("Task parameters:",task_params)
        self.proj_handler = proj_handler 
        self.task_name = task_name
        # if not task_configs:
        #     logging.critical(f"No directory to view. Add 'raw_data_files' to '{task_name}' task parameters.")

        ensemble_info = sigmond_util.get_ensemble_info(general_configs)
        self.L = ensemble_info.getLatticeXExtent()
        # path name needs to link toward the hdf5 file 
        if 'data_file' in task_params.keys(): # want to add if statement 
            file_path = task_params['data_file']

        print(file_path)
        print(self.L)
        self.dr = dr.LQCD_DATA_READER(file_path,self.L) # data reader is object dr
        self.data = self.dr.load_data()
        print("Data loaded")
        #print(self.data)

        # check that data from Hf file is real
        if not self.data:
            logging.critical(f"Ensure data has been generated for '{task_name}' to continue single-channel analysis.")

        self.alt_params = {
            'write_data': True,
            'create_pdfs': True,
            'plot': True,
            'figwidth':8,
            'figheight':6,
        }
         #Sarah
        self.irreps = self.alt_params['irreps']
        if type(self.irreps)!=dict:
            logging.error("""Incorrect setup for irreps. Expecting definition of form:
            irreps:
                PSQ0:
                - G1u
                - G1g
                - ...
                PSQ1:
                ...""")

        #hadron list
        self.single_hadron_list = np.array(self.dr.single_hadron_list())

        self.channel = task_params['channel']
        self.channel_1 = self.channel[0]
        self.channel_2 = self.channel[1]

        if np.any(np.isin(self.single_hadron_list,self.channel_1)) and np.any(np.isin(self.single_hadron_list,self.channel_2)):
            logging.info("Channel is confirmed to be in data file. Continuing analysis ...")
            
        #initialize your task, store default input in self.proj_dir_handler.log_dir() (basically, throw the full possible input with all parameters where all the assumed parameters have been filled in in there)
        with open( os.path.join(proj_handler.log_dir(), 'full_input.yml'), 'w+') as log_file:
            yaml.dump({"general":general_configs}, log_file)
            yaml.dump({"tasks":task_params}, log_file)


    def momentum_state(self,i): #
        # d = [0,0,0]
        if i == 'PSQ0':
            return np.array([0,0,0])
        # d = [0,0,1]
        elif i == 'PSQ1':
            return np.array([0,0,1])
        # d = [1,1,0]
        elif i == 'PSQ2':
            return np.array([1,1,0])
        # d = [1,1,1]
        elif i == 'PSQ3':
            return np.array([1,1,1])
        # d = [0,0,2]
        elif i == 'PSQ4':
            return np.array([0,0,2])
        else: 
            # Raise an exception for invalid input
            raise ValueError("Invalid value for 'i'. 'i' must be 0, 1, 2, or 3.")
    
    
    def run( self ):       
        log_path = os.path.join(self.proj_handler.log_dir(), 'luescher_log.yml') 
        # step 1, import the keys needed for analysis from single_hadron_list
        # first get the required channel masses
        ma_ref = self.dr.single_hadron_data(self.channel[0]) # make sure using ref masses
        mb_ref = self.dr.single_hadron_data(self.channel[1]) # 
        ref_mass = self.dr.single_hadron_data('ref')
        # next load in the data
        # lowest level from each data
        # if psq != 0, do next level instead
        psq_list = self.dr.load_psq()
        # collect irreps for each of the PSQ, which is the available data
        irreps_all = {}
        for psq in psq_list:
            irreps_all[psq] = []
            for key in self.dr.irrep_keys(psq):
                irreps_all[psq].append(key)

        irreps = self.irreps #{'PSQ0': ['G1u'],'PSQ1': ['G1'],'PSQ2': ['G'], 'PSQ3': ['G']}
        print(irreps)
        # irreps = {'PSQ0': ['G1u'],'PSQ1': ['G1'],'PSQ2': ['G'], 'PSQ3': ['G']}
        # print(irreps)
        psq_remove = []
        for psq in psq_list:
            if psq not in irreps:
                psq_remove.append(psq)
        for psq in psq_remove:
            psq_list.remove(psq)

        self.ecm_data = {} # save possible data
        for psq in psq_list:
            self.ecm_data[psq] = {}
            for irrep in irreps[psq]:
                if psq == 'PSQ0':
                    level = 0
                    # level = "ecm_0_ref" # level as an int
                else:
                    level = 1
                    # level = "ecm_1_ref" 
                # self.ecm_data[psq][irrep] = self.data.get(psq).get(irrep).get(level)[:]
                self.ecm_data[psq][irrep] = self.dr.ref_energy_data(psq,irrep,level)
        # now that we have ecm_data, lets do a fit and save results
        self.average_energies = []
        #first we want to use average data
        self.ecm_average_data = {}
        for psq in psq_list:
            self.ecm_average_data[psq] = {}
            for irrep in irreps[psq]:
                self.ecm_average_data[psq][irrep] = self.ecm_data[psq][irrep][0]
                self.average_energies.append(self.ecm_data[psq][irrep][0])

        self.ma_ave = ma_ref[0]
        self.mb_ave = mb_ref[0]
        self.ref_ave = ref_mass[0]
        #print(self.ref_ave)

        # get cov data
        self.ecm_bs = {}
        for psq in psq_list:
            self.ecm_bs[psq] = {}
            for irrep in irreps[psq]:
                self.ecm_bs[psq][irrep] = self.ecm_data[psq][irrep][1:]

        ecm_NN_bs_arr = []
        for psq in psq_list:
            for irrep in irreps[psq]:
                ecm_NN_bs_arr.append(self.ecm_bs[psq][irrep].tolist()) # the ordering should be kep tin which the calculation will go
        


        def energy_shift_data():

            ecm_data = ecm_NN_bs_arr

            #Sarah
            mref = np.array(self.dr.single_hadron_data('ref'))[1:]
            m1_ref = np.array(self.dr.single_hadron_data(self.channel_1))
            m2_ref = np.array(self.dr.single_hadron_data(self.channel_2))
            # mpi = np.array(self.dr.single_hadron_data('pi(0)'))[1:]
            # mS_ref = np.array(self.dr.single_hadron_data('S(0)_ref'))
            # mpi_ref = np.array(self.dr.single_hadron_data('pi(0)_ref'))
            # mk_ref = np.array(self.dr.single_hadron_data('k(0)_ref'))
            # mN_ref = np.array(self.dr.single_hadron_data('N(0)_ref'))
            
           # mapping for masses
            mass_map = {
                get_particle_name(self.channel_1): m1_ref[1:],
                get_particle_name(self.channel_2): m2_ref[1:],
                # 'pi': np.array(mpi_ref)[1:],
                # 'S': np.array(mS_ref)[1:],
                # 'k': np.array(mk_ref)[1:],
                # 'N': np.array(mN_ref)[1:],
            }
            
            def extract_values(input_str):
                # Find the position of the opening and closing parentheses
                open_paren = input_str.find('(')
                close_paren = input_str.find(')')
                
                if open_paren != -1 and close_paren != -1:
                    # Extract the part before the opening parenthesis
                    part_before_paren = input_str[:open_paren] #Mass from mass_map
                    
                    # Extract the number inside the parentheses and convert it to an integer
                    number_inside_paren = int(input_str[open_paren + 1:close_paren]) # d^2
                    
                    return part_before_paren, number_inside_paren
                else:
                    return None, None  # Return None for both values if parentheses are not found

            def deltaE(ecm,ma,mb,n,m,psq):#function to shift e_cm data to shifted data to free energy levels
                if psq == 0:
                    l = self.L*mref #Sarah
                    dE = ecm - np.sqrt(ma**2 + n*(2*math.pi/l)**2 ) - np.sqrt((mb)**2 + m*(2*math.pi/l)**2 )
                else:
                    l = self.L*mref #Sarah
                    elab  = np.sqrt((ma)**2 + n*(2*math.pi/l)**2) + np.sqrt((mb)**2 + m*(2*math.pi/l)**2)
                    ecmfree = np.sqrt(elab**2 - psq*(2*math.pi/(l))**2) 
                    dE = ecm - ecmfree
                return dE
            
            mom2 = ['PSQ0','PSQ1','PSQ2','PSQ3']
            # Define a mapping of mom values to irreps
            irreps = {
                'PSQ0': 'G1u',
                'PSQ1': 'G1',
                'PSQ2': 'G',
                'PSQ3': 'G',
            }
            lvls = [1,1,1,1]
            data_list = []
            total_label_num = -1  # Initialize total label number
            # now need to get the energy shift from two-particle free energy
            # these are given in hdf5 file attributes
            #print(hf['PSQ3/G'].attrs['free_levels']) example of free energy shifts
            for mom, lvl in zip(mom2, lvls):
                # Determine the associated irrep based on mom using the mapping
                irr = irreps.get(mom)

                # Determine the range of label numbers based on the level
                if mom == 'PSQ0':
                    # For PSQ0, use the full range including 0
                    label_number = 0
                else:
                    # For other mom values, skip 0 and start from 1
                    label_number = 1

                total_label_num += 1
                # get energy from ecm_data
                ma, n = extract_values(self.dr.free_levels(mom,irr,label_number)[0])
                mb, m = extract_values(self.dr.free_levels(mom,irr,label_number)[1])
                ma = mass_map[ma]
                mb = mass_map[mb]
                #data_array = self.energy_data(mom,irr,label)
                data_list.append(deltaE(ecm_data[total_label_num],ma,mb,n,m,int(mom[3:])))
            
            data = np.array(data_list)
            return data

        self.covariance_matrix = np.cov(np.array(ecm_NN_bs_arr))
        self.cov_de = np.cov(energy_shift_data())

        def determinant_condition(ecm,psq,a,b):
            #p = psq[3]
            return (kinematics.qcotd(ecm,self.L,psq,self.ma_ave,self.mb_ave,self.ref_ave) - parametrizations.ere_delta(ecm,self.ma_ave,self.mb_ave,a,b))

        def QC1(psq,irrep, a, b):
            func = lambda ecm: determinant_condition(ecm,psq, a, b)
            return fsolve(func,self.ecm_average_data[psq][irrep])[0] #guess is the energy going in

        def chi2(x):
            # if len(x) == 1:
            #     a = x
            # else:
            #     a,b = x
            a,b = x
            res = []
            for psq in psq_list:
                for irrep in irreps[psq]:
                    # print(psq)
                    # print(irrep)
                    #d = psq[3]
                    #print(self.ecm_average_data[psq][irrep])
                    diff = self.ecm_average_data[psq][irrep] - QC1(psq,irrep,a,b)
                    res.append(diff)
            value = np.array(res)@np.linalg.inv(self.cov_de)@np.array(res)
            return value

        print(f" Using Effective Range Expansion for single-channel:{self.channel_1} (m = {self.ma_ave}) ,{self.channel_2} (m = {self.mb_ave})")


        # next lets run a fit to check
        #logging.info(r" Minimizing ERE for average data set")
        def average_fit(): #~~~ returning the mean bootstrap sample fit 
            result = minimize(chi2,x0=[0.04,0.6],method='nelder-mead')
            print(result)
            return [result.x[0],result.x[1],result.fun]
        # def average_fit():#output_file=f"fit_results_channel_{self.channel_1}{self.channel_2}.txt"):
        #     # will do a 1 and 2 parameter fit and see which is best
        #     # best_a = None
        #     # best_b = None
        #     # best_chi2 = 
        #     try:
        #         # with open(output_file, 'w') as file: # 2 parameter ERE
        #         #     retry = True
        #         #     count = 0
        #         #     bigger_count = 0
        #         #     while retry:
        #         #         try:
        #                     # if count == 0:
        #                     #     a_guess = best_a_0
        #                     # else:
        #                     #     a_guess = np.random.uniform(-20,20)
        #         a_guess = np.random.uniform(-10,10)
        #         b_guess = np.random.uniform(-4,4)
        #         result = minimize(chi2, x0 = [a_guess , b_guess], method = 'nelder-mead')
        #         if result.success:
        #             best_a = result.x[0]
        #             best_b = result.x[1]
        #             chi2_res = result.fun
        #             # chi2_temp = result.fun
        #             # if count == 0:
        #             #     best_a = result.x[0]
        #             #     best_b = result.x[1]
        #             #     chi2_res = result.fun
        #             #     count += 1
        #             # else:
        #             #     if chi2_temp - chi2_res < 1e-2: # if chi2 not changing with different guess, its robust
        #             #         best_a = result.x[0]
        #             #         best_b = result.x[1]
        #             #         chi2_res = chi2_temp
        #             #         retry = False

        #             #     # elif chi2_temp < chi2_res : #if smaller chi2, try again
        #             #     #     best_a = a_guess
        #             #     #     best_b = b_guess
        #             #     #     chi2_res = chi2_temp
        #             #     #     count += 1
        #         else:
        #             pass
        #                         # retry with new parameters
        #                 # except Exception as e:
        #                 #     logging.error("Error:",e)

        #         print(f"Minimization for 2-parameter ERE: a = {best_a}, b = {best_b}, chi2 = {chi2_res} ")
        #         #file.write(f"{best_a} {best_b} {chi2_res}\n")

        #     except Exception as e:
        #         logging.error(f"Error writing to file:{e}")

        #     #logging.info(f" ERE fits written to file: {output_file}")                            
        #     return #[best_a, best_b, chi2_res]
            

        def deriv(n,i,x): # 2 parameter difference derivative
            if len(self.best_fit) == 1:
                a = self.best_fit
                b = 0
            else:
                a, b = self.best_fit

            eps = .001
            if n == 0:
                return (self.QC1(i,0,a-eps,b)-self.QC1(i,0,a+eps,b))/(2*eps)
            elif n ==1:
                return (self.QC1(i,0,a,b-eps)-self.QC1(i,0,a,b+eps))/(2*eps)
        
        def vij(self,x): #V_ij is error matrix for parameters, take _ii to get each error
            # v_nm = dp_i / dp_n
            # error estimation function
            # cov = self.cov = self.data.covariance_data()
            if len(self.best_fit) == 2:
                nint = [0,1]
            else: 
                nint = 0

            psq = [0,1,2,3]
            lmat = np.empty(0)
            for n in nint:
                if n == 0:
                    dl = []
                    for i in psq:
                        dl.append(self.deriv(n,i,x))

                    lmat = np.append(lmat,dl)
                else:
                    dl =[]
                    for i in psq:
                        dl.append(self.deriv(n,i,x))
                    
                    lmat = np.vstack([lmat,dl])
                
            Vnm = np.linalg.inv(lmat@np.linalg.inv(self.covdE)@np.transpose(lmat))
            
            return Vnm

        self.best_fit = average_fit()
        return self.best_fit#print(self.best_fit)
        # do the task, produce the data, data goes in self.proj_dir_handler.data_dir(), info/warning/errors about the process goes in self.proj_dir_handler.log_dir() (if any)

    def plot( self ):
        # x = []
        # x_range = []
        # y = []
        # y_range = []
        # for i in range(len(self.average_energies)):
        #     x.append(self.q2(self.average_energies[i], self.ma_ave, self.mb_ave))
        #     y.append(self.qcotd(self.average_energies[i], i, self.ma_ave, self.mb_ave))
        #     if i == 0:
        #         x_range_i = []
        #         y_range_i = []
        #         for en in np.linspace(self.average_energies[i] - 0.01, self.average_energies[i] + 0.01, 100):
        #             x_range_i.append(self.q2(en, self.ma_ave, self.mb_ave))
        #             y_range_i.append(self.qcotd(en, i, self.ma_ave, self.mb_ave))
        #         x_range.append(x_range_i)
        #         y_range.append(y_range_i)
        #     else:
        #         xp = []
        #         yp = []
        #         for en in np.linspace(self.average_energies[i] - 0.01, self.average_energies[i] + 0.01, 100):
        #             xp.append(self.q2(en,  self.ma_ave, self.mb_ave))
        #             yp.append(self.qcotd(en, i,  self.ma_ave, self.mb_ave))
        #         x_range.append(xp)
        #         y_range.append(yp)

        # ecm_values = np.linspace(min(self.average_energies)-2,max(self.average_energies)+2,1000)
        # q2_values = np.empty(0)
        # for en in ecm_values:
        #     q2_values = np.append(q2_values, self.q2(en,self.ma_ave, self.mb_ave))
        
        # virtual_state = []
        # for q2 in np.arange(-1,0.1,150):
        #     virtual_state.append(cmath.sqrt(-q2))

        # if len(self.best_fit) == 1:
        #     a = self.best_fit
        #     b = 0 
        # elif len(self.best_fit) == 2:
        #     a, b = self.best_fit
        # else:
        #     logging.critical("Best fit parameters are not generated")
        
        # def ere(ecm,a,b):
        #     return ((-1/a)+0.5*b*self.q2(ecm,self.ma_ave,self.mb_ave)) #in units of reference mass, usually mpi
        
        # best_fit_line = []
        # for q2 in q2_values:
        #     best_fit_line.append(ere(ecm_values,a,b))

        # plt.figure(figsize=(8,6))
        # shapes = ['o','s','D','v']
        # labels = ['$G_{1u}(0)$','$G_1 (1)$', '$G (2)$' , '$G(3)$']
        # legend_handles = []  # Create an empty list to store custom legend handles
        # for i in range(len(self.average_energies)):
        #     plt.plot(x[i], y[i], marker=shapes[i], color='blue', label=labels[i])
        #     plt.plot(x_range[i], y_range[i], color="blue", alpha=0.8)  # Plot the ranges with transparency
        #     # Add a custom legend handle (marker with no line)
        #     legend_handles.append(Line2D([0], [0], marker=shapes[i], color='w', markerfacecolor='blue', markersize=10, label=labels[i]))


        # plt.plot(q2_values,virtual_state,color='black',linestyle='--')
        # plt.plot(q2_values,best_fit_line , color='blue',linestyle='-.')
        # plt.axhline(y=0,color='black')
        # plt.axvline(x=0,color='black')
        # plt.ylim(0,0.6)
        # plt.xlim(-0.25,0.05)
        # # Customize the legend with custom handles (markers only)
        # legend = plt.legend(handles=legend_handles, loc='upper left', title='Legend', prop={'size': 12})
        # plt.xlabel("$q^{*2} / m_{\pi}^2$",fontsize=16, fontdict={'fontweight': 'bold', 'fontstyle': 'italic'})
        # plt.ylabel("$q^{*} / m_{\pi} \cot \delta $",fontsize=16, fontdict={'fontweight': 'bold', 'fontstyle': 'italic'})
        # plt.title(f' {self.channel_1} {self.channel_2} Scattering ',fontsize=16)
        # plt.savefig(f"{self.channel_1} {self.channel_2} Scattering.pdf")
        return "Done"
        
        # make the plots, store in self.proj_dir_handler.plot_dir(), again, any log/error warnings go in self.proj_dir_handler.log_dir() as well (if any)
def get_particle_name(particle_str):
    return particle_str.split("(")[0]