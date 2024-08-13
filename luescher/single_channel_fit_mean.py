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

        # retrieve data
        self.dr = dr.LQCD_DATA_READER(file_path,self.L) # data reader is object dr
        self.data = self.dr.load_data()
        logging.info("Data loaded")


        # check that data from Hf file is real
        if not self.data:
            logging.critical(f"Ensure data has been generated for '{task_name}' to continue single-channel analysis.")

        self.alt_params = {
            'write_data': True,
            'create_pdfs': True,
            'plot': True,
            'figwidth':8,
            'figheight':6,
            'ref_energies' : True, #trigger for using reference energies in analysis, currently default and support is TRUE
        }
        #Sarah
        self.channels_and_irreps = task_params['channels']#self.alt_params['irreps']
        # {['pi(0)_ref','S(0)_ref'] : {'PSQ0':['G1u'],'PSQ1':['G1'],'PSQ2':['G'],'PSQ3':['G']}}
        # if type(self.irreps)!=dict:
        #     logging.error("""Incorrect setup for irreps. Expecting definition of form:
        #     Isospin, Strangeness: #see data_reader
        #         PSQ0:
        #         - G1u:
        #             levels [0,1,...] 
        #         - G1g
        #         - ...
        #         PSQ1:
        #         ...""")

        #hadron list
        self.single_hadron_list = np.array(self.dr.single_hadron_list())
        # generate list of channels
        self.channel = []
        self.irreps = {}
        for channel in self.channels_and_irreps:
            self.channel.append(channel)
            self.irreps[channel] = self.channels_and_irreps[channel]
        print('channels',self.channel)
        print("irreps",self.irreps)
        for channel in self.channel:
            #print(channel)
            channel_1 = str(channel.split(',')[0])
            channel_2 = str(channel.split(',')[1])
            if np.any(np.isin(self.single_hadron_list,channel_1)) and np.any(np.isin(self.single_hadron_list,channel_2)):
                logging.info(f"scattering Channel {channel} is confirmed to be in data file. Continuing analysis ...")
            
        #initialize your task, store default input in self.proj_dir_handler.log_dir() (basically, throw the full possible input with all parameters where all the assumed parameters have been filled in in there)
        with open( os.path.join(self.proj_handler.log_dir(), 'full_input.yml'), 'w+') as log_file:
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
        self.m_ref_dict = {}
        self.ref_mass = {}
        self.psq_list = {}
        for i, channel in enumerate(self.channel):
            channel_1 = str(channel.split(',')[0])
            channel_2 = str(channel.split(',')[0])
            self.m_ref_dict[channel] = [self.dr.single_hadron_data(channel_1),self.dr.single_hadron_data(channel_2) ] #ma_ref, mb_ref in channel
            self.ref_mass[channel] = self.dr.single_hadron_data('ref')
            self.psq_list[channel] = self.dr.load_psq()
        # ma_ref = self.dr.single_hadron_data(self.channel[0]) # make sure using ref masses
        # mb_ref = self.dr.single_hadron_data(self.channel[1]) # 
        print("psq list",self.psq_list)
        # irreps_all = {}
        # irreps = {}
        # for channel in self.channel:
        #     irreps_all[channel] = {}
        #     irreps[channel] = {}
        #     for psq in psq_list[channel]:
        #         irreps_all[channel][psq] = []
        #         irreps[channel][psq] = []
        #         for key in self.dr.irrep_keys(psq):
        #             irreps_all[channel][psq].append(key)
        # print("irreps",irreps)
        
        # for channel in self.channel:
        #     irreps = self.irreps #{'PSQ0': ['G1u'],'PSQ1': ['G1'],'PSQ2': ['G'], 'PSQ3': ['G']}

        #logging.info(irreps )
        # irreps = {'PSQ0': ['G1u'],'PSQ1': ['G1'],'PSQ2': ['G'], 'PSQ3': ['G']}
        # print(irreps)
        # psq_remove = []
        # for psq in psq_list:
        #     if psq not in self.irreps:
        #         psq_remove.append(psq)
        # for psq in psq_remove:
        #     psq_list.remove(psq)

        self.ecm_data = {} # save possible data
        self.ecm_average_data = {}
        self.ecm_bootstrap_data = {}
        ecm_NN_bs_arr = {}
        for channel in self.channel:
            self.ecm_data[channel] = {} 
            self.ecm_average_data[channel] = {}
            self.ecm_bootstrap_data[channel] = {}
            ecm_NN_bs_arr[channel] = []
            for psq in self.psq_list[channel]:
                self.ecm_data[channel][psq] = {}
                self.ecm_average_data[channel][psq] = {}
                self.ecm_bootstrap_data[channel][psq] = {} 
                for irrep in self.irreps[channel][psq][0]: #need to add irrep list to data reader, [0] required to take out of list
                    self.ecm_data[channel][psq][irrep] = {}
                    self.ecm_average_data[channel][psq][irrep] = {} 
                    self.ecm_bootstrap_data[channel][psq][irrep] = {} 
                    for level in self.irreps[channel][psq][0][irrep]:
                        if self.alt_params['ref_energies']:
                            level_title = f"ecm_{level}_ref"
                            self.ecm_data[channel][psq][irrep][level_title] = self.dr.ref_energy_data(psq,irrep,level_title)
                            self.ecm_average_data[channel][psq][irrep][level_title] = self.ecm_data[channel][psq][irrep][level_title][0]
                            self.ecm_bootstrap_data[channel][psq][irrep][level_title] = self.ecm_data[channel][psq][irrep][level_title][1:]
                            ecm_NN_bs_arr[channel].append(self.ecm_bootstrap_data[channel][psq][irrep][level_title])
                        else:
                            level_title = f"ecm_{level}"
                            logging.critical("Need non-ref energies from data reader")


        def energy_shift_data(channel):
                
            ecm_data = self.ecm_bootstrap_data[channel]

            #Sarah
            mref = np.array(self.ref_mass[channel])[1:] #np.array(self.dr.single_hadron_data('ref'))[1:]
            m1_ref,m2_ref = self.m_ref_dict[channel]
            # m1_ref = np.array(self.dr.single_hadron_data(self.channel_1))
            # m2_ref = np.array(self.dr.single_hadron_data(self.channel_2))
            channel_1 = str(channel.split(',')[0])
            channel_2 = str(channel.split(',')[1])
           # mapping for masses
            mass_map = {
                get_particle_name(channel_1): m1_ref[1:],
                get_particle_name(channel_2): m2_ref[1:],
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
                    l = self.L*mref 
                    dE = ecm - np.sqrt(ma**2 + n*(2*math.pi/l)**2 ) - np.sqrt((mb)**2 + m*(2*math.pi/l)**2 )
                else:
                    l = self.L*mref #Sarah
                    elab  = np.sqrt((ma)**2 + n*(2*math.pi/l)**2) + np.sqrt((mb)**2 + m*(2*math.pi/l)**2)
                    ecmfree = np.sqrt(elab**2 - psq*(2*math.pi/(l))**2) 
                    dE = ecm - ecmfree
                return dE

            data_list = []
            for psq in self.psq_list[channel]:
                for irrep in self.irreps[channel][psq][0]:
                    for level in self.irreps[channel][psq][0][irrep]: # level is number
                        if self.alt_params['ref_energies']:
                            level_title = f"ecm_{level}_ref"
                        else:
                            level_title = f"ecm_{level}"
                        ma, n = extract_values(self.dr.free_levels(psq,irrep,level)[0])
                        mb, m = extract_values(self.dr.free_levels(psq,irrep,level)[1])
                        ma = mass_map[ma]
                        mb = mass_map[mb]
                        data_list.append(deltaE(ecm_data[psq][irrep][level_title] ,ma,mb,n,m,int(psq[3:])))

            # #total_label_num = -1  # Initialize total label number
            # # now need to get the energy shift from two-particle free energy
            # # these are given in hdf5 file attributes
            # #print(hf['PSQ3/G'].attrs['free_levels']) example of free energy shifts
            # for mom, lvl in zip(mom2, lvls):
            #     # Determine the associated irrep based on mom using the mapping
            #     irr = irreps.get(mom)

            #     # Determine the range of label numbers based on the level
            #     if mom == 'PSQ0':
            #         # For PSQ0, use the full range including 0
            #         label_number = 0
            #     else:
            #         # For other mom values, skip 0 and start from 1
            #         label_number = 1

            #     total_label_num += 1
            #     # get energy from ecm_data
            #     ma, n = extract_values(self.dr.free_levels(mom,irr,label_number)[0])
            #     mb, m = extract_values(self.dr.free_levels(mom,irr,label_number)[1])
            #     ma = mass_map[ma]
            #     mb = mass_map[mb]
            #     #data_array = self.energy_data(mom,irr,label)
            #     data_list.append(deltaE(ecm_data[total_label_num],ma,mb,n,m,int(mom[3:])))
            
            data = np.array(data_list)
            return data

        # set up covaraince matrixes for each channel
        self.covariance_matrix = {}
        self.cov_de = {}
        for channel in self.channel:
            self.covariance_matrix[channel] = np.cov(ecm_NN_bs_arr[channel])
            self.cov_de[channel] = np.cov(energy_shift_data(channel))
        

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


        # print(f" Using Effective Range Expansion for single-channel:{self.channel_1} (m = {self.ma_ave}) ,{self.channel_2} (m = {self.mb_ave})")


        # next lets run a fit to check
        #logging.info(r" Minimizing ERE for average data set")
        def average_fit(): #~~~ returning the mean bootstrap sample fit 
            result = minimize(chi2,x0=[0.04,0.6],method='nelder-mead')
            print(result)
            return [result.x[0],result.x[1]]
    
            
        def deriv(n,psq,irrep,x): # 2 parameter difference derivative
            if len(x) == 1:
                a = x
                b = 0
            else:
                a, b = x

            eps = .001
            if n == 0:
                return (QC1(psq,irrep,a-eps,b)-QC1(psq,irrep,a+eps,b))/(2*eps)
            elif n ==1:
                return (QC1(psq,irrep,a,b-eps)-QC1(psq,irrep,a,b+eps))/(2*eps)
        
        def vij(x): #V_ij is error matrix for parameters, take _ii to get each error
            # v_nm = dp_i / dp_n
            # error estimation function
            # cov = self.cov = self.data.covariance_data()
            if len(x) == 2:
                nint = [0,1]
            else: 
                nint = 0
            # irreps = {'PSQ0': ['G1u']
            psq = ['PSQ0','PSQ1','PSQ2','PSQ3']
            lmat = np.empty(0)
            for n in nint:
                if n == 0:
                    dl = []
                    for i in psq:
                        dl.append(deriv(n,i,self.irreps[i][0],x))

                    lmat = np.append(lmat,dl)
                else:
                    dl =[]
                    for i in psq:
                        dl.append(deriv(n,i,self.irreps[i][0],x))
                    
                    lmat = np.vstack([lmat,dl])
                
            Vnm = np.linalg.inv(lmat@np.linalg.inv(self.cov_de)@np.transpose(lmat))
            
            return Vnm

        def find_intersection(x, y1, y2):
            x_range = np.linspace(-0.10,-.05,100) # x range is q2 range
            f1 = np.interp(x_range,x,y1)
            f2 = np.interp(x_range,x,y2)
            # Find the indices where the two curves intersect (assuming they do)
            root = []
            for i in range(len(f1)):
                if f1[i] - f2[i] < .001:
                    root.append(x_range[i])
            return root

        def find_bound_state(x):
            q2_values = np.linspace(-0.15, -0.001, 300) # need to change q2 based on data (future)
            virtual_state = []
            for q2 in q2_values:
                virtual_state.append(cmath.sqrt(-q2))
            # x is best fit parameters
            a = x[0]
            b = x[1]
            best_fit_line = []
            for q2 in q2_values:
                ecm = np.sqrt(q2 +self.ma_ave**2 ) + np.sqrt(q2 +self.mb_ave**2 )
                best_fit_line.append(parametrizations.ere_delta(ecm,self.ma_ave,self.mb_ave,a,b))

            bound_mom_2 =  find_intersection(q2_values,virtual_state,best_fit_line)[0]

            # self.vnm_matrix = vij([a,b])
            # derivative errors 
            # g^T V g -> g is d(parametrrization)/d(param), so d(ERE_delta)/d(param)
            vec = lambda ecm:np.array([ecm,ecm*parametrizations.delta_Sp(ecm,self.ma_ave,self.mb_ave)]) 

            ecm_bound = np.sqrt(bound_mom_2 + self.ma_ave**2 ) + np.sqrt(bound_mom_2 +self.mb_ave**2 )
            sigma_f = np.sqrt(np.transpose(vec(ecm_bound))@self.vnm_matrix@vec(ecm_bound)) 
            # error is quadrature error
            def pEpk(q2):
                return np.sqrt(-q2)*(-q2 + self.ma_ave**2 )**(-1/2) + np.sqrt(-q2)*(-q2 + self.mb_ave**2 )**(-1/2)
            bound_error = np.sqrt((pEpk(bound_mom_2) * sigma_f)**2)
            return [ecm_bound,bound_error] #need to fix to include errors in Sigma data

        self.best_fit = average_fit()
        self.vnm_matrix = vij(self.best_fit)
        self.errors_in_parameters = np.diag(self.vnm_matrix)
        self.bound_state = find_bound_state(self.best_fit)
        return [self.best_fit,self.bound_state]#print(self.best_fit)
        # do the task, produce the data, data goes in self.proj_dir_handler.data_dir(), info/warning/errors about the process goes in self.proj_dir_handler.log_dir() (if any)

    def plot( self ):
        # data is here
        # self.ecm_average_data  for each irrep and psq
        psq_list = self.dr.load_psq()
        irreps = self.irreps 
        shapes = ['o','s','D','v'] #automate to make shapes for different psq
        labels = ['$G_{1u}(0)$','$G_1 (1)$', '$G (2)$' , '$G(3)$']
        legend_handles = []  # Create an empty list to store custom legend handles
        shapes_dict = {}
        labels_dict = {}
        enu = 0
        for psq in psq_list: 
            shapes_dict[psq] = shapes[enu]
            labels_dict[psq] = labels[enu]
            enu += 1
        #irrep labels
        # q2 for values near bound state condition
        q2_values = np.linspace(-0.25, 0.01, 150)
        virtual_state = []
        for q2 in q2_values:
            virtual_state.append(cmath.sqrt(-q2))

        plt.figure(figsize=(10,6.134))
        for psq in psq_list:
            for irrep in irreps[psq]:
                x = kinematics.q2(self.ecm_average_data[psq][irrep], self.ma_ave,self.mb_ave)
                y = kinematics.qcotd(self.ecm_average_data[psq][irrep],self.L,psq,self.ma_ave,self.mb_ave,self.ref_ave)
                plt.plot(x, y, marker=shapes_dict[psq], color='blue', label=labels_dict[psq])
                legend_handles.append(Line2D([0], [0], marker=shapes_dict[psq], color='w', markerfacecolor='blue', markersize=10, label=labels_dict[psq]))
        
        for psq in psq_list:
            x_range = []
            y_range = []
            for irrep in irreps[psq]:
                std_deviation = np.std(self.ecm_bs[psq][irrep])
                for en in np.linspace(self.ecm_average_data[psq][irrep]-std_deviation , self.ecm_average_data[psq][irrep] +std_deviation , 100):
                    x_range.append(kinematics.q2(en, self.ma_ave,self.mb_ave))
                    y_range.append(kinematics.qcotd(en,self.L,psq,self.ma_ave,self.mb_ave,self.ref_ave))
            plt.plot(x_range, y_range, color="blue", alpha=0.8)  
        
        a,b = self.best_fit
        best_fit_line = []
        for q2 in q2_values:
            ecm = np.sqrt( q2 + self.ma_ave**2) + np.sqrt( q2 + self.mb_ave**2 )
            best_fit_line.append(parametrizations.ere_delta(ecm,self.ma_ave, self.mb_ave,a,b))

        vij = self.vnm_matrix # using best_fit
        # vec for error estimation is derivative in each parameter of paramerization (ERE_Delta)
        vec = lambda ecm:np.array([ecm,ecm*parametrizations.delta_Sp(ecm,self.ma_ave, self.mb_ave)]) #np.array([deriv(ecm,a,b,0,.001),deriv(ecm,a,b,1,.001)]) #np.array([ecm,ecm*self.fit.ere_delta(ecm,0,a,b)])
        sigma_f = [np.sqrt(np.transpose(vec(kinematics.q2toecm(q2,self.ma_ave, self.mb_ave)))@vij@vec(kinematics.q2toecm(q2,self.ma_ave, self.mb_ave))) for q2 in q2_values]
        upper = np.array(best_fit_line) + np.array(sigma_f)
        lower = np.array(best_fit_line) - np.array(sigma_f)
        #fk_values = [self.qc.q2(ecm,0) for ecm in ecm_values]
        #bound_mom = self.bound_state[0]   
        plt.plot(q2_values,best_fit_line , color='blue',linestyle='-.')
        #plt.fill_between(fk_values, lower_bounds, upper_bounds, color='lightblue', alpha=0.5)
        #plt.plot(bound_mom,parametrizations.ere_delta(kinematics.q2toecm(bound_mom,self.ma_ave, self.mb_ave),self.ma_ave, self.mb_ave, self.best_fit[0], self.best_fit[1]),'r*',markersize=10)
        plt.fill_between(q2_values,lower,upper,alpha = 0.5, color = 'lightblue')
        plt.axhline(y=0,color='black')
        plt.axvline(x=0,color='black')
        legend = plt.legend(handles=legend_handles, loc='upper left', title='Legend', prop={'size': 12})
        plt.xlabel("$q^{*2} / m_{\pi}^2$",fontsize=16)
        plt.ylabel("$q^{*} / m_{\pi} \cot \delta $",fontsize=16)
        plt.title(f'{self.channel_1},{self.channel_2}  Scattering ',fontsize=16)
        plt.savefig(f'{self.channel_1},{self.channel_2} _Scattering.pdf')
        #legend.set_title('Legend', prop={'size': 12})  # Set legend title and font size

        plt.show()
        return print(f"Plotting for {self.channel_1},{self.channel_2} Complete")
        
        # make the plots, store in self.proj_dir_handler.plot_dir(), again, any log/error warnings go in self.proj_dir_handler.log_dir() as well (if any)
def get_particle_name(particle_str):
    return particle_str.split("(")[0]