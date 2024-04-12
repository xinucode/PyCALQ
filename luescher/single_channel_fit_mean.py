import logging
import general.data_reader as dr
import general.plotting_handler as ph
import fvspectrum.sigmond_util as sigmond_util

############################################################################################
##      Required packages
import cmath                    #math library (complex)
import csv                      # read csv files
import json
import math                     # math library, mostly use it for math.pi, math.pow()
import matplotlib.pyplot as plt #plotting library            
import numpy as np              #basic functions, linear algebra, etc. 
import random
import scipy.special as sp      #special functions
import scipy.integrate as integrate #needed to do numerical integration, alt?
import yaml
import os
##
from matplotlib.lines import Line2D  # Import Line2D for custom legend handles
#from tqdm import tqdm  # Import the tqdm library
from scipy import optimize,integrate
from scipy.optimize import fsolve,minimize
from scipy.integrate import quad
from luescher.zeta import Z
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
        self.proj_handler = proj_handler 
        self.task_name = task_name
        # if not task_configs:
        #     logging.critical(f"No directory to view. Add 'raw_data_files' to '{task_name}' task parameters.")

        ensemble_info = sigmond_util.get_ensemble_info(general_configs)
        self.L = ensemble_info.getLatticeXExtent()
        #print(task_params)
        #print(self.L)
        # path name needs to link toward the hdf5 file 
        if 'data_file' in task_params.keys():
            file_path = task_params['data_file'][0]
            #print(file_path)

        self.dr = dr.LQCD_DATA_READER(file_path,self.L) # data reader is object dr
        self.data = self.dr.load_data()

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

        #hadron list
        self.single_hadron_list = np.array(self.dr.single_hadron_list())
        # irreps used in fit
        self.irreps = task_params['irreps']
        print(self.irreps)

        self.channel = task_params['channel']
        self.channel_1 = self.channel[0]
        self.channel_2 = self.channel[1]
        #print(self.single_hadron_list)
        # also check the ref mass
        self.ref_mass = self.dr.single_hadron_data('Pion')[0]
        #print(self.ref_mass)

        if np.any(np.isin(self.single_hadron_list,self.channel_1)) and np.any(np.isin(self.single_hadron_list,self.channel_2)):
            logging.info("Channel is confirmed to be in data file. Continuing analysis ...")
        else:
            logging.critical("Channels not in data file")

        
        #print(task_params['irreps']['PSQ0'])
        #initialize your task, store default input in self.proj_dir_handler.log_dir() (basically, throw the full possible input with all parameters where all the assumed parameters have been filled in in there)
        with open( os.path.join(proj_handler.log_dir(), 'full_input.yml'), 'w+') as log_file:
            yaml.dump({"general":general_configs}, log_file)
            yaml.dump({"tasks":task_params}, log_file)
            print("Initialization complete.")


    def momentum_state(self,i):
        # d = [0,0,0]
        if i == 0:
            return np.array([i,i,i])
        # d = [0,0,1]
        elif i == 1:
            return np.array([0,0,i])
        # d = [1,1,0]
        elif i == 2:
            return np.array([i-1,i-1,0])
        # d = [1,1,1]
        elif i == 3:
            return np.array([i/3,i/3,i/3])
        # d = [0,0,2]
        elif i == 4:
            return np.array([0,0,2])
        else: 
            # Raise an exception for invalid input
            raise ValueError("Invalid value for 'i'. 'i' must be 0, 1, 2, or 3.")

    def momentum_state_string_input(self,i):
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
    
    def q2(self,ecm,ma,mb):
        #ecm is the energy data
        #ma,mb comes from the energies in the channel
        q2 = ecm**2 / 4 - (ma**2 + mb**2) / 2 + ((ma**2 - mb**2)**2) / (4*ecm**2)
        return q2
    
    def msplit(self,ecm,ma,mb): #if ma,mb are degenerate its 1
        return 1 + ((ma**2 - mb**2)/ecm**2 ) 
    # delta is distance to threshold in energy, for fitting
    def delta_ab(self,ecm,ma,mb):
        #mpi, mS = self.data.sigma_pi_masses()
        return (ecm**2 - (ma+mb)**2 ) / (ma+ mb)**2
    
     #gamma, lorentz factor
    def gamma(self,ecm,d):
        d_vec = self.momentum_state_string_input(d) 
        L_ref = self.L*self.ref_mass
        E = math.sqrt(ecm**2 + (((2*math.pi)/L_ref)**2)*np.dot(d_vec,d_vec))
        return E/ecm
    
    # s-wave luscher quantization condition
    def qcotd(self,ecm,psq,ma,mb):
        L_ref = self.L*self.ref_mass
        d_vec = self.momentum_state_string_input(psq) #0,1,2,3
        c = 2 / (self.gamma(ecm,psq)*L_ref*math.sqrt(math.pi))
        #print(f'g={self.gamma(ecm,psq)},L = {L_ref},q2={self.q2(ecm,ma,mb)},ms={self.msplit(ecm,ma,mb)}')
        return c*Z(self.q2(ecm,ma,mb)*((L_ref/(2*math.pi))**2),gamma=self.gamma(ecm,psq),l=0,m=0,d=d_vec,m_split=self.msplit(ecm,ma,mb),precision=1e-11).real


    def run( self ):       
        #log_path = os.path.join(self.proj_handler.log_dir(), 'luescher_log.yml') 
        # step 1, import the keys needed for analysis from single_hadron_list
        # first get the required channel masses
        ma = self.dr.single_hadron_data(self.channel[0]) # make sure using ref masses
        mb = self.dr.single_hadron_data(self.channel[1]) # 
        # next load in the data
        # lowest level from each data
        # if psq != 0, do next level instead
        psq_list = self.dr.load_psq()
        #print(psq_list)
        # collect irreps for each of the PSQ, which is the available data
        # irreps = {}
        # for psq in psq_list:
        #     irreps[psq] = []
        #     for key in self.data.get(psq):
        #         irreps[psq].append(key)
        #print(self.irreps)
        # this section gets our data together to do covariance matrix
        # data is obtained for specific irreps
        self.ecm_data = {} # save possible data
        self.ecm_average_data = {} #saving only the average data
        self.ecm_bs = {}
        ecm_bs_arr = []
        for psq in psq_list:
            self.ecm_data[psq] = {}
            self.ecm_average_data[psq] = {}
            self.ecm_bs[psq] = {}
            irrep = self.irreps[psq]
            #print(irrep)
            if isinstance(irrep, list):
                for irrs in irrep:
                        if psq == 'PSQ0':
                            level = "ecm_0_ref" # level as an int
                        else:
                            level = "ecm_1_ref"
                        self.ecm_data[psq][irrs] = self.data.get(psq).get(irrs).get(level)[:]
                        self.ecm_average_data[psq][irrep] = self.ecm_data[psq][irrep][0]
                        self.ecm_bs[psq][irrep] = self.ecm_data[psq][irrep][1:]
                        ecm_bs_arr.append(self.ecm_bs[psq][irrep])
            else:
                if psq == 'PSQ0':
                    level = "ecm_0_ref" # level as an int
                else:
                    level = "ecm_1_ref"
                self.ecm_data[psq][irrep] = self.data.get(psq).get(irrep).get(level)[:]
                self.ecm_average_data[psq][irrep] = self.ecm_data[psq][irrep][0]
                self.ecm_bs[psq][irrep] = self.ecm_data[psq][irrep][1:]
                ecm_bs_arr.append(self.ecm_bs[psq][irrep])
  
        self.covariance_matrix = np.cov(np.array(ecm_bs_arr))
        #print(self.covariance_matrix)
        # now that we have ecm_data, lets do a fit and save results
        #self.average_energies = []
        #first we want to use average data
        # self.ecm_average_data = {}
        # for psq in psq_list:
        #     self.ecm_average_data[psq] = {}
        #     for irrep in self.irreps[psq]:
        #         self.ecm_average_data[psq][irrep] = self.ecm_data[psq][irrep][0]
        #         #self.average_energies.append(self.ecm_data[psq][irrep][0])
        

        with open( os.path.join(self.proj_handler.log_dir(), 'energies.json') , 'w+') as log_file:
            json.dump(self.ecm_average_data, log_file)

        # with open( os.path.join(self.proj_handler.log_dir(), 'cov.json') , 'w+') as log_file:
        #     json.dump(self.covariance_matrix, log_file)
        




        self.ma_ave = ma[0]
        self.mb_ave = mb[0]
        #print(self.ma_ave)

        # # get cov data
        # self.ecm_bs = {}
        # for psq in psq_list:
        #     self.ecm_bs[psq] = {}
        #     for irrep in self.irreps[psq]:
        #         self.ecm_bs[psq][irrep] = self.ecm_data[psq][irrep][1:].tolist()

        # ecm_NN_bs_arr = []
        # for psq in psq_list:
        #     for irrep in self.irreps[psq]:
        #         ecm_NN_bs_arr.append(self.ecm_bs[psq][irrep]) # the ordering should be kep tin which the calculation will go

        # with open( os.path.join(self.proj_handler.log_dir(), 'energies_bs.json') , 'w+') as log_file:
        #     json.dump(self.ecm_bs, log_file)



        # fit ere expansion
        def ere(ecm,a,b):
            return ((-1/a)+0.5*b*self.q2(ecm,self.ma_ave,self.mb_ave)) #in units of reference mass, usually mpi

        def deter(ecm,psq,a,b):
            #p = psq[3]
            return self.qcotd(ecm,psq,self.ma_ave,self.mb_ave) - ere(ecm,a,b)

        def QC1(psq,irrep, a, b):
            func = lambda ecm: deter(ecm,psq, a, b)
            return fsolve(func,self.ecm_average_data[psq][irrep])[0]

        #print('ere=',ere(6.7,1,1))
        #print('deter=',deter(6.7,'PSQ0',1,1))
        #print('QC1=',QC1('PSQ0','G1u',1,1))

        # def chi2(x):
        #     res = []
        #     for psq in psq_list:
        #         irrep = self.irreps[psq]
        #         for irrep in self.irreps[psq]:
        #             #d = psq[3]
        #             res.append(self.ecm_average_data[psq][irrep] - QC1(psq,irrep,a,b))
        #     value = np.array(res)@np.linalg.inv(self.covariance_matrix)@np.array(res)
        #     return value

        def chi2(x,b=0):
            if isinstance(x, (int,float)):
                a = x
            else:
                a,b = x

            res = []
            for psq in psq_list:
                irreps = self.irreps[psq]
                if isinstance(irreps, list):
                    for irrep in irreps:
                        #d = psq[3]
                        res.append(self.ecm_average_data[psq][irrep] - QC1(psq,irrep,a,b))
                else:
                    res.append(self.ecm_average_data[psq][irreps] - QC1(psq,irreps,a,b))

            value = np.array(res)@np.linalg.inv(self.covariance_matrix)@np.array(res)
            return value
        
        #print(chi2(1,1))
        logging.info(f" Using Effective Range Expansion for single-channel:{self.channel_1} (m = {self.ma_ave} $m_pi$) ,{self.channel_2} (m = {self.mb_ave} $m_pi$)")
        print(f" Using Effective Range Expansion for single-channel:{self.channel_1} (m = {self.ma_ave} $m_pi$) ,{self.channel_2} (m = {self.mb_ave} $m_pi$)")
        # next lets run a fit to check
        #logging.info(r" Minimizing ERE for average data set")
        def average_fit(output_file=os.path.join(self.proj_handler.log_dir(), f'fit_results_channel_{self.channel_1}{self.channel_2}.txt')):
            # will do a 1 and 2 parameter fit and see which is best
            # best_a = None
            # best_b = None
            # best_chi2 = 
            try:
                with open(output_file, 'w') as file:
                    print("Fit starting:")
                    for fit_number in range(2):
                        if fit_number == 0: # fitting just with a
                            print("1 parameter fit")
                            retry = True
                            count = 0
                            bigger_count = 0
                            while retry:
                                try:
                                    a_guess = np.random.uniform(-10,10)
                                    print(f"trying a={a_guess}")
                                    result = minimize(chi2, x0 = a_guess , args = 0 , method = 'nelder-mead')
                                    if result.success:
                                        print(result.fun)
                                        best_a_0 = a_guess
                                        chi2_temp = result.fun
                                        if count == 0:
                                            chi2_res_1 = chi2_temp
                                            count += 1

                                        else:
                                            if chi2_temp - chi2_res_1 < 1e-2: # if chi2 not changing with different guess, its robust
                                                retry = False
                                                #break

                                            elif chi2_temp < chi2_res_1: #if smaller chi2, try again
                                                chi2_res_1 = chi2_temp
                                                count += 1

                                            elif chi2_temp > chi2_res_1: #if bigger chi2, try again
                                                bigger_count += 1 
                                                if bigger_count > 3:
                                                    retry = False
                                                    #break
                                                #chi2_res_1 = chi2_temp
                                    else:
                                        print("retrying")
                                        # retry with new parameters
                                except Exception as e:
                                    logging.error("Error with fit.")

                            logging.info(f"Minimization for 1-parameter ERE: a = {best_a_0}, chi2 = {chi2_res_1} ")
                            print(f"Minimization for 1-parameter ERE: a = {best_a_0}, chi2 = {chi2_res_1} ")
                            file.write(f"{fit_number+1} {best_a_0} {chi2_res_1}\n")

                        elif fit_number == 1: # 2 parameter ERE
                            retry = True
                            count = 0
                            bigger_count = 0
                            while retry:
                                try:
                                    if count == 0:
                                        a_guess = best_a_0
                                    else:
                                        a_guess = np.random.uniform(-20,20)

                                    b_guess = np.random.uniform(-5,5)
                                    result = minimize(chi2, x0 = [a_guess , b_guess], method = 'nelder-mead')
                                    if result.success:
                                        best_a = a_guess
                                        best_b = b_guess
                                        chi2_temp = result.fun
                                        if count == 0:
                                            chi2_res_2 = chi2_temp
                                            count += 1
                                        else:
                                            if chi2_temp - chi2_res_2 < 1e-3: # if chi2 not changing with different guess, its robust
                                                retry = False

                                            elif chi2_temp < chi2_res_2: #if smaller chi2, try again
                                                chi2_res_2 = chi2_temp
                                                count += 1
                                            
                                            elif chi2_temp > chi2_res_2: #if bigger chi2, try again
                                                bigger_count += 1 
                                                if bigger_count > 3:
                                                    retry = False
                                    else:
                                        pass
                                        # retry with new parameters
                                except Exception as e:
                                    print("Error with fit, trying again, tried:")
                                    logging.error("Error with fit.")

                            logging.info(f"Minimization for 2-parameter ERE: a = {best_a}, b = {best_b}, chi2 = {chi2_res_2} ")
                            file.write(f"{fit_number+1} {best_a} {best_b} {chi2_res_2}\n")

            except Exception as e:
                logging.error(f"Error writing to file:{e}")

            logging.info(f" ERE fits written to file: {output_file}")                            
            if chi2_res_1 < chi2_res_2:
                return best_a_0
            else: 
                return [best_a, best_b]
            

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
        
        def vij(x): #V_ij is error matrix for parameters, take _ii to get each error
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

        #self.best_fit = average_fit()
        return print("fit is done")
        # do the task, produce the data, data goes in self.proj_dir_handler.data_dir(), info/warning/errors about the process goes in self.proj_dir_handler.log_dir() (if any)


    def plot( self ):
        # x = []
        # x_range = []
        # y = []
        # y_range = []
        # average_data = []
        # for keys in self.ecm_average_data:
        #     average_data.append(self.ecm_average_data[keys])

        # print(average_data)

        # for i in range(len(average_data)):
        #     x.append(self.q2(average_data[i], self.ma_ave, self.mb_ave))
        #     y.append(self.qcotd(average_data[i], i, self.ma_ave, self.mb_ave))
        #     if i == 0:
        #         x_range_i = []
        #         y_range_i = []
        #         for en in np.linspace(average_data[i] - 0.01, average_data[i] + 0.01, 100):
        #             x_range_i.append(self.q2(en, self.ma_ave, self.mb_ave))
        #             y_range_i.append(self.qcotd(en, i, self.ma_ave, self.mb_ave))
        #         x_range.append(x_range_i)
        #         y_range.append(y_range_i)
        #     else:
        #         xp = []
        #         yp = []
        #         for en in np.linspace(average_data[i] - 0.01, average_data[i] + 0.01, 100):
        #             xp.append(self.q2(en,  self.ma_ave, self.mb_ave))
        #             yp.append(self.qcotd(en, i,  self.ma_ave, self.mb_ave))
        #         x_range.append(xp)
        #         y_range.append(yp)

        # ecm_values = np.linspace(min(average_data)-2,max(average_data)+2,1000)
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
        # for i in range(len(average_data)):
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
        
        return 
        # make the plots, store in self.proj_dir_handler.plot_dir(), again, any log/error warnings go in self.proj_dir_handler.log_dir() as well (if any)