import logging
import general.data_reader as dr
import general.plotting_handler as ph
import fvspectrum.sigmond_util as sigmond_util

############################################################################################
##      Required packages
import cmath                    #math library (complex)
import csv                      # read csv files
import math                     # math library, mostly use it for math.pi, math.pow()
import datetime
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

        self.ensemble_info = sigmond_util.get_ensemble_info(general_configs)
        self.L = self.ensemble_info.getLatticeXExtent()
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
            'figwidth':10,
            'figheight':6.132,
            'ref_energies' : True, #trigger for using reference energies in analysis, currently default and support is TRUE
            'delta_E_covariance':True,
            'error_estimation': False,
            'chi2_energy_compare': True 
        }
        self.channels_and_irreps = task_params['channels']#self.alt_params['irreps']

        #hadron list
        self.single_hadron_list = np.array(self.dr.single_hadron_list())
        print(self.single_hadron_list)
        # generate list of channels
        self.channel = []
        self.irreps = {}
        self.fit_parametrization = {}
        for channel in self.channels_and_irreps:
            self.channel.append(channel)
            self.irreps[channel] = self.channels_and_irreps[channel]
            if self.alt_params.get('parametrization') or task_params.get('parametrization'):
                self.fit_parametrization[channel] = task_params['parametrization']
            else:
                # automatic param is delta ERE for now
                self.fit_parametrization[channel] = 'ERE_delta'
                logging.info('Parametrization not chosen. Default is Effective Range Expansion (ERE)')
        print('channels',self.channel)
        print("irreps",self.irreps)
        for channel in self.channel:
            #print(channel)
            channel_1 = str(channel.split(',')[0])
            channel_2 = str(channel.split(',')[1])
            if np.any(np.isin(self.single_hadron_list,channel_1)) and np.any(np.isin(self.single_hadron_list,channel_2)):
                logging.info(f"scattering Channel {channel} is confirmed to be in data file. Continuing analysis ...")
            else:
                logging.critical(f"Scattering Channel {channel} not found in data. Must be hadrons including {self.single_hadron_list}")

        # check parametrization
        
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
         
        # step 1, import the keys needed for analysis from single_hadron_list
        # first get the required channel masses
        self.log_path = {}
        self.m_ref_dict = {}
        self.ref_mass = {}
        self.psq_list = {}
        self.irrep_list = {}
        for i, channel in enumerate(self.channel):
            channel_1 = str(channel.split(',')[0])
            channel_2 = str(channel.split(',')[1])
            self.m_ref_dict[channel] = [self.dr.single_hadron_data(channel_1),self.dr.single_hadron_data(channel_2) ] #ma_ref, mb_ref in channel
            self.ref_mass[channel] = self.dr.single_hadron_data('ref')
            self.psq_list[channel] = self.dr.load_psq()
            self.log_path[channel] = os.path.join(self.proj_handler.log_dir(), f'luescher_{channel_1}_{channel_2}_log.txt')

        self.ecm_data = {} # save possible data
        self.ecm_average_data = {}
        self.ecm_bootstrap_data = {}
        ecm_bs_arr = {}
        for channel in self.channel:
            self.ecm_data[channel] = {} 
            self.ecm_average_data[channel] = {}
            self.ecm_bootstrap_data[channel] = {}
            ecm_bs_arr[channel] = []
            self.irrep_list[channel] = {}
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
                            ecm_bs_arr[channel].append(self.ecm_bootstrap_data[channel][psq][irrep][level_title])
                        else:
                            level_title = f"ecm_{level}"
                            logging.critical("Need non-ref energies from data reader")

        def energy_shift_data(channel):
                
            ecm_data = self.ecm_bootstrap_data[channel]

            #Sarah
            mref = np.array(self.ref_mass[channel])[1:] #np.array(self.dr.single_hadron_data('ref'))[1:]
            # m1_ref,m2_ref = self.m_ref_dict[channel]
            # m1_ref = m1_ref[1:]
            # m2_ref = m2_ref[1:]
            # m1_ref = np.array(self.dr.single_hadron_data(self.channel_1))
            # m2_ref = np.array(self.dr.single_hadron_data(self.channel_2))
            channel_1 = str(channel.split(',')[0])
            channel_2 = str(channel.split(',')[1])
           # mapping for masses
        #    mass_map = {}
        #    for 
        #     mass_map = {
        #         get_particle_name(channel_1): m1_ref[1:],
        #         get_particle_name(channel_2): m2_ref[1:],
        #         #get_particle_name(channel_2): m2_ref[1:],
        #     }
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
                        if self.alt_params['ref_energies']:
                            hadron_title_a = f'{ma}(0)_ref'
                            hadron_title_b = f'{mb}(0)_ref'
                        else:
                            hadron_title_a = f'{ma}(0)'
                            hadron_title_b = f'{mb}(0)'
                        ma = self.dr.single_hadron_data(hadron_title_a)[1:]
                        mb = self.dr.single_hadron_data(hadron_title_b)[1:]

                        data_list.append(deltaE(ecm_data[psq][irrep][level_title] ,ma,mb,n,m,int(psq[3:])))

            data = np.array(data_list)
            return data

        # set up covaraince matrixes for each channel
        self.covariance_matrix = {}
        for channel in self.channel:
            if self.alt_params['delta_E_covariance']:
                self.covariance_matrix[channel] =  np.cov(energy_shift_data(channel))
            else:
                self.covariance_matrix[channel] = np.cov(ecm_bs_arr[channel])

        def determinant_condition(ecm,psq,ma,mb,ref,fit_parametrization, fit_params):
            # p is priors, a, b , ... for fits
            #p = psq[3]
            # parametrizations.ere_delta(ecm,self.ma_ave,self.mb_ave,a,b)
            return (kinematics.qcotd(ecm,self.L,psq,ma,mb,ref) - parametrizations.output(ecm,ma,mb,fit_parametrization, fit_params) )
        
        def QC1(energy,psq,ma,mb,ref,fit_parametrization, fit_params):
            #self.ecm_average_data[psq][irrep]
            func = lambda ecm: determinant_condition(ecm,psq,ma,mb,ref,fit_parametrization, fit_params)
            # energy is the expected energy level
            return fsolve(func,energy)[0] #guess is the energy going in

        def chi2(fit_params,channel):
            res = []
            for psq in self.psq_list[channel]:
                for irrep in self.irreps[channel][psq][0]:
                    for level in self.irreps[channel][psq][0][irrep]:
                        level_title = f"ecm_{level}_ref"
                        ma, mb = self.m_ref_dict[channel]  
                        ma = ma[0]
                        mb = mb[0]
                        ref = self.ref_mass[channel][0]            
                        diff = self.ecm_average_data[channel][psq][irrep][level_title] - QC1(self.ecm_average_data[channel][psq][irrep][level_title],psq,ma,mb,ref,self.fit_parametrization[channel],fit_params)
                        res.append(diff)
            value = np.array(res)@np.linalg.inv(self.covariance_matrix[channel])@np.array(res)
            return value

        def average_fit(channel): #~~~ returning the mean bootstrap sample fit 
            result = minimize(chi2,x0=[0.01,0.7],args=(channel),method='nelder-mead')
            #print(result)
            return result#[result.x[0],result.x[1]]
    
        def deriv(n,energy,psq,ma,mb,ref, fit_param,fit_params):  # Generalized derivative function
            # first order difference, need to add more orders 
            eps = 0.001  # Small perturbation
            # QC1(energy,psq,ma,mb,ref,fit_params)
            x_eps = fit_params.copy()  # Create a copy of x to perturb
            # Apply the perturbation to the nth parameter
            x_eps[n] -= eps
            QC1_minus = QC1(energy,psq,ma,mb,ref,fit_param, x_eps)

            x_eps = fit_params.copy()
            x_eps[n] +=  eps  # Reset and apply in the other direction
            QC1_plus = QC1(energy,psq,ma,mb,ref,fit_param, x_eps)

            # Calculate the derivative
            return (QC1_minus - QC1_plus) / (2 * eps)   

        
        def vij(channel, fit_params): #V_ij is error matrix for parameters, take _ii to get each error
            # v_nm = dp_i / dp_n
            # error estimation function
            # Determine the number of parameters
            num_params = len(fit_params)
            nint = list(range(num_params))  # Generate a list [0, 1, 2, ..., num_params-1]
            # deriv(n,energy,psq,ma,mb,ref, fit_params)
            lmat = []
            ma, mb = self.m_ref_dict[channel]
            ref = self.ref_mass[channel]
            for n in nint:
                dl = []
                for psq in self.psq_list[channel]:
                    for irrep in self.irreps[channel][psq][0]:
                        for level in self.irreps[channel][psq][0][irrep]:
                            level_title = f"ecm_{level}_ref"
                            dl.append(deriv(n, self.ecm_average_data[channel][psq][irrep][level_title], psq, ma[0],mb[0],ref[0], self.fit_parametrization[channel] ,fit_params))
                
                lmat.append(np.array(dl))
            lmat = np.array(lmat)
    
            Vnm = np.linalg.inv(lmat@np.linalg.inv(self.covariance_matrix[channel])@np.transpose(lmat))
            
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
        
        self.fit_results = {}
        if self.alt_params['error_estimation']:
            self.vnm_matrix = {}
        for channel in self.channel:
            logging.info(f"Fit results in {self.log_path[channel]}")
            average_fit_results = average_fit(channel)
            self.fit_results[channel] = list(average_fit_results.x)
            if self.alt_params['error_estimation']:
                self.vnm_matrix[channel] = vij(channel,self.fit_results[channel])
            current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            with open( self.log_path[channel], 'w+') as log_file:
                log_file.write(f"Log date and time: {current_time}\n")
                log_file.write(f"Ensemble: {self.ensemble_info}\n")
                log_file.write(f"Fit results for Scattering channel: {channel}\n")
                log_file.write(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
                log_file.write(f"Irreps analyzed: {self.irreps[channel]} \n")
                log_file.write(f"Average data: {self.ecm_average_data[channel]} \n")
                log_file.write(f"Parametrization used: {self.fit_parametrization[channel]}\n")
                log_file.write(f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
                log_file.write(f"{average_fit_results }\n")
                log_file.write(f"\n")
                log_file.write(f"Covariance Matrix: {self.covariance_matrix[channel]} \n")
                if self.alt_params['error_estimation']:
                    log_file.write(f"V_nm (estimated uncertainty in parameters): {self.vnm_matrix[channel]} \n")

        #self.best_fit = average_fit()
        # self.vnm_matrix = vij(self.best_fit)
        # self.errors_in_parameters = np.diag(self.vnm_matrix)
        # self.bound_state = find_bound_state(self.best_fit)
        #[self.best_fit,self.bound_state]#print(self.best_fit)
        # do the task, produce the data, data goes in self.proj_dir_handler.data_dir(), info/warning/errors about the process goes in self.proj_dir_handler.log_dir() (if any)

    def plot( self ):
        if self.alt_params['plot']:
            logging.info(f"Saving plots to directory {self.proj_handler.log_dir()}...")
        else:
            logging.info(f"No plots requested.")
            return
        fig_params = [self.alt_params['figwidth'],self.alt_params['figheight']]
        #log_path = self.proj_handler.log_dir()
        # data is here
        # self.ecm_average_data  for each irrep and psq
        
        def determinant_condition(ecm,psq,ma,mb,ref,fit_parametrization, fit_params):
            # p is priors, a, b , ... for fits
            #p = psq[3]
            # parametrizations.ere_delta(ecm,self.ma_ave,self.mb_ave,a,b)
            return (kinematics.qcotd(ecm,self.L,psq,ma,mb,ref) - parametrizations.output(ecm,ma,mb,fit_parametrization, fit_params) )
        
        def QC1(energy,psq,ma,mb,ref,fit_parametrization, fit_params):
            #self.ecm_average_data[psq][irrep]
            func = lambda ecm: determinant_condition(ecm,psq,ma,mb,ref,fit_parametrization, fit_params)
            # energy is the expected energy level
            return fsolve(func,energy)[0] #guess is the energy going in


        for channel in self.channel:
            ma, mb = self.m_ref_dict[channel]  
            ma_ave = ma[0]
            mb_ave = mb[0]
            ref_ave = self.ref_mass[channel][0]
            psq_list = self.dr.load_psq()
            irreps = self.irreps[channel] 
            x = {}
            y = {}
            x_range = {}
            y_range = {}
            e_vals = []
            if self.alt_params['chi2_energy_compare']:
                chi2_energy = {}
            for psq in psq_list:
                x[psq] = {}
                y[psq] = {}
                x_range[psq] = {}
                y_range[psq] = {}
                if self.alt_params['chi2_energy_compare']:
                    chi2_energy[psq] = {} 
                for irrep in self.irreps[channel][psq][0]:
                    x[psq][irrep] = {}
                    y[psq][irrep] = {}
                    x_range[psq][irrep] = {}
                    y_range[psq][irrep] = {}
                    if self.alt_params['chi2_energy_compare']:
                        chi2_energy[psq][irrep] = {} 
                    for level in self.irreps[channel][psq][0][irrep]:
                        level_title = f"ecm_{level}_ref"
                        std_deviation = np.std(self.ecm_bootstrap_data[channel][psq][irrep][level_title])
                        x[psq][irrep][level] = kinematics.q2(self.ecm_average_data[channel][psq][irrep][level_title], ma_ave,mb_ave)
                        if self.alt_params['chi2_energy_compare']:
                            chi2_energy[psq][irrep][level] = QC1(self.ecm_average_data[channel][psq][irrep][level_title],psq,ma_ave,mb_ave,ref_ave,self.fit_parametrization[channel],self.fit_results[channel])
                        e_vals.append(self.ecm_average_data[channel][psq][irrep][level_title])
                        y[psq][irrep][level]  = kinematics.qcotd(self.ecm_average_data[channel][psq][irrep][level_title],self.L,psq,ma_ave,mb_ave,ref_ave)   
                        x_range[psq][irrep][level] = []
                        y_range[psq][irrep][level] = []
                        for en in np.linspace(self.ecm_average_data[channel][psq][irrep][level_title]-std_deviation , self.ecm_average_data[channel][psq][irrep][level_title] +std_deviation , 100):
                            x_range[psq][irrep][level].append(kinematics.q2(en, ma_ave,mb_ave))
                            y_range[psq][irrep][level].append(kinematics.qcotd(en,self.L,psq,ma_ave,mb_ave,ref_ave))
                        #plt.plot(x, y, marker=shapes_dict[psq], color='blue', label=labels_dict[psq])
                        # legend_handles.append(Line2D([0], [0], marker=shapes_dict[psq], color='w', markerfacecolor='blue', markersize=10, label=labels_dict[psq]))
            x_in = [x,x_range]
            y_in = [y,y_range]
            ph.PlottingHandler().single_channel_plot( fig_params, channel, irreps, x_in, y_in)
            ph.PlottingHandler().save_pdf(os.path.join(self.proj_handler.log_dir(), f'{channel}_Scattering_Data.pdf'), transparent=True)
            # find min and max of data in energy
            # ecm_min_value = min(
            #     self.ecm_average_data[channel][psq][irrep][level]
            #     for psq in psq_list
            #     for irrep in self.ecm_average_data[channel][psq]
            #     for level in self.ecm_average_data[channel][psq][irrep]
            # )
            # ecm_max_value = max(
            #     self.ecm_average_data[channel][psq][irrep][level]
            #     for psq in psq_list
            #     for irrep in self.ecm_average_data[channel][psq]
            #     for level in self.ecm_average_data[channel][psq][irrep]
            # )
            ecm_fit_values = np.linspace(min(e_vals)-0.08,max(e_vals)+0.08, 100)            
            # fit parametrization on top
            q2_for_fit = []
            best_fit_line = []
            for e in ecm_fit_values: # need to change to make q2 automatic
                q2_for_fit.append( kinematics.q2(e, ma_ave,mb_ave))
                # ecm = np.sqrt( q2 + ma_ave**2) + np.sqrt( q2 + mb_ave**2 )
                best_fit_line.append( parametrizations.output(e,ma_ave,mb_ave,self.fit_parametrization[channel],self.fit_results[channel]))
            plt.plot( q2_for_fit,best_fit_line, color='blue', lw=2,ls='--')
            ph.PlottingHandler().save_pdf(os.path.join(self.proj_handler.log_dir(), f'{channel}_Scattering_fit.pdf'), transparent=True)
            
            if self.alt_params['error_estimation']:
                pass

            if self.alt_params['chi2_energy_compare']:
                # add plotting with error here
                chi_in = [self.ecm_average_data,chi2_energy]
                ph.PlottingHandler().chi2_energies_compare_plot(fig_params, channel, irreps, chi_in)
                ph.PlottingHandler().save_pdf(os.path.join(self.proj_handler.log_dir(), f'{channel}_chi2_energy_compare.pdf'), transparent=True)
            else:
                pass
            
            
                # for now no error bars

            # vij = self.vnm_matrix[channel] # using best_fit
            # # vec for error estimation is derivative in each parameter of paramerization (ERE_Delta)
            # vec = lambda ecm:np.array([ecm,ecm*parametrizations.delta_Sp(ecm,ma_ave,mb_ave)]) #np.array([deriv(ecm,a,b,0,.001),deriv(ecm,a,b,1,.001)]) #np.array([ecm,ecm*self.fit.ere_delta(ecm,0,a,b)])
            # sigma_f = [np.sqrt(np.transpose(vec(kinematics.q2toecm(q2,ma_ave, mb_ave)))@vij@vec(kinematics.q2toecm(q2,ma_ave, mb_ave))) for q2 in q2_values]
            # upper = np.array(best_fit_line) + np.array(sigma_f)
            # lower = np.array(best_fit_line) - np.array(sigma_f)
            # #fk_values = [self.qc.q2(ecm,0) for ecm in ecm_values]
            # #bound_mom = self.bound_state[0]   
            # plt.plot(q2_values,best_fit_line , color='blue',linestyle='-.')
            # #plt.fill_between(fk_values, lower_bounds, upper_bounds, color='lightblue', alpha=0.5)
            # #plt.plot(bound_mom,parametrizations.ere_delta(kinematics.q2toecm(bound_mom,self.ma_ave, self.mb_ave),self.ma_ave, self.mb_ave, self.best_fit[0], self.best_fit[1]),'r*',markersize=10)
            # plt.fill_between(q2_values,lower,upper,alpha = 0.5, color = 'lightblue')
            # plt.axhline(y=0,color='black')
            # plt.axvline(x=0,color='black')
            # legend = plt.legend(handles=legend_handles, loc='upper left', title='Legend', prop={'size': 12})
            # plt.xlabel("$q^{*2} / m_{\pi}^2$",fontsize=16)
            # plt.ylabel("$q^{*} / m_{\pi} \cot \delta $",fontsize=16)
            # plt.title(f'{channel_1},{channel_2}  Scattering ',fontsize=16)     
            # plt.savefig(os.path.join(self.proj_handler.log_dir(), f'{channel}_Scattering.pdf') )
            #legend.set_title('Legend', prop={'size': 12})  # Set legend title and font size

            #plt.show()
        return print(f"Plotting  Complete")
        
        # make the plots, store in self.proj_dir_handler.plot_dir(), again, any log/error warnings go in self.proj_dir_handler.log_dir() as well (if any)
def get_particle_name(particle_str):
    return particle_str.split("(")[0]