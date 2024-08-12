############################################################################################
##                 Date: 10/01/2023
##          
##          Purpose: Read HDF5 Data for LQCD
#           Reads Lattice Size from input, and file path to hdf5
#           Assume HDF5 has irreps with data in format
#           PSQ#,IRREP, LEVEL
#           NEEDS: to be generalized, ASAP
############################################################################################
##      Required packages
import h5py as h5               # library to open lattice data
import logging                  # incorporate logging for errors
#ends code when run logging.error(message) or logging.critical(message)
# logging.warning(message) and logging.debug(message) won't end code but will output at certain verbosity settings
import math                     # math library, mostly use it for math.pi, math.pow()          
import numpy as np              #basic functions, linear algebra, etc. 
import os


############################################################################################
# Set up logging
# logging.basicConfig(level=logging.CRITICAL)

# doc = '''
# H5Data - a task a read in h5py from correlator analysis
# Must have h5 from correlator analysis in project directory 
# required inputs
# -------------
# general:
#   ensemble_id: cls21_c103
#   project_dir
#   path to h5 file
# assuming hdf5 file is coming in like traditional sigmond data
# main directory is psq and single_hadrons

# '''

class LQCD_DATA_READER:
    
    def __init__(self, file_path,L, isospin=None, strangeness=None):
        self.file_path = file_path # file path that leads directo
        self.data = None  # You can initialize data here if needed
        self.L = L
        # check if file_path exists
        if not file_path:
            logging.critical("Ensure file path is in project directory")
        
        if isospin==None and strangeness==None:
            self.channel = None
        else:
            strange_str = f"S{strangeness}"
            if strangeness<0:
                strange_str = f"Sm{-strangeness}"

            self.channel = f"iso{isospin}_{strange_str}"


        # if not os.path.exists(self.file_path):
        #     logging.critical(f"The HDF5 file does not exist at the specified path: {self.file_path}")
        # else:
        #     pass
    # Continue with reading the file using h5py

    def load_data(self):
        # Open the H5py file and load the data
        self.data = h5.File(self.file_path,'r')

        #check the channel name, retrieve if one is available
        if self.channel==None:
            available_channels = list(self.data.keys())
            available_channels.remove('single_hadrons')
            if len(available_channels)==1:
                self.channel = available_channels[0]
            elif len(available_channels)>1:
                logging.critical(f"Too many channels in input file, please define a channel using 'isospin' and 'strangeness' parameters.")
            else:
                logging.critical(f"No channels available in input data file.")
        else:
            if self.channel not in self.data.keys():
                logging.critical(f"Channel {self.channel} not found in file {self.file_path}. Please correct name or data file.")

        return self.data

    def load_psq(self):
        if self.data is None:
            self.load_data()
        
        psq_keys = []
        for key in self.data[self.channel]:
            if key[0:3] == 'PSQ':
                psq_keys.append(key)
        return psq_keys

    def irrep_keys(self,psq):
        if self.data is None:
            self.load_data()

        return list(self.data[self.channel][psq].keys())

    def single_hadron_list(self):
        if self.data is None:
            self.load_data()
        
        single_hadron_list = [] #create a list with keys of hadrons in this analysis
        for key in self.data.get('single_hadrons'):
            single_hadron_list.append(key)

        return single_hadron_list
    
    def single_hadron_data(self,hadron):
        #returns full data for hadron in list
        #average data is 0, while BS data is 1:
        if self.data is None:
            self.load_data()  # Load data if not already loaded
        # term is a key for h5py files, can be accessed through .keys()  
        # single hadron has hadrons
        # L (Lambda), N (Nucleon), S (sigma), X (?), k (kaon), pi (pion)
        return self.data.get('single_hadrons')[hadron][:]
    
    # need to change, [0] data is the average, rest is bootstrap
    def ref_energy_data(self,PSQ,Irrep,level): 
        if self.data is None:
            self.load_data()
        # PSQ can be PSQ0,1,2,3
        # Irrep is the FV irrep 'G1u',...
        #level = f"ecm_{level}_ref"
        return self.data[self.channel].get(PSQ)[Irrep].get(level)[:]

    def free_levels(self,mom,irr,level):
        # particle_mom = self.data[self.channel][mom][irr].attrs['free_levels'][level] #pi(0)
        # particle, mom = particle_mom.split('(') #pi, 0)
        # mom = mom[:-1] #0
        return self.data[self.channel][mom][irr].attrs['free_levels'][level]
    
    # def energy_data_bootstrap_samples(self,PSQ,Irrep,level): 
    #     if self.data is None:
    #         self.load_data()
    #     # PSQ can be PSQ0,1,2,3
    #     # Irrep is the FV irrep 'G1u',...
    #     return self.data.get(PSQ)[Irrep].get(level)[1:]
    
    # def energy_data_load(self):
    #     if self.data is None:
    #         self.load_data()

    #     mom2 = ['PSQ0','PSQ1','PSQ2','PSQ3']

    #     # Define a mapping of mom values to irreps
    #     irreps = {
    #         'PSQ0': 'G1u',
    #         'PSQ1': 'G1',
    #         'PSQ2': 'G',
    #         'PSQ3': 'G',
    #     }

    #     lvls = [4,4,3,4]
    #     data_list = []
    #     for mom, lvl in zip(mom2, lvls):
    #         # Determine the associated irrep based on mom using the mapping
    #         irr = irreps.get(mom)

    #         # Determine the range of label numbers based on the level
    #         if mom == 'PSQ0':
    #             # For PSQ0, use the full range including 0
    #             label_numbers = range(lvl)
    #         else:
    #             # For other mom values, skip 0 and start from 1
    #             label_numbers = range(1, lvl + 1)

    #         # Iterate through the label_numbers range
    #         for label_num in label_numbers:
    #             # Construct the label based on mom, irr, and label_num
    #             label = f"ecm_{label_num}_ref"
    #             data_array = self.energy_data(mom,irr,label)
    #             data_list.append(data_array)

    #     data = np.array(data_list)
    #     return data

    # def energy_shift_data(self):
    #     if self.data is None:
    #         self.load_data()

    #     ecm_data = self.sigma_pi_boot_data() #15 

    #     mpi = np.array(self.pi_masses())[1:]
    #     mpi_ref , mS_ref = self.sigma_pi_ref_masses()
    #     mk_ref, mN_ref = self.nucleon_kaon_ref_masses()
        
    #      # mapping for masses
    #     mass_map = {
    #         'pi': np.array(mpi_ref)[1:],
    #         'S': np.array(mS_ref)[1:],
    #         'k': np.array(mk_ref)[1:],
    #         'N': np.array(mN_ref)[1:],
    #     }
        
    #     def extract_values(input_str):
    #         # Find the position of the opening and closing parentheses
    #         open_paren = input_str.find('(')
    #         close_paren = input_str.find(')')
            
    #         if open_paren != -1 and close_paren != -1:
    #             # Extract the part before the opening parenthesis
    #             part_before_paren = input_str[:open_paren] #Mass from mass_map
                
    #             # Extract the number inside the parentheses and convert it to an integer
    #             number_inside_paren = int(input_str[open_paren + 1:close_paren]) # d^2
                
    #             return part_before_paren, number_inside_paren
    #         else:
    #             return None, None  # Return None for both values if parentheses are not found

    #     def deltaE(ecm,ma,mb,n,m,psq):#function to shift e_cm data to shifted data to free energy levels
    #         if psq == 0:
    #             l = self.L*mpi
    #             dE = ecm - np.sqrt(ma**2 + n*(2*math.pi/l)**2 ) - np.sqrt((mb)**2 + m*(2*math.pi/l)**2 )
    #         else:
    #             l = self.L*mpi
    #             elab  = np.sqrt((ma)**2 + n*(2*math.pi/l)**2) + np.sqrt((mb)**2 + m*(2*math.pi/l)**2)
    #             ecmfree = np.sqrt(elab**2 - psq*(2*math.pi/(l))**2) 
    #             dE = ecm - ecmfree
    #         return dE
        
    #     mom2 = ['PSQ0','PSQ1','PSQ2','PSQ3']
    #     # Define a mapping of mom values to irreps
    #     irreps = {
    #         'PSQ0': 'G1u',
    #         'PSQ1': 'G1',
    #         'PSQ2': 'G',
    #         'PSQ3': 'G',
    #     }
    #     lvls = [1,1,1,1]
    #     data_list = []
    #     total_label_num = -1  # Initialize total label number
    #     # now need to get the energy shift from two-particle free energy
    #     # these are given in hdf5 file attributes
    #     #print(hf['PSQ3/G'].attrs['free_levels']) example of free energy shifts
    #     for mom, lvl in zip(mom2, lvls):
    #         # Determine the associated irrep based on mom using the mapping
    #         irr = irreps.get(mom)

    #         # Determine the range of label numbers based on the level
    #         if mom == 'PSQ0':
    #             # For PSQ0, use the full range including 0
    #             label_number = 0
    #         else:
    #             # For other mom values, skip 0 and start from 1
    #             label_number = 1

    #         total_label_num += 1
    #         # get energy from ecm_data
    #         ma, n = extract_values(self.data[mom][irr].attrs['free_levels'][label_number][0])
    #         mb, m = extract_values(self.data[mom][irr].attrs['free_levels'][label_number][1])
    #         ma = mass_map[ma]
    #         mb = mass_map[mb]
    #         #data_array = self.energy_data(mom,irr,label)
    #         data_list.append(deltaE(ecm_data[total_label_num],ma,mb,n,m,int(mom[3:])))
        
    #     data = np.array(data_list)
    #     return data

    
# first a class to access the h5py data
# class H5Data:
#     def __init__(self, file_path,L):
#         self.file_path = file_path
#         self.data = None  # You can initialize data here if needed
#         self.L = L

#     def load_data(self):
#         # Open the H5py file and load the data
#         self.data = h5.File(self.file_path)
    
#     def single_hadron_data(self,term):
#         if self.data is None:
#             self.load_data()  # Load data if not already loaded
#         # term is a key for h5py files, can be accessed through .keys()  
#         # single hadron has hadrons
#         # L (Lambda), N (Nucleon), S (sigma), X (?), k (kaon), pi (pion)
#         return self.data.get('single_hadrons')[term][:]
    
#     # need to change, [0] data is the average, rest is bootstrap
#     def energy_data(self,PSQ,Irrep,label): 
#         if self.data is None:
#             self.load_data()
#         # PSQ can be PSQ0,1,2,3
#         # Irrep is the FV irrep 'G1u',...
#         return self.data.get(PSQ)[Irrep].get(label)[:]
    
#     def pi_masses(self):
#         if self.data is None:
#             self.load_data()
#         # the masses required for calculation of reduced L
#         mpi = self.single_hadron_data('pi(0)').tolist() #ref is in units of mpi
#         return mpi
    
#     def sigma_pi_ref_masses(self):
#         if self.data is None:
#             self.load_data()
#         # the masses required for calculations including SigmaPi
#         mpi = self.single_hadron_data('pi(0)_ref').tolist() #ref is in units of mpi
#         mS = self.single_hadron_data('S(0)_ref').tolist() #ref is in units of mpi
#         return mpi, mS
    
#     def sigma_pi_data(self):
#         if self.data is None:
#             self.load_data()
#         # output the energy values needed, _ref means its E/mpi
#         E1 = self.energy_data('PSQ0','G1u','ecm_0_ref').tolist()
#         E2 = self.energy_data('PSQ1','G1','ecm_1_ref').tolist()
#         E3 = self.energy_data('PSQ2','G','ecm_1_ref').tolist()
#         E4 = self.energy_data('PSQ3','G','ecm_1_ref').tolist()
#         # combine them for energy_data
#         energy_data = [E1,E2,E3,E4]
#         return energy_data 
    
#     def sigma_pi_boot_data(self):
#         if self.data is None:
#             self.load_data()
#         # output the energy values needed, _ref means its E/mpi
#         E1 = self.energy_data('PSQ0','G1u','ecm_0_ref').tolist()[1:]
#         E2 = self.energy_data('PSQ1','G1','ecm_1_ref').tolist()[1:]
#         E3 = self.energy_data('PSQ2','G','ecm_1_ref').tolist()[1:]
#         E4 = self.energy_data('PSQ3','G','ecm_1_ref').tolist()[1:]
#         # combine them for energy_data
#         energy_data = [E1,E2,E3,E4]
#         return energy_data 
    
#     def covariance_data(self):
#         if self.data is None:
#             self.load_data()
#         # covariance matrix of bootstrap samples
#         sp_data = np.array(self.sigma_pi_data())
#         return np.cov(sp_data[:,1:])