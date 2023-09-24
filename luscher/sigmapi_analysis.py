############################################################################################
##                 Date: 08/01/2023
##              D250 LQCD Data Analysis: Sigma-Pi Scattering in m+pi = 200 MeV
##          Purpose: Resonance analysis for L(1405), this is single-scattering channel 
############################################################################################
##      Required packages
import cmath                    #math library (complex)
import csv                      # read csv files
import h5py as h5               # library to open lattice data
import math                     # math library, mostly use it for math.pi, math.pow()
import matplotlib.pyplot as plt #plotting library            
import numpy as np              #basic functions, linear algebra, etc. 
import random
import scipy.special as sp      #special functions
import scipy.integrate as integrate #needed to do numerical integration, alt?
##
from matplotlib.lines import Line2D  # Import Line2D for custom legend handles
from tqdm import tqdm  # Import the tqdm library
from scipy import optimize,integrate
from scipy.optimize import fsolve,minimize
from scipy.integrate import quad
from zeta import Z
## Z function that is fast
#%run create_momentum_array #only need to run this once, might be broken in newest python
#%run zeta.py
############################################################################################
# first a class to access the h5py data

class H5Data:
    def __init__(self, file_path,L=64):
        self.file_path = file_path
        self.data = None  # You can initialize data here if needed
        self.L = 64

    def load_data(self):
        # Open the H5py file and load the data
        self.data = h5.File(self.file_path)
    
    def single_hadron_data(self,term):
        if self.data is None:
            self.load_data()  # Load data if not already loaded
        # term is a key for h5py files, can be accessed through .keys()  
        # single hadron has hadrons
        # L (Lambda), N (Nucleon), S (sigma), X (?), k (kaon), pi (pion)
        return self.data.get('single_hadrons')[term][:]
    
    # need to change, [0] data is the average, rest is bootstrap
    def energy_data(self,PSQ,Irrep,label): 
        if self.data is None:
            self.load_data()
        # PSQ can be PSQ0,1,2,3
        # Irrep is the FV irrep 'G1u',...
        return self.data.get(PSQ)[Irrep].get(label)[:]
    
    def pi_masses(self):
        if self.data is None:
            self.load_data()
        # the masses required for calculation of reduced L
        mpi = self.single_hadron_data('pi(0)').tolist() #ref is in units of mpi
        return mpi
    
    def sigma_pi_ref_masses(self):
        if self.data is None:
            self.load_data()
        # the masses required for calculations including SigmaPi
        mpi = self.single_hadron_data('pi(0)_ref').tolist() #ref is in units of mpi
        mS = self.single_hadron_data('S(0)_ref').tolist() #ref is in units of mpi
        return mpi, mS
    
    def sigma_pi_data(self):
        if self.data is None:
            self.load_data()
        # output the energy values needed, _ref means its E/mpi
        E1 = self.energy_data('PSQ0','G1u','ecm_0_ref').tolist()
        E2 = self.energy_data('PSQ1','G1','ecm_1_ref').tolist()
        E3 = self.energy_data('PSQ2','G','ecm_1_ref').tolist()
        E4 = self.energy_data('PSQ3','G','ecm_1_ref').tolist()
        # combine them for energy_data
        energy_data = [E1,E2,E3,E4]
        return energy_data 
    
    def sigma_pi_boot_data(self):
        if self.data is None:
            self.load_data()
        # output the energy values needed, _ref means its E/mpi
        E1 = self.energy_data('PSQ0','G1u','ecm_0_ref').tolist()[1:]
        E2 = self.energy_data('PSQ1','G1','ecm_1_ref').tolist()[1:]
        E3 = self.energy_data('PSQ2','G','ecm_1_ref').tolist()[1:]
        E4 = self.energy_data('PSQ3','G','ecm_1_ref').tolist()[1:]
        # combine them for energy_data
        energy_data = [E1,E2,E3,E4]
        return energy_data 
    
    def covariance_data(self):
        if self.data is None:
            self.load_data()
        # covariance matrix of bootstrap samples
        sp_data = np.array(self.sigma_pi_data())
        return np.cov(sp_data[:,1:])



class QC:
    # these are for the single-channel case
    def __init__(self,data):
        self.data = data
        self.mpi_ref, self.mS_ref = self.data.sigma_pi_ref_masses()
        self.mpi = self.data.pi_masses()
    # function to define momentum state for Irrep. d^2 = 1 is [0,0,1], ex
    # by symmetry, [1,0,0] should give the same result
    def momentum_state(self,i):
        if i == 0:
            return np.array([i,i,i])
        elif i == 1:
            return np.array([0,0,i])
        elif i == 2:
            return np.array([i-1,i-1,0])
        elif i == 3:
            return np.array([i/3,i/3,i/3])
        else: 
            # Raise an exception for invalid input
            raise ValueError("Invalid value for 'i'. 'i' must be 0, 1, 2, or 3.")
    # each entry gives E_cm, energy in the center of mass, / mpi
    # need a function to convert to q^2 
    # q^2 is in COM frame, and depends on masses
    def q2(self,ecm,ind):
        #ecm is the energy data
        #mpi, mS = self.data.sigma_pi_masses()
        q2 = ecm**2 / 4 - (self.mpi_ref[ind]**2 + self.mS_ref[ind]**2) / 2 + ((self.mpi_ref[ind]**2 - self.mS_ref[ind]**2)**2) / (4*ecm**2)
        return q2

    # delta is distance to threshold in energy, for fitting
    def delta_Sp(self,ecm,ind):
        #mpi, mS = self.data.sigma_pi_masses()
        return (ecm**2 - (self.mpi_ref[ind]+self.mS_ref[ind])**2 )/ (self.mpi_ref[ind]+self.mS_ref[ind])**2
    
    #gamma, lorentz factor
    def gamma(self,ecm,d,ind):
        d_vec = self.momentum_state(d) 
        L_ref = self.data.L*self.mpi[ind]
        E = math.sqrt(ecm**2 + (((2*math.pi)/L_ref)**2)*np.dot(d_vec,d_vec))
        return E/ecm
    # mass different for non-identical particles
    def msplit(self,ecm,ind):
        return (1 + (self.mS_ref[ind]**2 - self.mpi_ref[ind]**2)/ecm**2 )
    
    def qcotd(self,ecm,d,ind):
        L_ref = self.data.L*self.mpi[ind]
        d_vec = self.momentum_state(d) 
        c = 2 / (self.gamma(ecm,d,ind)*L_ref*math.sqrt(math.pi))
        return c*Z(self.q2(ecm,ind)*((L_ref/(2*math.pi))**2),gamma=self.gamma(ecm,d,ind),l=0,m=0,d=d_vec,m_split=self.msplit(ecm,ind),precision=1e-11).real

    

class Chi2Fit:
    def __init__(self,data):
        # Create an instance of the H5Data class
        self.data = data
        #self.data = H5Data(data_path)
        #self.energy_data_full = self.data.sigma_pi_data() #bootstrap samples, 800
        self.average_data = np.array(self.data.sigma_pi_data())[:,:1]
        self.energy_data = self.data.sigma_pi_boot_data()
        self.cov = self.data.covariance_data()
        self.qc = QC(self.data)
        self.mpi = self.qc.mpi
        self.mpi_ref = self.qc.mpi_ref
        self.mS_ref = self.qc.mS_ref
        self.sample_length = 100
        self.x0 = [20,5]

    def ere(self,ecm,ind,a,b):
        return ((1/a)+0.5*b*self.qc.q2(ecm,ind))
    
    def ere_delta(self,ecm,ind,a,b):
        return (ecm)*(a+b*self.qc.delta_Sp(ecm,ind))

    def deter(self,ecm,d,ind,a,b):
        return self.qc.qcotd(ecm,d,ind) - self.ere(ecm,ind,a,b)
    
    
    def QC1(self, d, ind, a, b):
        func = lambda ecm: self.deter(ecm, d, ind, a, b)
        return fsolve(func,self.energy_vec[d])[0]
    
    # def QC1(self, d, ind, a, b): #gives worse results, for some reason
    #     root = []
    #     tol = 1e-2
    #     ecm_min = 6.5
    #     ecm_max = 7.1
    #     num_points = 50
    #     ecm_guesses = np.linspace(ecm_min, ecm_max, num_points)

    #     previous_sign = np.sign(self.deter(ecm_guesses[0], d, ind, a, b))

    #     for i in range(1, len(ecm_guesses)):
    #         current_sign = np.sign(self.deter(ecm_guesses[i], d, ind, a, b))

    #         if current_sign != previous_sign:
    #             # The function crosses zero between ecm_guesses[i-1] and ecm_guesses[i]
    #             # Use fsolve at the midpoint of the interval
    #             ecm_midpoint = (ecm_guesses[i-1] + ecm_guesses[i]) / 2.0
    #             def det_func(ecm):
    #                 return self.deter(ecm, d, ind, a, b)
    #             res = fsolve(det_func, ecm_midpoint)
    #             root.append(res[0])  # Append the root to the list
    #             break  # Exit the loop after finding one root

    #         previous_sign = current_sign

    #     return np.array(root)
    
    # def energies(self):
    #     for j in range(len(self.energy_data)):
    #         self.energy_samples[j] = self.energy_data[j][0:self.sample_length]
    def energy_truncate(self, num_samples):
        self.energy_samples = [[] for _ in range(len(self.energy_data))]

        # Calculate the indices for the first data point (j=0)
        first_data = self.energy_data[0]
        self.selected_indices = random.sample(range(len(first_data)), min(num_samples, len(first_data)))

        # Create temporary lists to store shortened vectors
        mpi_ref_shortened = []
        mS_ref_shortened = []
        mpi_shortened = []

        for j in range(len(self.energy_data)):
            data = self.energy_data[j]

            # Extract the corresponding data points using the indices from the first data point
            self.energy_samples[j] = [data[i] for i in self.selected_indices]

            # Append the shortened vectors to the temporary lists
            if j == 0:
                mpi_ref_shortened.extend([self.mpi_ref[i] for i in self.selected_indices])
                mS_ref_shortened.extend([self.mS_ref[i] for i in self.selected_indices])
                mpi_shortened.extend([self.mpi[i] for i in self.selected_indices])

        # Assign the shortened vectors to their respective attributes outside the loop
        self.mpi_ref = mpi_ref_shortened
        self.mS_ref = mS_ref_shortened
        self.mpi = mpi_shortened

    # def energy_truncate(self, num_samples):
    #     self.energy_samples = [[] for _ in range(len(self.energy_data))]

    #     # Create temporary lists to store shortened vectors
    #     mpi_ref_shortened = []
    #     mS_ref_shortened = []
    #     mpi_shortened = []

    #     for j in range(len(self.energy_data)):
    #         data = self.energy_data[j]
    #         mean_value = self.average_data[j]

    #         # Calculate the absolute differences between each data point and the mean value
    #         abs_diff = np.abs(data - mean_value)

    #         # Sort the indices based on the absolute differences
    #         sorted_indices = np.argsort(abs_diff)

    #         # Select the 200 closest indices (if available)
    #         selected_indices = sorted_indices[:num_samples]

    #         # Extract the corresponding data points and assign them to self.energy_samples[j]
    #         self.energy_samples[j] = [data[i] for i in selected_indices]

    #         # Append the shortened vectors to the temporary lists
    #         if j == 0:
    #             mpi_ref_shortened.extend([self.mpi_ref[i] for i in selected_indices])
    #             mS_ref_shortened.extend([self.mS_ref[i] for i in selected_indices])
    #             mpi_shortened.extend([self.mpi[i] for i in selected_indices])

    #     # Assign the shortened vectors to their respective attributes outside the loop
    #     self.mpi_ref = mpi_ref_shortened
    #     self.mS_ref = mS_ref_shortened
    #     self.mpi = mpi_shortened

    # minimizing chi2_ave took 1min25s, giving a=22.371,b=.2250
    def chi2_ave(self,x):
        a,b = x
        res = np.zeros(len(self.energy_data))
        self.energy_vec = self.average_data
        for irrep in range(len(self.energy_data)):
            res[irrep] = self.average_data[irrep] - self.QC1(irrep,0,a,b)
        value = res@np.linalg.inv(self.cov)@res
        return value

    def chi2(self,x):
        a,b = x
        res = np.zeros(len(self.energy_data))
        for irrep in range(len(self.energy_data)):
            res[irrep] = self.energy_vec[irrep] - self.QC1(irrep,self.ind,a,b)
        value = res@np.linalg.inv(self.cov)@res
        return value
    # 
    def bootstrap_fit(self, output_file='fitresults_100_imp.txt'):
        self.energy_truncate(self.sample_length)

        a_results = []
        b_results = []
        chi2_results = []

        try:
            with open(output_file, 'w') as file:
                # Create a tqdm progress bar
                progress_bar = tqdm(range(self.sample_length))

                for i in progress_bar:
                    self.ind = i
                    self.energy_vec = [item[i] for item in self.energy_samples]

                    # Flag to indicate whether we need to retry with adjusted parameters
                    retry = True

                    while retry:
                        try:
                            # Create a Gaussian distribution around the parameters you choose
                            #std_deviation = [5, 2.5]  # Calculate the width from data
                            #a_guess = random.uniform(14, 23)
                            #b_guess = random.uniform(0.5, 8)
                            # Use the calculated standard deviations as the errors for Gaussian distributions
                            a_guess = np.random.normal(22.5, np.std(self.energy_samples))
                            b_guess = np.random.normal(2.5, np.std(self.energy_samples))
                            # Generate random values from the Gaussian distribution
                            #self.x0g = np.random.normal(self.x0, std_deviation)

                            result = minimize(self.chi2, x0=[a_guess,b_guess], method='nelder-mead')

                            if result.success:
                                a_results.append(result.x[0])
                                b_results.append(result.x[1])
                                chi2_results.append(result.fun)

                                # Write the results to the file
                                file.write(f"{result.x[0]} {result.x[1]} {self.selected_indices[self.ind]}\n")
                                file.flush()

                                # If we reached here without exceptions and success is True, no need to retry
                                retry = False
                            else:
                                # If minimization was not successful, retry with new parameters
                                print(f"Minimization was not successful. Retrying with different parameters.")
                        except Exception as e:
                            print(f"An error occurred: {e}")

        except Exception as e:
            # Handle any exceptions (e.g., file write errors)
            print(f"An error occurred: {e}")
        

        return [a_results, b_results, chi2_results]

    def deriv(self,n,i,x): # 2 parameter difference derivative
        self.energy_vec = self.average_data
        a, b = x 
        eps = .001
        if n == 0:
            return (self.QC1(i,0,a-eps,b)-self.QC1(i,0,a+eps,b))/(2*eps)
        elif n ==1:
            return (self.QC1(i,0,a,b-eps)-self.QC1(i,0,a,b+eps))/(2*eps)
        
    def vij(self,x): #V_ij is error matrix for parameters, take _ii to get each error
        # v_nm = dp_i / dp_n
        # error estimation function
        # cov = self.cov = self.data.covariance_data()
        nint = [0,1]
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
            
        
        return lmat
            



  

class Plotter:
    def __init__(self, data):
        self.data = data


    def phase_shift_plot(self):
        self.fit = Chi2Fit(self.data)
        average_energies = self.fit.average_data
        x = []
        x_range = []
        y = []
        y_range = []
        for i in range(len(average_energies)):
            x.append(self.fit.qc.q2(average_energies[i], 0))
            y.append(self.fit.qc.qcotd(average_energies[i], i, 0))
            if i == 0:
                x_range_i = []
                y_range_i = []
                for en in np.linspace(average_energies[i] - 0.01, average_energies[i] + 0.01, 100):
                    x_range_i.append(self.fit.qc.q2(en, 0))
                    y_range_i.append(self.fit.qc.qcotd(en, i, 0))
                x_range.append(x_range_i)
                y_range.append(y_range_i)
            else:
                xp = []
                yp = []
                for en in np.linspace(average_energies[i] - 0.01, average_energies[i] + 0.01, 100):
                    xp.append(self.fit.qc.q2(en, 0))
                    yp.append(self.fit.qc.qcotd(en, i, 0))
                x_range.append(xp)
                y_range.append(yp)

        #ecm_values = np.linspace(min(average_energies),max(average_energies),100)
        q2_values = np.linspace(-0.25, 0.1, 150)
        virtual_state = []
        for q2 in q2_values:
            virtual_state.append(cmath.sqrt(-q2))
        
        x_list = [array[0] for array in x]
        y_list = [array[0] for array in y]
        fit= np.polyfit(np.array(x_list),np.array(y_list),1)

        line = []
        #fit = [-(1/3.30),2*1.582] #with ecm factor
        for q2 in q2_values:
            line.append(np.polyval(fit,q2))
        results_file = './fitresults_100.txt'
        # load data from the bootstrap file, currently only using 100 
        data_bs = np.loadtxt(results_file) 

        # Extract the parameters (a_results and b_results) from the loaded data
        a_results = data_bs[:, 0]
        b_results = data_bs[:, 1]

        # Calculate the covariance matrix from the parameters
        parameter_matrix = np.array([a_results, b_results])
        covariance_matrix = np.cov(parameter_matrix)
        vij = self.fit.vij([3.30,1.582])
        errs = vij@np.linalg.inv(covariance_matrix)@np.transpose(vij)
        
        fit_min = np.polyval(fit - 0.5*errs,q2_values)
        fit_max = np.polyval(fit + 0.5*errs,q2_values)

        plt.figure(figsize=(8,6))
        shapes = ['o','s','D','v']
        labels = ['$G_{1u}(0)$','$G_1 (1)$', '$G (2)$' , '$G(3)$']
        legend_handles = []  # Create an empty list to store custom legend handles
        for i in range(len(average_energies)):
            plt.plot(x[i], y[i], marker=shapes[i], color='blue', label=labels[i])
            plt.plot(x_range[i], y_range[i], color="blue", alpha=0.8)  # Plot the ranges with transparency
            # Add a custom legend handle (marker with no line)
            legend_handles.append(Line2D([0], [0], marker=shapes[i], color='w', markerfacecolor='blue', markersize=10, label=labels[i]))


        plt.plot(q2_values,virtual_state,color='black',linestyle='--')
        plt.plot(q2_values,line, color='blue',linestyle='-.')
        plt.fill_between(q2_values,fit_min,fit_max,alpha = 0.5, color = 'lightblue')
        plt.axhline(y=0,color='black')
        plt.axvline(x=0,color='black')
        plt.ylim(0,0.6)
        plt.xlim(-0.22,0.05)
        # Customize the legend with custom handles (markers only)
        legend = plt.legend(handles=legend_handles, loc='lower right', title='Legend', prop={'size': 12})
        #legend.set_title('Legend', prop={'size': 12})  # Set legend title and font size

        plt.show()
        
        


    def phase_shift_plot_error(self,results_file, best_fit_a, best_fit_b, num_points=100, confidence_level=95):
        # load data from the bootstrap file, currently only using 100 
        data_bs = np.loadtxt(results_file) 

        # Extract the parameters (a_results and b_results) from the loaded data
        a_results = data_bs[:, 0]
        b_results = data_bs[:, 1]

        # Calculate the covariance matrix from the parameters
        parameter_matrix = np.array([a_results, b_results])
        covariance_matrix = np.cov(parameter_matrix)

        # Number of points for the error bands
        q2_values = np.linspace(-0.25, 0.1, num_points)

        # Generate deviations from the best-fit parameters using the covariance matrix
        np.random.seed(42)  # for reproducibility
        parameter_deviations = np.random.multivariate_normal([0, 0], covariance_matrix, num_points)

        # Calculate the best-fit line without parameter variations
        best_fit_line = (-1 / best_fit_a) + 0.5 * best_fit_b * q2_values
        # Calculate y_values using the equation with perturbed 'a' values
        y_values = (-1 / (best_fit_a + parameter_deviations[:, 0])) + 0.5 * (best_fit_b + parameter_deviations[:, 1]) * q2_values
        # Plot the best-fit line and error bands
        plt.figure(figsize=(8, 6))
        plt.plot(q2_values, best_fit_line, color='blue', linestyle="-." ,linewidth=2)
        plt.fill_between(q2_values, np.percentile(y_values, 16, axis=0), np.percentile(y_values, 84, axis=0),
                        color='lightblue', alpha=0.5)
        plt.xlabel('q^2')
        plt.ylabel('y')
        plt.legend()
        plt.title('Best Fit Line with 1 Sigma (68%) Error Bands')
        plt.grid(True)
        plt.show()


    def plot_params_bootstrap(self):
        results_file = './fitresults_100.txt'
        # load data from the bootstrap file, currently only using 100 
        data_bs = np.loadtxt(results_file) 

        # Extract the parameters (a_results and b_results) from the loaded data
        a_results = data_bs[:, 0]
        b_results = data_bs[:, 1]

        plt.figure(figsize=(8, 6))
        plt.hist(a_results, bins = 20, color = 'red', alpha = 0.5, label = "a")
        plt.hist(b_results, bins = 20, color = 'blue', alpha = 0.5, label = "b")
        plt.title('Bootstrap Parameter Samples: 100')
        plt.legend(fontsize=14)
        plt.show()
    
    def plot_params_bootstrap_points(self):
        results_file = './fitresults_100.txt'
        # load data from the bootstrap file, currently only using 100 
        data_bs = np.loadtxt(results_file) 

        # Extract the parameters (a_results and b_results) from the loaded data
        a_results = data_bs[:, 0]
        b_results = data_bs[:, 1]

        plt.figure(figsize=(8, 6))
        plt.plot( b_results,'o')
        #plt.hist(b_results, bins = 20, color = 'blue', alpha = 0.5, label = "b")
        plt.title('Bootstrap Parameter Samples: 100')
        plt.legend(fontsize=14)
        plt.show()

    def plot_energy_histograms(self,separate_plots=False):
        E1, E2, E3, E4 = self.data.sigma_pi_data()
        if separate_plots:
            plt.figure(figsize=(12, 6))  # Set the figure size

            plt.subplot(2, 2, 1)
            plt.hist(E1, bins=50, color='blue', alpha=0.7)
            plt.title('$G_{1u} (0)$')

            plt.subplot(2, 2, 2)
            plt.hist(E2, bins=50, color='green', alpha=0.7)
            plt.title('$G_1 (1)$')

            plt.subplot(2, 2, 3)
            plt.hist(E3, bins=50, color='red', alpha=0.7)
            plt.title('$G(2)$')

            plt.subplot(2, 2, 4)
            plt.hist(E4, bins=50, color='purple', alpha=0.7)
            plt.title('$G(3)$')
            plt.tight_layout()  # Adjust subplot layout for better spacing
        else:
            # Create a single plot with all histograms
            plt.figure(figsize=(8, 6))  # Set the figure size

            plt.hist(E1, bins=50, color='blue', alpha=0.5, label='$G_{1u} (0)$')
            plt.hist(E2, bins=50, color='green', alpha=0.5, label='$G_1 (1)$')
            plt.hist(E3, bins=50, color='red', alpha=0.5, label='$G(2)$')
            plt.hist(E4, bins=50, color='purple', alpha=0.5, label='$G(3)$')

            plt.title('Bootstrap Energy Samples')
            #plt.xlabel('Value')
            #plt.ylabel('Frequency')
            plt.legend(fontsize=14)
        # Show the plots
        plt.show()

    def plot_energy_histograms_1sigma(self, separate_plots=True):
        E1, E2, E3, E4 = self.data.sigma_pi_data()

        if separate_plots:
            # Separate plots for each dataset
            datasets = [(E1, '$G_{1u} (0)$', 'blue'),
                    (E2, '$G_1 (1)$', 'green'),
                    (E3, '$G(2)$', 'red'),
                    (E4, '$G(3)$', 'purple')]

            plt.figure(figsize=(12, 6))

            for i, (data, label, color) in enumerate(datasets, start=1):
                plt.subplot(2, 2, i)
                plt.hist(data[1:], bins=50, color=color, alpha=0.7, label='Full Distribution')
                plt.title(label)
            
                mean_value = data[0]
                std_deviation = np.std(data[1:])
                lower_bound = mean_value - std_deviation
                upper_bound = mean_value + std_deviation
            
                sigma_data = [x for x in data if lower_bound <= x <= upper_bound]
            
                plt.hist(sigma_data, bins=50, color='black', alpha=0.5, label='1 Sigma Distribution')
                plt.axvline(mean_value, color='black', linestyle='--')  # Add vertical dashed line at mean
                plt.legend()

            plt.tight_layout()

        else:
            # Create a single plot with all histograms
            plt.figure(figsize=(8, 6))

            datasets = [(E1, '$G_{1u} (0)$', 'blue'),
                    (E2, '$G_1 (1)$', 'green'),
                    (E3, '$G(2)$', 'red'),
                    (E4, '$G(3)$', 'purple')]

            for data, label, color in datasets:
                mean_value = np.mean(data)
                std_deviation = np.std(data)
                lower_bound = mean_value - std_deviation
                upper_bound = mean_value + std_deviation

                plt.hist(data, bins=50, color=color, alpha=0.5, label=label)
                plt.hist([x for x in data if lower_bound <= x <= upper_bound], bins=50, color=color, alpha=0.7)
                plt.axvline(mean_value, color=color, linestyle='--')  # Add vertical dashed line at mean
                plt.legend(fontsize=14)

            plt.title('Bootstrap Energy Samples')
    
    # Show the plots
        plt.show()






            
if __name__ == "__main__" :     
    data_path = './Data/rebin10_fer.hdf5'
    data = H5Data(data_path)
    pl = Plotter(data)
    result_file = './fitresults_100.txt'
    #a, b = [22.371,.2250]
    #pl.plot_params_bootstrap_points()
    pl.phase_shift_plot()
    #print(chi.chi2_ave([22.371,.2250])) # these are the mean values best fit
    #result = chi.bootstrap_fit() # run this one to generate text file with bootstrap results




#Plotter(Data)
# if __name__ == "__main__":
#     data_path = './Data/rebin10_fer.hdf5'  # Replace with the path to your data file
#     chi2_fit = Chi2Fit(data_path)
#     chi2_fit.QC1(0,0,20,0.7)
    #chi2_fit.energies()
    #chi2_fit.cov = np.cov(chi2_fit.energy_data)
    #chi2_fit.ind = 0
    #chi2_fit.energy_vec = [item[chi2_fit.ind ] for item in chi2_fit.energy_data]
    #print(chi2_fit.chi2([10,-0.5]))
    # if x0 is too small, gives overflow error! Need to find way to fix.
    #result = minimize(chi2_fit.chi2,x0=chi2_fit.x0,method = 'nelder-mead')
    #print(result)
    #results = chi2_fit.bootstrap_fit()

# #      Save the results to a file without any separation
#     with open("results.txt", "w") as file:
#         for a, b, chi2 in zip(results[0], results[1], results[2]):
#             file.write(f"{a}\t{b}\t{chi2}\n")