## PLotting for Single-Channel scattering only 

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
