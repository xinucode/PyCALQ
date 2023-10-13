single_had_key = 'single_hadrons'

#designates colors and markers for plotting
colors = ['green','blue','orange','purple','brown','gray','darkblue','red','pink',
          'yellow','lightblue','lightgreen','#1f77b4', '#ff7f0e', 
          '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
markers = ['o','s','D','v','^','*','x','+']
zigzag_shifts = [-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1]
stripe_shifts = [0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1]
stripe4_shifts = [0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5]
stripe5_shifts = [0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1]

#list of files that are jackknife samples rather than bootstrap
jackknife_sampling_methods = ['isoquartet_nonstrange_fermionic\\qsqr_samplings_isoquartet_nonstrange_fermionic_colin_rebin20_jackknife.hdf5','isoquartet_nonstrange_fermionic\\qsqr_samplings_isoquartet_nonstrange_fermionic_colin_rebin20_jackknife_6-12.hdf5']

#Gives a numerical order for the irreps so that I can control how they show up in the plots
alphabetical = {'A1g':0.01107,'F1':0.0601,'F2':0.0602,'G': 0.07,'G1': 0.0701,'G1g': 0.07107,'G1u': 0.07121,'G2': 0.0702, 'Hg': 0.0807, 'Hu': 0.0821, 'T1g':0.20107, 'A1g':0.01107, 'A1gm':0.0110713, 'A2g':0.01207, 'A1u':0.01121, 'A2u':0.01221, 'Eg':0.0507, 'Eu':0.0521, 'T1u':0.20121, 'T2g':0.20207, 'T2u':0.20221, 'A1':0.011, 'A2':0.012, 'B1':0.021, 'B2':0.022, 'E':0.05}

#latex format of all irreps and particles
latex_format = {'A1g':r"$A_{1g}$", 'A1gm':r"$A_{1g}^{-}$", 'F1':r"$F_1$", 'F2':r"$F_2$", 'G': r'$G$', 'G1':r"$G_1$", 'G1g':r"$G_{1g}$", 'G1u':r"$G_{1u}$", 'G2':r"$G_2$", 'Hg':r"$H_{g}$", 'Hu':r"$H_{u}$", 'N': r'$N$', 'pi':r"$\pi$", 'T1g':r"$T_{1g}$", 'A1u':r"$A_{1u}$", 'A2u':r"$A_{2u}$", 'Eg':r"$E_{g}$", 'Eu':r"$E_{u}$", 'T1u':r"$T_{1u}$", 'T2g':r"$T_{2g}$", 'T2u':r"$T_{2u}$", 'A1':r"$A_1$", 'A2':r"$A_2$", 'B1':r"$B_1$", 'B2':r"$B_2$", 'E':r"$E$", 'S':r"$\Sigma$", 'k':r"$\bar{K}$", 'L': r"$\Lambda$", 'X': r"$\Xi$"}

#sigmond scripts tmin plot end tags
tmin_file_tags = {'0.agr':'single-exponential fit','4.agr':'double-exponential fit','8.agr':'geometric fit'}
dtmin_file_tags = {'0_D.agr':'single-exponential fit','4_D.agr':'double-exponential fit','8_D.agr':'geometric fit'}

fit_nicknames = {'single-exponential fit':'single','double-exponential fit':'double','geometric fit':'geometric','single-exponential ratio fit':'singleR','double-exponential ratio fit':'doubleR','geometric ratio fit':'geometricR'}
fit_longnames = {'single':'single-exponential fit','double':'double-exponential fit','geometric':'geometric fit','singleR':'single-exponential ratio fit','doubleR':'double-exponential ratio fit','geometricR':'geometric ratio fit'}

########## xmgrace ################

#xmgrace label format for particles #insert link here to help
xmgrace_format = {'N': 'N', 'pi': '\\xp'}

#settings for xmgrace plots
#legend_setting = "0.25, 0.92" #string of x,y location where legend should be located on graph
legend_setting = "0.7, 0.92" #string of x,y location where legend should be located on graph
standard_graph_width = 11.64 #graphs produced by sigmond scripts have these widths
standard_y_min = 2.33 #graphs produced by sigmond typically have this size of gap below them for axis labels and such
standard_graph_height = 12.42 #graphs produced by sigmond scripts have these heights
standard_xmax_gap = 15.52-standard_graph_width-3.1 #graphs produced by sigmond scripts have this right hand vertical gap
max_number_of_ticks = 6 #max number of labels on the y axis

spectrum_color_indices = [4,6,9,10,11,12,15,5,3,7]
y_label = '"aE\\s\Qfit"'
