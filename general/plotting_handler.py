
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pickle
import pylatex
import os
import logging
import numpy as np

import sigmond
import fvspectrum.sigmond_util as sigmond_util
import fvspectrum.spectrum_plotting_settings.settings as psettings
from sigmond_scripts.analysis.utils import util as utils
from sigmond_scripts.analysis.sigmond_info import fit_info


ctext_x = 0.3
ctext_y = 0.7
blank = patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                 lw=0, alpha=0)
def shift_levels( indexes, vals, errors, shifted_array=np.array([]), index=0 ):
    if len(indexes):
        if not shifted_array.any():
            shifted_array = np.array([0.0]*len(indexes))
            
        if index==len(indexes):
            return shifted_array
            
        these_indexes = indexes[index:]
        this_remaining_irrep_index = these_indexes[0]
        this_remaining_irrep_indexes = np.where( these_indexes==this_remaining_irrep_index )[0]+index
        
        if len(this_remaining_irrep_indexes)==1:
            return shift_levels( indexes, vals, errors, shifted_array, index=index+1 )
        else:
            shift = 1
            this_val_upper = vals[index]+errors[index]
            this_val_lower = vals[index]-errors[index]
            this_cluster = [index]
            for i in this_remaining_irrep_indexes[1:]:
                overlap = False  
                compare_upper = vals[i]+errors[i]
                compare_lower = vals[i]-errors[i]
                if shifted_array[i]!=0.0:
                    continue
                if compare_lower<=this_val_lower and compare_upper>=this_val_upper: #set new bounds on this_value
                    overlap=True
                elif compare_lower>=this_val_lower and compare_upper<=this_val_upper:
                    overlap=True
                elif compare_lower<=this_val_upper and compare_lower>=this_val_lower and compare_upper>=this_val_upper:
                    overlap=True
                elif compare_lower<=this_val_lower and compare_upper<=this_val_upper and compare_upper>=this_val_lower:
                    overlap=True
                if overlap:
                    shifted_array[i] = shift
                    this_cluster.append(i)
                    if shift>0:
                        shift = -shift
                    else:
                        shift = abs(shift)+1
            
            if len(this_cluster)%2==0 and len(this_cluster):
                for i in this_cluster:
                    shifted_array[i]-=0.5
            
            new_index = np.where(shifted_array[index+1:]==0.0)[0]
            if new_index.any():
                return shift_levels( indexes, vals, errors, shifted_array, index=new_index[0]+index+1 )
            else:
                return shifted_array
                        
    return 0.0

class PlottingHandler:

    def __init__(self):

        # Flag to check if LaTeX is available for plot labels
        self.latex = True
        #checks if latex compler is present, if not sets this variable to false
            #use this for latex labels, and make an option for non-latex labels.
        self.latex = sigmond_util.set_latex_in_plots(plt.style)
        self.doc = []

    #############################
    ##### matplotlib actions ####
    #############################

    def create_fig(self, figwidth, figheight):
        """Create a new Matplotlib figure with specified width and height."""
        self.fig = plt.figure(figsize=(figwidth, figheight)) 
        self.figheight = figheight
        self.figwidth = figwidth
        return self.fig
    
    def clf(self):
        """Clear the current Matplotlib figure."""
        plt.clf()
        plt.gcf().set_size_inches(self.figwidth, self.figheight)

    def set_figsize(self, figwidth, figheight):
        self.figheight = figheight
        self.figwidth = figwidth
        plt.gcf().set_size_inches(figwidth, figheight)

    def moving_textbox(self, labels):
        plt.legend([blank]*len(labels), labels, frameon=False, fancybox=False, framealpha=0.0, shadow=None)

    def correlator_plot(self,df, ptype=0, op1=None, op2=None, color_index = 0):
        """Generate a correlator plot using Matplotlib."""
        plt.errorbar(x=df["aTime"],y=df["FullEstimate"],yerr=df["SymmetricError"], linewidth=0.0, elinewidth=1.5, capsize=5, color=psettings.colors[color_index], marker=psettings.markers[color_index] )
        if ptype==0:
            plt.ylabel(r"$C(t)$") #not all dollar signs require a latex compiler, for example this is okay.
        else:
            if self.latex:
                plt.ylabel(r"$a_tE_{\textup{lab}}$") #but the use of "\textup{}" command requires a latex compiler
            else: 
                plt.ylabel(r"$a_tE_{lab}$")

        plt.xlabel(r"$t/a_t$")

        # Annotate the plot with additional information if provided
        labels = []
        if op1:
            labels.append(f"snk: {str(op1)}") #double check that I'm not messing up sink and source
        if op2:
            labels.append(f"src: {str(op2)}")
        if labels:
            self.moving_textbox(labels)

        plt.tight_layout()

    def sigmond_corrfit_plot(self,df, fit_result_info, Nt, ptype=0, op1=None, sh_index = 0, color_index = 0):
        """Generate a correlator plot with fit using Matplotlib."""
        labels = []
        plt.errorbar(x=df["aTime"],y=df["FullEstimate"],yerr=df["SymmetricError"], linewidth=0.0, elinewidth=1.5, capsize=5, color=psettings.colors[color_index], marker=psettings.markers[color_index], zorder=1)

        if ptype==0:
            if fit_result_info["success"]:
                x=np.linspace(fit_result_info["info"].tmin, fit_result_info["info"].tmax, 100)
                model = fit_result_info["info"].model.sigmond_object(Nt)
                params = []
                for estimate in fit_result_info["estimates"]:
                    params.append(estimate.getFullEstimate())
                y = [model.eval(params, xi) for xi in x]
                plt.plot(x,y, color="black", ls="--")
            if fit_result_info["info"].ratio:
                plt.ylabel(r"$R(t)$") #not all dollar signs require a latex compiler, for example this is okay.
            else:
                plt.ylabel(r"$C(t)$") #not all dollar signs require a latex compiler, for example this is okay.
        else:
            if fit_result_info["success"]:
                #for single -> use eval for double. 
                energy_index = fit_result_info["info"].energy_index
                if fit_result_info["info"].model.short_name!="1-exp":
                    x=np.linspace(fit_result_info["info"].tmin, fit_result_info["info"].tmax, fit_result_info["info"].tmax-fit_result_info["info"].tmin+1)
                    model = fit_result_info["info"].model.sigmond_object(Nt)
                    if fit_result_info["info"].sim_fit: #make distinct between Deg/2-3exponential
                        if sh_index==0:
                            min_index = 0
                            max_index = fit_result_info["info"].num_params
                            params = [estimate.getFullEstimate() for estimate in fit_result_info["estimates"][min_index:max_index]] 
                        if sh_index==1:
                            model=fit_info.FitModel.TimeForwardTwoExponential.sigmond_object(Nt)
                            min_index = fit_result_info["info"].num_params
                            indexes = [min_index,min_index+1,2,min_index+2]
                            params = [fit_result_info["estimates"][i].getFullEstimate() for i in indexes]
                        if sh_index==2:
                            model=fit_info.FitModel.TimeForwardTwoExponentialForCons.sigmond_object(Nt)
                            min_index = fit_result_info["info"].num_params+3
                            indexes = [min_index,min_index+1,2,3,min_index+2]
                            params = [fit_result_info["estimates"][i].getFullEstimate() for i in indexes]
                        energy_index = min_index #correspond to fit model
                    else:
                        params = []
                        for estimate in fit_result_info["estimates"]:
                            params.append(estimate.getFullEstimate())
                    y = [-np.log(model.eval(params, x[i+1])/model.eval(params, x[i])) for i in range(len(x)-1)]
                    plt.plot(x[:-1]+0.5,y, color="black", ls="--")
                energy_result = fit_result_info["estimates"][energy_index].getFullEstimate()
                energy_err = fit_result_info["estimates"][energy_index].getSymmetricError()
                plt.hlines(energy_result, fit_result_info["info"].tmin, fit_result_info["info"].tmax, color="black", zorder=2)
                plt.gca().add_patch(patches.Rectangle((fit_result_info["info"].tmin, energy_result-energy_err), fit_result_info["info"].tmax-fit_result_info["info"].tmin, 2.0*energy_err, zorder=0, color="gray"))
            if self.latex:
                if fit_result_info["info"].ratio:
                    yscale = r"$a_t \delta E_{\textup{lab}}$" #but the use of "\textup{}" command requires a latex compiler
                else:
                    yscale = r"$a_tE_{\textup{lab}}$" #but the use of "\textup{}" command requires a latex compiler
            else: 
                if fit_result_info["info"].ratio:
                    yscale = r"$a_t \delta E_{lab}$" #unsure if \delta needs latex
                else:
                    yscale = r"$a_tE_{lab}$"
                    
            plt.ylabel(yscale)
            if fit_result_info["success"]:
                nice_value = utils.nice_value(energy_result,energy_err)
                labels.append(f"{yscale}={nice_value}")

        plt.xlabel(r"$t/a_t$")

        # Annotate the plot with additional information if provided
        if op1:
            labels.insert(0,f"corr: {str(op1)}") #double check that I'm not messing up sink and source
        if fit_result_info["success"]:
            labels.append(f"$\chi^2/dof={format(float(fit_result_info['chisqrdof']),'.2f')}$")

        self.moving_textbox(labels)

        plt.tight_layout()

    def add_tmin_fits(self, t,e,de,color_index, filled=True, model=None): #incorporate pvalue and chosen fit
        if filled:
            marker_color = psettings.colors[color_index]
        else:
            marker_color = "white"
        plt.errorbar(x=t,y=e,yerr=de, linewidth=0.0, elinewidth=1.5, capsize=5, color=psettings.colors[color_index], 
                     marker=psettings.markers[color_index], label=model, mfc=marker_color)

    def add_chosen_fit(self, energyval, energyerr, label="chosen"):
        plt.axhspan(energyval-energyerr, energyval+energyerr, color="gray")
        plt.axhline(energyval,color="black",ls="--", label=label)

    def finalize_tmin_plot(self, title=None, ratio=False):
        if self.latex:
            if ratio:
                yscale = r"$a_t \delta E_{\textup{lab}}$" #but the use of "\textup{}" command requires a latex compiler
            else:
                yscale = r"$a_tE_{\textup{lab}}$" #but the use of "\textup{}" command requires a latex compiler
        else: 
            if ratio:
                yscale = r"$a_t \delta E_{lab}$" #unsure if \delta needs latex
            else:
                yscale = r"$a_tE_{lab}$"

        if self.latex:
            plt.xlabel(r"$t_{\textup{min}}/a_t$")
        else:
            plt.xlabel(r"$t_{min}/a_t$")
        plt.ylabel(yscale)
        plt.legend(title=title)
        plt.tight_layout()

    def summary_plot(self,indexes,levels,errs,xticks, reference=None, thresholds=[]):
        color_index = 0
        indexes = np.array(indexes)
        levels = np.array(levels)
        errs = np.array(errs)
        shifted_array = 0.25*shift_levels(indexes,levels,errs)
        plt.errorbar(x=indexes+shifted_array, y=levels, yerr=errs,linewidth=0.0, elinewidth=1.5, capsize=5, color=psettings.colors[color_index], marker=psettings.markers[color_index])
        plt.xlim(min(indexes)-1.0,max(indexes)+1.0)

        dd = 0.005
        for line in thresholds:
            plt.axhline(line[1], color="black", ls="--")
            line_label = ""
            for particle in line[0]:
                line_label+=psettings.latex_format[particle]
            minx,maxx = plt.xlim()
            miny,maxy = plt.ylim()
            plt.text(minx+dd*(maxx-minx), line[1]+dd*(maxy-miny), line_label )
        xticks = [f"{psettings.latex_format[irrep]}({mom})" for (irrep, mom) in xticks]
        plt.xticks(list(range(len(xticks))),xticks)
        if self.latex:
            if reference:
                latex_rest_mass = psettings.latex_format[reference].replace('$',"")
                plt.ylabel(rf"$E_{{\textup{{cm}}}}/E_{{{latex_rest_mass}}}$") #change to ref
            else:
                plt.ylabel(r"$a_t E_{\textup{cm}}$") #change to ref
        else:
            if reference:
                latex_rest_mass = psettings.latex_format[reference].replace('$',"")
                plt.ylabel(rf"$E_{{cm}}/E_{{{latex_rest_mass}}}$") #change to ref
            else:
                plt.ylabel(r"$a_t E_{cm}$")

    def plot_operator_overlaps(self, values, errors, opname):
        plt.bar(range(len(values)),values)
        plt.errorbar(x=range(len(values)), y=values, yerr=errors,linewidth=0.0, elinewidth=1.5, capsize=5, marker=None, color="black")
        plt.xlabel("rotate level")
        plt.ylabel("$|Z|^2$")
        self.moving_textbox([opname])

    #these pickle can only be opened on same system as it was generated 
    # (they are written with bytes because that was the only it would let me)
    def save_pickle( self, filename):
        """Save the current Matplotlib figure as a pickle file."""
        # fig = copy.deepcopy(self.fig)
        with open(filename, "wb") as f:
            pickle.dump(self.fig, f)

    def save_pdf( self, filename, transparent=True):
        """Save the current Matplotlib figure as a PDF file."""
        plt.savefig( filename, transparent=transparent ) 

    ###########################
    ##### pylatex actions #####
    ###########################
    
    def create_summary_doc(self, title):
        """Create a new LaTeX document with a specified title."""
        self.doc.append(utils.create_doc(title))

    def append_section(self, title, index = 0):
        """Append a new section to the LaTeX document."""
        self.doc[index].append(pylatex.Command("section",title))

    def append_subsection(self, title, index = 0):
        """Append a new section to the LaTeX document."""
        self.doc[index].append(pylatex.Command("subsection",title))

    def append_subsubsection(self, title, index = 0):
        """Append a new section to the LaTeX document."""
        self.doc[index].append(pylatex.Command("subsubsection",title))

    def add_correlator_subsection(self,corrname, leftplotfile, rightplotfile, index = 0): #add table of estimates?, list of correlators?
        """Add a subsection with two plots to the LaTeX document."""
        if not os.path.exists(leftplotfile) and not os.path.exists(rightplotfile):
            logging.warning(f"Unable to include {corrname} in summary pdf.")
            return

        with self.doc[index].create(pylatex.Subsubsection(corrname)):
            self.include_additional_plots(leftplotfile,rightplotfile, index)
    
    def include_additional_plots(self, leftplotfile, rightplotfile, index = 0): #add table of estimates?, list of correlators?
        """Add a subsection with two plots to the LaTeX document."""

        rlinewidth = r'\linewidth'
        llinewidth = r'\linewidth'
        if not os.path.exists(leftplotfile) and not os.path.exists(rightplotfile):
            return

        with self.doc[index].create(pylatex.Figure(position='H')) as thisfig:
            if os.path.exists(leftplotfile) and os.path.exists(rightplotfile):
                with self.doc[index].create(pylatex.SubFigure(
                    position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as left_fig:
                    left_fig.add_image(leftplotfile,width=pylatex.NoEscape(llinewidth), placement=pylatex.NoEscape("\centering") )
                with self.doc[index].create(pylatex.SubFigure(
                    position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as right_fig:
                    right_fig.add_image(rightplotfile,width=pylatex.NoEscape(rlinewidth), placement=pylatex.NoEscape("\centering") )
            elif os.path.exists(leftplotfile):
                thisfig.add_image(leftplotfile,width=pylatex.NoEscape(r'0.5\linewidth'), placement=pylatex.NoEscape("\centering") )
            elif os.path.exists(rightplotfile):
                thisfig.add_image(rightplotfile,width=pylatex.NoEscape(r'0.5\linewidth'), placement=pylatex.NoEscape("\centering") )

    def add_single_plot(self, plotfile, index=0):
        if not os.path.exists(plotfile):
            logging.warning(f"Unable to include {plotfile} in summary pdf.")
            return
        
        with self.doc[index].create(pylatex.Figure(position='H')) as thisfig:
            thisfig.add_image(plotfile,width=pylatex.NoEscape(r'\linewidth'), placement=pylatex.NoEscape("\centering") )


    def summary_table(self, reference, headers, data, title = "", index = 0):
        if reference:
            latex_rest_mass = psettings.latex_format[reference].replace('$',"")
            headers = [header.replace("latex_rest_mass",latex_rest_mass) for header in headers]
        headers = [pylatex.NoEscape(header) for header in headers]
        with self.doc[index].create(pylatex.Center()) as centered:
            # if title:
                # self.doc[index].append(pylatex.Command("caption",pylatex.NoEscape(pylatex.utils.bold(pylatex.NoEscape(title)))))
                # self.doc[index].append(pylatex.NoEscape("\n\n"))
            with self.doc[index].create(pylatex.Table(position='h!')) as table:
                table.append(pylatex.Command("centering"))
                if title:
                    table.add_caption(pylatex.NoEscape(pylatex.utils.bold(pylatex.NoEscape(title))))
                with table.create(pylatex.Tabular("|".join(['c']*len(headers)))) as current_table:
                    current_table.add_row(headers)
                    current_table.add_hline()
                    for line in data:
                        line = [psettings.latex_format[col] if col in psettings.latex_format.keys() else col for col in line]
                        line = [pylatex.NoEscape(col) for col in line]
                        current_table.add_row(line)

    def add_operator_overlaps(self, files, index=0):
        with self.doc[index].create(pylatex.Subsubsection("Operator Overlaps")):
            self.add_plot_series(files,index)
            # for x,y in zip(files[::2],files[1::2]):
            #     self.include_additional_plots(x,y, index)
            # if len(files)%2:
            #     self.include_additional_plots(files[-1],files[-1]+"2")

    def add_plot_series(self, files, index=0):
        for x,y in zip(files[::2],files[1::2]):
            self.include_additional_plots(x,y, index)
        if len(files)%2:
            self.include_additional_plots(files[-1],files[-1]+"2")

    def compile_pdf(self, filename, index = 0):
        """Compile the LaTeX document into a PDF file."""
        utils.compile_pdf(self.doc[index],filename) #detect compiler? -> if self.latex, check for pylatex?