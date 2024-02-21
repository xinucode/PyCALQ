
import matplotlib.pyplot as plt
import pickle
import pylatex
import copy
import os
import logging

import fvspectrum.sigmond_util as sigmond_util
import fvspectrum.spectrum_plotting_settings.settings as psettings
import fvspectrum.sigmond_utils.util as utils


ctext_x = 0.3
ctext_y = 0.7

class PlottingHandler:

    doc = []

    def __init__(self):

        # Flag to check if LaTeX is available for plot labels
        self.latex = True
        #checks if latex compler is present, if not sets this variable to false
            #use this for latex labels, and make an option for non-latex labels.
        self.latex = sigmond_util.set_latex_in_plots(plt.style)

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

    def correlator_plot(self,df, ptype=0, op1=None, op2=None, color_index = 0):
        """Generate a correlator plot using Matplotlib."""
        plt.errorbar(x=df["aTime"],y=df["FullEstimate"],yerr=df["SymmetricError"], linewidth=0.0, elinewidth=1.5, capsize=5, color=psettings.colors[color_index], marker=psettings.markers[color_index] )
        if ptype==0:
            plt.ylabel(r"$C(t)$") #not all dollar signs require a latex compiler, for example this is okay.
        else:
            if self.latex:
                plt.ylabel(r"$a_tE_{\textup{fit}}$") #but the use of "\textup{}" command requires a latex compiler
            else: 
                plt.ylabel(r"$a_tE_{fit}$")

        plt.xlabel(r"$t/a_t$")

        # Annotate the plot with additional information if provided
        if op1:
            plt.figtext(ctext_x,ctext_y+0.1,f"snk: {str(op1)}") #double check that I'm not messing up sink and source
        if op2:
            plt.figtext(ctext_x,ctext_y,f"src: {str(op2)}")

        plt.tight_layout()

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

    def add_correlator_subsection(self,corrname, leftplotfile, rightplotfile, index = 0): #add table of estimates?, list of correlators?
        """Add a subsection with two plots to the LaTeX document."""

        rlinewidth = r'\linewidth'
        llinewidth = r'\linewidth'
        if not os.path.exists(leftplotfile):
            leftplotfile = rightplotfile
            llinewidth = r'0.1\linewidth'
        if not os.path.exists(rightplotfile):
            rightplotfile = leftplotfile
            rlinewidth = r'0.1\linewidth'
        if not os.path.exists(leftplotfile) and not os.path.exists(rightplotfile):
            logging.warning(f"Unable to include {corrname} in summary pdf.")
            return

        with self.doc[index].create(pylatex.Subsection(corrname)):
            with self.doc[index].create(pylatex.Figure(position='H')) as thisfig:
                with self.doc[index].create(pylatex.SubFigure(
                    position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as left_fig:
                    left_fig.add_image(leftplotfile,width=pylatex.NoEscape(llinewidth), placement=pylatex.NoEscape("\centering") )
                with self.doc[index].create(pylatex.SubFigure(
                    position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as right_fig:
                    right_fig.add_image(rightplotfile,width=pylatex.NoEscape(rlinewidth), placement=pylatex.NoEscape("\centering") )
    
    def compile_pdf(self, filename, index = 0):
        """Compile the LaTeX document into a PDF file."""
        utils.compile_pdf(self.doc[index],filename) #detect compiler? -> if self.latex, check for pylatex?