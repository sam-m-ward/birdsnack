"""
Plotting Functions Module

Module containing PLOTTER class
Helpful functions for plotting, contains arguments such as fontsize, choice to show/save plot etc.


Contains:
--------------------
PLOTTER class:
	inputs: plotting_parameters, plotpath

	Methods are:
		finish_plot(,xlab,ylab,savename=None,legend=False,invert=False)
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

from miscellaneous import ensure_folders_to_file_exist
import matplotlib.pyplot as pl
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

class PLOTTER:

	def __init__(self, plotting_parameters, plotpath):
		self.choices  = plotting_parameters
		self.plotpath = plotpath

		for key,value in self.choices.items():
			setattr(self,key,value)

	def finish_plot(self,xlab,ylab,savename=None,legend=False,invert=False):
		"""
		Finish Plot Function

		A generic function to finish off a plot

		Parameters
		----------
		xlab: string
			the x-axis label

		ylab: string
			the y-axis label

		savename : str


		invert: bool (optional; default=False)
			if True, flip the y-axis

		Returns
		----------
			Either saves plot, or shows plot, or both
		"""
		pl.ylabel(ylab,fontsize=self.FS)
		pl.xlabel(xlab,fontsize=self.FS)
		if legend:
			pl.legend(fontsize=self.FS)
		pl.tick_params(labelsize=self.FS)
		if invert:
			pl.gca().invert_yaxis()
		pl.tight_layout()
		if self.save and savename is not None:
			ensure_folders_to_file_exist(savename)
			pl.savefig(savename,bbox_inches='tight')
		if self.show:
			pl.show()
