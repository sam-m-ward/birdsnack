"""
Plot Simulation-Based Calibration
--------------------
Simple script that uses SBC fits to:
	plot up recovery of input dust hyperparameters
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

import sys
sys.path.append('model_files/')
from SBC import *
from sbc_plot_functions import *
import argparse,yaml
import matplotlib.pyplot as pl

parser = argparse.ArgumentParser(description="SBC Input Dust Hyperparameters")
parser.add_argument("--loop_par",	default='muRV',		  help='Set parameter to loop over')
parser.add_argument("--muRV",		default=2.5, 	      help='Float for non-loop RV Pop. Mean')
parser.add_argument("--sigRV",		default=0.5,		  help='Float for non-loop RV Gaussian Dist.')
parser.add_argument("--tauA",		default=0.5, 	      help='Float for non-loop AV Dust Hyp')
parser.add_argument("--loop_muRV",  default='1.5,2.5,3.5',help='Parameter values to loop over')
parser.add_argument("--loop_sigRV", default='0.1,0.5,1.0',help='Parameter values to loop over')
parser.add_argument("--loop_tauA",  default='0.2,0.5',    help='Parameter values to loop over')

args = parser.parse_args().__dict__
parnames      = ['muRV','sigRV','tauA']
parlabels     = dict(zip(parnames,['\mu_{R_V}','\sigma_{R_V}','\\tau_A']))
dfpars        = dict(zip(parnames,['mu_RV','sig_RV','tauA']))
loop_par      = args['loop_par']
loop_par_dict = {loop_par:[float(s) for s in args[f"loop_{args['loop_par']}"].split(',')]}
non_loop_pars = {par:args[par] for par in parnames if par!=args['loop_par']}

Nsim_keep = 20
Rhat_threshold = 1.05

if __name__ == "__main__":
	print (f"Plotting SBC using input dust hyperparameters: {loop_par_dict};{non_loop_pars}")
	with open('sbc.yaml') as f:
		sbc_choices = yaml.load(f, Loader=yaml.FullLoader)

	GLOB_FITS = {}
	for parvalue in loop_par_dict[loop_par]:
		#Current Dust Hyps
		dust_hyps = {**non_loop_pars,**{loop_par:parvalue}}
		#Load fits with dust hyperparameters set to these values
		edit_dict = {'simulate_parameters':dust_hyps}
		#Get SBC_CLASS
		sbc  = SBC_CLASS(sbc_choices,edit_dict)
		#Get FITS
		FITS = sbc.get_fits()
		#Trim FITS to remove large Rhat
		FITS = trim_to_KEEPERS({parvalue:FITS},get_KEEPERS({parvalue:FITS},Nsim_keep,Rhat_threshold,loop_par,dfpars[loop_par]))[parvalue]
		#Store
		GLOB_FITS[parvalue] = FITS

	#Plot SBC
	fig,axs = pl.subplots(len(GLOB_FITS),1,figsize=(8,6*len(GLOB_FITS)),sharex=False)
	ax      = fig.axes
	for iax,true_par in enumerate(GLOB_FITS):
		FITS    = GLOB_FITS[true_par]
		plotter = SBC_FITS_PLOTTER(iax,ax,[true_par,loop_par,dfpars[loop_par],parlabels[loop_par]],FITS,sbc.bs.choices['analysis_parameters'],sbc.path_to_birdsnack_rootpath)
		plotter.plot_sbc_panel()
	title_tuples = [l for key,value in non_loop_pars.items() for l in [parlabels[key],value]]
	ax[0].set_title('Fits to SED-Integrated Simulated Data;\nTrue Simulation Parameters: '+r'$%s = %s ; %s = %s\,$mag'%(title_tuples[0],title_tuples[1],title_tuples[2],title_tuples[3]),fontsize=plotter.FS)
	ax[-1].set_xlabel(r'$%s$'%parlabels[loop_par],fontsize=plotter.FS)
	fig.text(0.01, 0.5, 'Posterior Densities', ha='center', va='center', rotation='vertical',fontsize=plotter.FS)
	pl.tight_layout()
	pl.show()
