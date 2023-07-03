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
parser.add_argument("--plot_par",   default='nu',          help='Set parameter to plot')
parser.add_argument("--loop_par",   default='S',           help='Set parameter to change in each panel')
parser.add_argument("--loop_S",     default='65,125,250',  help='NSNe on each panel')
parser.add_argument("--save",	default=True,		      help='Save plot')
parser.add_argument("--show",	default=False,		      help='Show plot')
parser.add_argument("--quantilemode",	default=True,	  help='If True, annotate with 16,84, False, use sample std.')

args = parser.parse_args().__dict__
plot_par = args['plot_par']
true_plot_par = 1
loop_par = args['loop_par']
loop_par_dict = {loop_par:[int(s) for s in args["loop_S"].split(',')]}
parlabels = {'nu':'\\nu_{A}'}
dfpars    = {'nu':'nu'}

Nsim_keep = 5
Rhat_threshold = 1.1
load_file = 'Fiducial_CensoredCut1.0'#FILE USED FOR SIMULATING; FIT TO REAL DATA
rec_file  = 'Fiducial_AVGamma_CensoredCut1.0'#FILE USED FOR RECOVERY; FIT TO REAL DATA
#Choices for applying HBM to simulated data
BIRDSNACK_EDIT_DICT = {'analysis_parameters':
						{'HBM_savekey':'PPC_CensoredCut1.0_ExpFitGamma',
						'CensoredData':True,'CensoredCut':1.0,
						'AVprior':'Gamma','n_warmup':1000,'n_sampling':2000,'n_thin':1000}}

if __name__ == "__main__":
	print (f"Plotting PPC of {loop_par_dict};")
	with open('ppc.yaml') as f:
		sbc_choices = yaml.load(f, Loader=yaml.FullLoader)

	GLOB_FITS = {}
	for parvalue in loop_par_dict[loop_par]:
		#Upload these fits
		additional = {loop_par:parvalue,'pre_defined_hyps':{'load_file':load_file}}
		edit_dict  = update_edit_dict_for_ppc(sbc_choices,{'simulate_parameters':additional})
		#Get SBC_CLASS
		sbc  = SBC_CLASS(sbc_choices,edit_dict)
		#Get FITS
		FITS = sbc.get_fits(edit_dict=BIRDSNACK_EDIT_DICT)
		#Trim FITS to remove large Rhat
		FITS = trim_to_KEEPERS({'dummy':FITS},get_KEEPERS({'dummy':FITS},Nsim_keep,Rhat_threshold,plot_par,dfpars[plot_par]))['dummy']
		#Store
		GLOB_FITS[parvalue] = FITS

	#Real Data Fit Samples
	with open(sbc.bs.FITSpath+f"FIT{rec_file}.pkl",'rb') as f:
		rec_samps = pickle.load(f)['df'][dfpars[plot_par]]

	#Plot SBC
	fig,axs = pl.subplots(len(GLOB_FITS),1,figsize=(8,14),sharex=True,gridspec_kw={'height_ratios': [1.5, 1, 1]})
	for iax,true_loop_par in enumerate(GLOB_FITS):
		FITS    = GLOB_FITS[true_loop_par]
		plotter = SBC_FITS_PLOTTER(iax,fig.axes,[true_plot_par,plot_par,dfpars[plot_par],parlabels[plot_par]],FITS,sbc.bs.choices['analysis_parameters'],sbc.path_to_birdsnack_rootpath,quantilemode=args['quantilemode'])
		plotter.plot_sbc_panel(Ncred=False,Parcred=True,annotate_true=False,real_data_samps=rec_samps if iax==0 else False)
		fig.axes[iax].annotate(f'No. of Simulated SNe = {true_loop_par}',	xy=(0.95,0.475+0.02-0.05),xycoords='axes fraction',fontsize=plotter.FS,ha='right',weight='bold')

	fig.axes[-1].set_xlabel(r'$%s$'%parlabels[plot_par],fontsize=plotter.FS)
	fig.axes[0].set_ylabel('Posterior Densities',fontsize=plotter.FS,color='white')#For spacing
	fig.text(0.01, 0.5, 'Posterior Densities', ha='center', va='center', rotation='vertical',fontsize=plotter.FS)
	fig.axes[0].set_title(r'Free-$\nu_A$ Fits to Data Simulated with $\nu_A=1$'+'\n'+r'True Simulation Parameters: $\mu_{R_V}=%s$ ; $\sigma_{R_V} = %s$ ; $\tau_A = %s\,$mag'%(round(sbc.muRV,2),round(sbc.sigRV,2),round(sbc.tauA,2)),fontsize=plotter.FS)
	pl.tight_layout()
	if args['save']:
		pl.savefig(f"{sbc.plotpath}PPC_plotpar{plot_par}_looppars{loop_par_dict}_quantilemode{args['quantilemode']}.pdf",bbox_inches='tight')
	if args['show']:
		pl.show()
