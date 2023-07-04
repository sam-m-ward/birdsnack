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
parser.add_argument("--loop_pars",     default='muRV,sigRV',  help='Parameters to inspect on each panel')
parser.add_argument("--quantilemode",	default=True,	  help='If True, annotate with 16,84, False, use sample std.')
parser.add_argument("--save",	default=True,		      help='Save plot')
parser.add_argument("--show",	default=False,		      help='Show plot')

args = parser.parse_args().__dict__
loop_pars = [s for s in args["loop_pars"].split(',')]
parnames,dfpars,parlabels = get_pars()

Nsim_keep = 5
Rhat_threshold = 1.05
rec_file = 'AVExp_Cens1.0'#FILE USED FOR SIMULATING; FIT TO REAL DATA
#Choices for applying HBM to simulated data
BIRDSNACK_EDIT_DICT = {'analysis_parameters':
						{'HBM_savekey':'PPC_LowBVwCens_DevFitDev',
						'CensoredData':True}}

edit_dict = {'simulate_parameters':{'S':250,'pre_defined_hyps':{'load_file':'AVExp_Cens1.0'}}}

if __name__ == "__main__":
	print (f"Plotting PPC of {loop_pars};")
	with open('ppc.yaml') as f:
		sbc_choices = yaml.load(f, Loader=yaml.FullLoader)

	#Upload these fits
	edit_dict  = update_edit_dict_for_ppc(sbc_choices,edit_dict)
	#Get SBC_CLASS
	sbc  = SBC_CLASS(sbc_choices,edit_dict)
	#Get FITS
	FITS = sbc.get_fits(edit_dict=BIRDSNACK_EDIT_DICT)
	#Trim FITS to remove large Rhat
	for par in loop_pars:
		FITS = trim_to_KEEPERS({'dummy':FITS},get_KEEPERS({'dummy':FITS},Nsim_keep,Rhat_threshold,par,dfpars[par]))['dummy']

	#Real Data Fit Samples
	with open(sbc.bs.FITSpath+f"FIT{rec_file}.pkl",'rb') as f:
		rec_samps = pickle.load(f)['df']

	#Plot SBC
	fig,axs = pl.subplots(len(loop_pars),1,figsize=(8,6*len(loop_pars)))
	for iax,par in enumerate(loop_pars):
		plotter = SBC_FITS_PLOTTER(iax,fig.axes,[sbc.__dict__[par],par,dfpars[par],parlabels[par]],FITS,sbc.bs.choices['analysis_parameters'],sbc.path_to_birdsnack_rootpath,quantilemode=args['quantilemode'])
		plotter.plot_sbc_panel(Ncred=False,Parcred=False,annotate_true=False,plot_true=False,include_pmedian=False,real_data_samps=rec_samps[dfpars[par]])
		fig.axes[iax].set_xlabel(r'$%s$'%parlabels[par],fontsize=plotter.FS)
		fig.axes[iax].plot(sbc.__dict__[par]*np.ones(2),[0,fig.axes[iax].get_ylim()[1]],c='C3',linewidth=2,linestyle='-')

	fig.axes[0].set_ylabel('Posterior Densities',fontsize=plotter.FS,color='white')#For spacing
	fig.text(0.01, 0.5, 'Posterior Densities', ha='center', va='center', rotation='vertical',fontsize=plotter.FS)
	fig.axes[0].set_title(r'Recovery of $R_V$ Hyperparameters from fits to %s Simulated SNe'%sbc.S+'\n'+ \
			r'True Simulation Parameters: $\mu_{R_V}=%s$ ; $\sigma_{R_V} = %s$ ; $\tau_A = %s\,$mag'%(round(sbc.muRV,2),round(sbc.sigRV,2),round(sbc.tauA,2)),fontsize=plotter.FS)
	pl.tight_layout()
	if args['save']:
		pl.savefig(f"{sbc.plotpath}PPC_looppars{loop_pars}_quantilemode{args['quantilemode']}_HBM{sbc.bs.choices['analysis_parameters']['HBM_savekey']}.pdf",bbox_inches='tight')
	if args['show']:
		pl.show()
