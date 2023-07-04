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
parser.add_argument("--plot_par",   default='muRV',          help='Set parameter to plot')
parser.add_argument("--quantilemode",	default=False,	  help='If True, annotate with 16,84, False, use sample std.')
parser.add_argument("--save",	default=True,		      help='Save plot')
parser.add_argument("--show",	default=False,		      help='Show plot')

args = parser.parse_args().__dict__
plot_par  = args['plot_par']
parlabels = {'muRV':'\\mu_{R_V}','tauA':'\\tau_{A}'}
dfpars    = {'muRV':'mu_RV','tauA':'tauA'}

Nsim_keep = 100
Rhat_threshold = 1.1
#Choices for applying HBM to simulated data
BIRDSNACK_EDIT_DICT1 = {'analysis_parameters':
						{'HBM_savekey':'PPC_muRV2.5sigRV0.1tauA0.5_LowBVNoCens',
						'CensoredData':True,'CensoredCut':0.3}}
BIRDSNACK_EDIT_DICT2 = {'analysis_parameters':
						{'HBM_savekey':'PPC_muRV2.5sigRV0.1tauA0.5_WithCensoredData',
						'CensoredData':True,'CensoredCut':'inf'}}

#edit_dict = {'simulate_parameters':{'S':100,'tauA':0.5,'muRV':2.5,'sigRV':0.5,'PredefinedExtrinsicHyps':False,
#			'pre_defined_hyps':{'load_file':'AVExp_Cens1.0'}}}
edit_dict = {'simulate_parameters':{'S':100,'tauA':0.5,'muRV':2.5,'sigRV':0.1,'PredefinedExtrinsicHyps':False,
			'pre_defined_hyps':{'load_file':'AVExp_Cens1.0'}}}

BS_editdicts = dict(zip([False,True],[BIRDSNACK_EDIT_DICT1,BIRDSNACK_EDIT_DICT2]))

if __name__ == "__main__":
	print (f"Plotting PPC of Censored Data for {plot_par};")
	with open('ppc.yaml') as f:
		sbc_choices = yaml.load(f, Loader=yaml.FullLoader)

	GLOB_FITS = {}
	for include_cens in BS_editdicts:
		#Upload these fits
		edit_dict  = update_edit_dict_for_ppc(sbc_choices,edit_dict)
		#Get SBC_CLASS
		sbc  = SBC_CLASS(sbc_choices,edit_dict)
		#Get FITS
		FITS = sbc.get_fits(edit_dict=BS_editdicts[include_cens])
		#Trim FITS to remove large Rhat
		FITS = trim_to_KEEPERS({'dummy':FITS},get_KEEPERS({'dummy':FITS},Nsim_keep,Rhat_threshold,plot_par,dfpars[plot_par]))['dummy']
		#Store
		GLOB_FITS[include_cens] = FITS

	#True parameter value
	true_plot_par = sbc.__dict__[plot_par]

	#Plot SBC
	fig,axs = pl.subplots(len(GLOB_FITS),1,figsize=(8,12),sharex=True)
	for iax,include_cens in enumerate(GLOB_FITS):
		Lside   = False#True if include_cens and plot_par=='tauA' else False
		FITS    = GLOB_FITS[include_cens]
		plotter = SBC_FITS_PLOTTER(iax,fig.axes,[true_plot_par,plot_par,dfpars[plot_par],parlabels[plot_par]],FITS,sbc.bs.choices['analysis_parameters'],sbc.path_to_birdsnack_rootpath,quantilemode=args['quantilemode'])
		plotter.plot_sbc_panel(annotate_true=False,color={0:'C0',1:'indigo'}[iax],Lside=Lside,FAC=int(Nsim_keep/8))
		if include_cens:
			fig.axes[iax].annotate(f'Include Censored Data',xy=(0.95-(0.95-0.0225)*Lside,0.5+0.02),xycoords='axes fraction',fontsize=plotter.FS,ha={True:'left',False:'right'}[Lside],weight='bold')


	fig.axes[-1].set_xlabel(r'$%s$'%parlabels[plot_par]+' (mag)'*(plot_par=='tauA'),fontsize=plotter.FS)
	fig.axes[0].set_ylabel('Posterior Densities',fontsize=plotter.FS,color='white')#For spacing
	fig.text(0.01, 0.5, 'Posterior Densities', ha='center', va='center', rotation='vertical',fontsize=plotter.FS)
	fig.axes[0].set_title(r'Fits to Simulated Data with $|B-V|<0.3\,$mag Cut Applied'+'\n'+r'True Simulation Parameters: $\mu_{R_V}=%s$ ; $\sigma_{R_V} = %s$ ; $\tau_A = %s\,$mag'%(round(sbc.muRV,2),round(sbc.sigRV,2),round(sbc.tauA,2)),fontsize=plotter.FS)
	pl.tight_layout()
	if args['save']:
		pl.savefig(f"{sbc.plotpath}PPCCensoredData_plotpar{plot_par}_quantilemode{args['quantilemode']}.pdf",bbox_inches='tight')
	if args['show']:
		pl.show()
