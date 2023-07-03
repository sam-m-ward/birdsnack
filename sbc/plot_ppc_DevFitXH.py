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
parlabels = {'muRV':'\\mu_{R_V}'}
dfpars    = {'muRV':'mu_RV'}

Nsim_keep = 10
Rhat_threshold = 1.1

#For X-H Sims (X-H model applied to Deviations-simulated SNe)
BIRDSNACK_EDIT_DICT = {'analysis_parameters':
						{'HBM_savekey':'PPC_CensoredCut1.0_DevFitXH',
						'CensoredData':True,'CensoredCut':1.0,
						'IntrinsicModel':'X-H','DataTransformation' : 'X-H'}}

#Choices for simulating data based on previous stan fit with BirdSnack
edit_dict1 = {'simulate_parameters':{'Nsims':20,'S':250,'pre_defined_hyps':{'load_file':'Fiducial_CensoredCut1.0'}}}
edit_dict2 = {'simulate_parameters':{'Nsims':20,        'pre_defined_hyps':{'load_file':'Fiducial_CensoredCut1.0'}}}
editdicts = dict(zip([65,250],[edit_dict1,edit_dict2]))

if __name__ == "__main__":
	print (f"Plotting PPC of X-H Model applied to Deviations-Simulated Data for {plot_par};")
	with open('ppc.yaml') as f:
		sbc_choices = yaml.load(f, Loader=yaml.FullLoader)

	GLOB_FITS = {}
	for S in editdicts:
		#Upload these fits
		edit_dict  = update_edit_dict_for_ppc(sbc_choices,editdicts[S])
		#Get SBC_CLASS
		sbc  = SBC_CLASS(sbc_choices,edit_dict)
		#Get FITS
		FITS = sbc.get_fits(edit_dict=BIRDSNACK_EDIT_DICT)
		#Trim FITS to remove large Rhat
		FITS = trim_to_KEEPERS({'dummy':FITS},get_KEEPERS({'dummy':FITS},Nsim_keep,Rhat_threshold,plot_par,dfpars[plot_par]))['dummy']
		#Store
		GLOB_FITS[S] = FITS

	#True parameter value
	true_plot_par = sbc.__dict__[plot_par]

	#Plot SBC
	pl.figure()
	for iim,S in enumerate(GLOB_FITS):
		FITS    = GLOB_FITS[S]
		plotter = SBC_FITS_PLOTTER(0,[pl.gca()],[true_plot_par,plot_par,dfpars[plot_par],parlabels[plot_par]],FITS,sbc.bs.choices['analysis_parameters'],sbc.path_to_birdsnack_rootpath,quantilemode=args['quantilemode'])
		plotter.plot_sbc_panel(Ncred=False,annotate_true=False,plot_ind=False,plot_true=False,plot_medians=False,dress_figure=False,fill_between=False,color=f"C{1-iim}",linestyle=['--','-'][iim],line_sap_title='Sim. Posterior')

	#Real Data Fit Samples
	XH_noLCfile = f"XHCols_Cens1.0"
	XH_wLCfile  = f"XHCols_wLCshape_Cens1.0"
	with open(sbc.bs.FITSpath+f"FIT{XH_noLCfile}.pkl",'rb') as f:	XH_noLC = pickle.load(f)['df'][dfpars[plot_par]]
	with open(sbc.bs.FITSpath+f"FIT{XH_wLCfile}.pkl",'rb') as f:	XH_wLC  = pickle.load(f)['df'][dfpars[plot_par]]
	#Plot real data fits
	import sys
	sys.path.append(self.path_to_birdsnack_rootpath+'model_files/')
	from posterior_plotter import PARAMETER
	samps = PARAMETER(XH_noLC,dfpars[plot_par],parlabels[plot_par],plotter.lims[plot_par],plotter.bounds[plot_par],0,0,{})
	samps.get_xgrid_KDE()
	pl.plot(samps.xgrid,samps.KDE,alpha=0.5,linewidth=3,color='C0',linestyle='-.',
				label=r"Real Posterior"+'\n'+r"(w/ $\Delta m_{15}(B)$ term)"+'\n'+r'$%s = %s \pm %s$'%(parlabels[par],XH_wLC.quantile(0.5).round(2),XH_wLC.std().round(2)))
	samps = PARAMETER(XH_wLC,dfpars[plot_par],parlabels[plot_par],plotter.lims[plot_par],plotter.bounds[plot_par],0,0,{})
	samps.get_xgrid_KDE()
	pl.plot(samps.xgrid,samps.KDE,alpha=0.5,linewidth=3,color='C0',linestyle=':',
				label=r"Real Posterior"+ '\n'+r'$%s = %s \pm %s$'%(parlabels[par],XH_noLC.quantile(0.5).round(2),XH_noLC.std().round(2)))

	#Finish figure
	for iim,S in enumerate(editdicts):
		pl.annotate(r'$N_{SNe}=%s$'%S,xy=(1.75,pl.gca().get_ylim()[1]*(0.825-0.1*iim)),color=f'C{1-iim}',fontsize=plotter.FS)
	pl.title(r'Recovery of $%s=%s$ under'%(parlabels[plot_par],round(true_plot_par,2))+'\n' r'$X-H$ Intrinsic Colour Model',fontsize=plotter.FS)
	pl.legend(fontsize=plotter.FS-5)
	pl.xlim([1.7,5])
	pl.ylim([0,None])
	pl.plot([true_plot_par,true_plot_par],[0,pl.gca().get_ylim()[1]],color='black')
	pl.yticks([])
	pl.xlabel(r"$%s$"%parlabels[plot_par],fontsize=plotter.FS)
	pl.ylabel('Simulation-Averaged Posterior',fontsize=plotter.FS-1)
	pl.tick_params(labelsize=plotter.FS)
	pl.tight_layout()
	if args['save']:
		pl.savefig(f"{sbc.plotpath}XHmodelrecovery_plotpar{plot_par}_quantilemode{args['quantilemode']}.pdf",bbox_inches='tight')
	if args['show']:
		pl.show()
