"""
Simulation-Based Calibration Plotting Functions

Module contains functions useful for plotting SBC

Contains
----------
SBC_FITS_PLOTTER class
	inputs : self,iax,ax,PAR,FITS,analysis_choices,path_to_birdsnack_rootpath,quantilemode=True,Quantiles=[0,0.025,0.05,0.16,0.5,0.84,0.95,0.975,1],FS=18

	Methods are:
		get_SAMPS()
		get_QUANTILES()
		plot_sbc_panel(Ncred=True,Parcred=False,annotate_true=True,real_data_samps=False,plot_ind=True,plot_true=True,plot_medians=True,include_pmedian=True,dress_figure=True,fill_between=True,color='C0',linestyle='-',Lside=False,FAC=None,line_sap_title=None,line_rec_title=None)

Functions are:
	get_KEEPERS(GLOB_FITS,Nsim_keep,Rhat_threshold,loop_par,dfpar)
	trim_to_KEEPERS(GLOB_FITS,KEEPERS)
	update_edit_dict_for_ppc(sbc_choices, edit_dict)

--------------------

Written by Sam M. Ward: smw92@cam.ac.uk

"""
import pandas as pd
import numpy as np
from contextlib import suppress
import pickle

class SBC_FITS_PLOTTER:

	def __init__(self,iax,ax,PAR,FITS,analysis_choices,path_to_birdsnack_rootpath,quantilemode=True,Quantiles=[0,0.025,0.05,0.16,0.5,0.84,0.95,0.975,1],FS=18):
		self.iax = iax
		self.ax  = ax
		true_par,par,dfpar,parlabel = PAR[:]
		self.true_par = true_par
		self.par      = par
		self.dfpar    = dfpar
		self.parlabel = parlabel
		self.FITS                       = FITS
		self.analysis_choices           = analysis_choices
		self.path_to_birdsnack_rootpath = path_to_birdsnack_rootpath
		self.quantilemode = quantilemode
		self.Quantiles    = Quantiles
		self.FS = 18


	def get_SAMPS(self):
		"""
		Get Samples

		Method to get pandas df of samples across simulations

		Returns
		----------
		SAMPS : pandas df
			each row is one simulation; samples for the (hyper)parameter
		"""
		SAMPS = pd.DataFrame(data={ISIM:self.FITS[ISIM]['df'][self.dfpar].values for ISIM in self.FITS})
		return SAMPS

	def get_QUANTILES(self):
		"""
		Get Quantiles

		Method to get pandas df of quantiles across simulations

		Returns
		----------
		QUANTILES : pandas df
			quantiles across simulations
		"""
		QUANTILES = pd.DataFrame(data={q:[self.FITS[ISIM]['df'][self.dfpar].quantile(q) for ISIM in self.FITS] for q in self.Quantiles})
		return QUANTILES

	def plot_sbc_panel(self,Ncred=True,Parcred=False,annotate_true=True,real_data_samps=False,plot_ind=True,plot_true=True,plot_medians=True,include_pmedian=True,dress_figure=True,fill_between=True,color='C0',linestyle='-',Lside=False,FAC=None,line_sap_title=None,line_rec_title=None):
		"""
		Plot SBC Panel

		Method to plot a single panel of hyperparameter recovery

		Parameters
		----------
		Ncred : bool (optional; default=True)
			if True, plot N68, N95
		Parcred : bool (optional; default=False)
			if True, plot 68 and 95 quantiles of parameter
		annotate_true : bool (optional; default=True)
			if True, plot 'True par = value'
		real_data_samps : array (optional; default is None)
			if True, plot up samples from real-data fit
		plot_ind : bool (optional; default=True)
			if True, plot faint line posteriors, one for each sim
		plot_true: bool (optional; default=True)
			if True, plot line for True parameter value
		plot_medians : bool (optional; default=True)
			plot medians and include in legend
		include_pmedian : bool (optional; default=True)
			if True include p(median<Truth) line in legend
		dress_figure : bool (optional; default=True)
			apply e.g. lims, set yticks etc.
		fill_between : bool (optional; default=True)
			if True, include shaded regions for 68% credible interval in simulation averaged posterior
		color : str
			color of panel
		linestyle : str
			linestyle for simulation-averaged posterior
		Lside : bool (optional; default=False)
			if True, put annotations on LHS of panel
		FAC : float (optional; default=None)
			factor to reduce KDE grid for simulation-averaed posterior by compared to No.of samples
		line_sap_title : str (optional; default=None)
			string used in legend for simulation-averaged posterior, defaults to 'Simulation-Averaged Posterior'
		line_rec_title : str (optional; default=None)
			string used in legend for real-data posterior, defaults to 'Real-Data Posterior'

		End Product(s)
		----------
		ax[iax] panel with:
			faint lines for per-Sim posterior
			thick line for simulation-averaged posterior
			orange band for simulation medians
			legend with median and sap summaries
			line for true parameter and annotation
			N68,N95 annotations
		"""
		#Imports from BirdSnack
		import sys
		sys.path.append(self.path_to_birdsnack_rootpath+'model_files/')
		from posterior_plotter import PARAMETER
		from plotting_functions import get_parlabels
		analysis_choices = self.analysis_choices

		true_par  = self.true_par
		loop_par  = self.par
		dfpar     = self.dfpar
		parlabel  = self.parlabel
		FITS      = self.FITS
		Quantiles = self.Quantiles
		iax,ax    = self.iax,self.ax
		FS = self.FS
		HA = {True:'left',False:'right'}[Lside]

		#Get Parameter Bounds
		pars,dfparnames,parlabels,bounds  = get_parlabels({'analysis_parameters':analysis_choices})
		bounds = {loop_par:bounds[dfparnames.index(dfpar)]}

		#Get SAMPS, QUANTILES
		SAMPS     = self.get_SAMPS()
		QUANTILES = self.get_QUANTILES()
		lims      = {loop_par:[SAMPS.min().min(),SAMPS.max().max()]}
		self.lims   = lims
		self.bounds = bounds

		#Plot N68, N95
		N68   = abs(SAMPS.median()-true_par)/(SAMPS.quantile(0.84)-SAMPS.quantile(0.16)) < 1 ; N95   = abs(SAMPS.median()-true_par)/(SAMPS.quantile(0.975)-SAMPS.quantile(0.025)) < 1
		N68   = SAMPS.transpose()[N68.values].shape[0] 										 ; N95   = SAMPS.transpose()[N95.values].shape[0]
		print(round(100*(N68/SAMPS.shape[1]),1),'%',round(100*(N95/SAMPS.shape[1]),1),'%')
		if Ncred:
			line68 = r'$N_{68} = %s$'%round(100*(N68/SAMPS.shape[1]),1)+'%' 						 ;
			line95 = r'$N_{95} = %s$'%round(100*(N95/SAMPS.shape[1]),1)+'%'
		elif Parcred:
			line68 = r'$%s;\,68}=%s ^{+%s}_{-%s}$'%(parlabel.split('}')[0],SAMPS.quantile(0.68).median().round(2), round(SAMPS.quantile(0.68).quantile(0.84)-SAMPS.quantile(0.68).quantile(0.5),2),round(SAMPS.quantile(0.68).quantile(0.5)-SAMPS.quantile(0.68).quantile(0.16),2))
			line95 = r'$%s;\,95}=%s ^{+%s}_{-%s}$'%(parlabel.split('}')[0],SAMPS.quantile(0.95).median().round(2), round(SAMPS.quantile(0.95).quantile(0.84)-SAMPS.quantile(0.95).quantile(0.5),2),round(SAMPS.quantile(0.95).quantile(0.5)-SAMPS.quantile(0.95).quantile(0.16),2))
		else:
			line68,line95='',''
		print (line68+line95)
		ax[iax].annotate(line68,xy=(0.95-(0.95-0.0225)*Lside,0.425-0.08*Parcred),xycoords='axes fraction',fontsize=FS,ha=HA)
		ax[iax].annotate(line95,xy=(0.95-(0.95-0.0225)*Lside,0.375-Parcred*(0.13-0.025*(iax==0))),xycoords='axes fraction',fontsize=FS,ha=HA)

		#Real-Data Posterior Fit
		if real_data_samps is not False:
			if line_rec_title is None: line_rec_title   = 'Real-Data Posterior'
			samps = PARAMETER(real_data_samps,dfpar,parlabel,lims[loop_par],bounds[loop_par],-1,iax,{})
			samps.get_xgrid_KDE()
			ax[iax].plot(samps.xgrid,samps.KDE,alpha=1,linewidth=3,color='C3',label=f"{line_rec_title} \n"+r'$%s = %s ^{+%s}_{-%s}$'%(parlabel,real_data_samps.quantile(0.5).round(2),round(real_data_samps.quantile(0.84)-real_data_samps.quantile(0.5),2),round(real_data_samps.quantile(0.5)-real_data_samps.quantile(0.16),2)))
			simavheight = samps.KDE[np.argmin(np.abs(real_data_samps.quantile(0.5)-samps.xgrid))]
			ax[iax].plot(real_data_samps.quantile(0.5)*np.ones(2),[0,simavheight],c='C3',linewidth=2)

		#Plot per-Sim faint posteriors
		if plot_ind:
			KDEmax = 0
			for ISIM in SAMPS:
				samps = PARAMETER(SAMPS[ISIM],dfpar,parlabel,lims[loop_par],bounds[loop_par],FITS[ISIM]['fitsummary'].loc[dfpar]['r_hat'],iax,{})
				samps.get_xgrid_KDE()
				KDEmax = max(KDEmax,np.amax(samps.KDE))
				ax[iax].plot(samps.xgrid,samps.KDE,alpha=0.08,color=color)

		#Plot True Parameter Value, and Annotate
		if plot_true is True:
			ax[iax].plot(true_par*np.ones(2),[0,KDEmax],c='black',linewidth=2,linestyle='--')
		if annotate_true:
			ax[iax].annotate(r'True $%s=%s$'%(parlabel,true_par),xy=(0.95-(0.95-0.0225)*Lside,0.5+0.02),xycoords='axes fraction',fontsize=FS,ha='left' if ('tau' in loop_par and iax>0) else 'right')


		###Plot and Simulation Averaged Posterior
		samps = PARAMETER(SAMPS.stack(),dfpar,parlabel,lims[loop_par],bounds[loop_par],None,iax,{})
		sap_chain = samps.chain
		#Plot and label
		if FAC is None:	samps.Nsamps /= 10
		else:			samps.Nsamps /= FAC
		samps.get_xgrid_KDE()
		if line_sap_title is None: line_sap_title = 'Simulation-Averaged Posterior'
		if self.quantilemode:	line_sap_summary = r'$%s = %s ^{+%s}_{-%s}$'%(parlabel,sap_chain.quantile(0.5).round(2),round(sap_chain.quantile(0.84)-sap_chain.quantile(0.5),2),round(sap_chain.quantile(0.5)-sap_chain.quantile(0.16),2))
		else:					line_sap_summary = r'$%s = %s \pm %s$'%(parlabel,sap_chain.quantile(0.5).round(2),sap_chain.std().round(2))
		print (line_sap_title+line_sap_summary)
		ax[iax].plot(samps.xgrid,samps.KDE,alpha=1,color=color,linewidth=3,label='\n'.join([line_sap_title,line_sap_summary]),linestyle=linestyle)
		#Fill between with quantiles
		if fill_between:
			for qlo,qhi in zip([0.16,0.025],[0.84,0.975]):
				siglo = samps.get_KDE_values(value=sap_chain.quantile(qlo), return_only_index=True)
				sigup = samps.get_KDE_values(value=sap_chain.quantile(qhi), return_only_index=True)+1
				ax[iax].fill_between(samps.xgrid[siglo:sigup],np.zeros(sigup-siglo),samps.KDE[siglo:sigup],color=color,alpha=0.2)
			simavheight = samps.KDE[np.argmin(np.abs(sap_chain.quantile(0.5)-samps.xgrid))]
			ax[iax].plot(sap_chain.quantile(0.5)*np.ones(2),[0,simavheight],c=color,linewidth=2,linestyle=linestyle)

		#Plot and Annotate Medians
		if plot_medians:
			if self.quantilemode:	line_median  = r'Median-$%s=%s^{+%s}_{-%s}$'%(parlabel,QUANTILES[0.5].quantile(0.5).round(2),round(QUANTILES[0.5].quantile(0.84)-QUANTILES[0.5].quantile(0.5),2),round(QUANTILES[0.5].quantile(0.5)-QUANTILES[0.5].quantile(0.16),2))
			else:					line_median  = r'Median-$%s=%s\pm%s$'%(parlabel,QUANTILES[0.5].quantile(0.5).round(2),QUANTILES[0.5].std().round(2))
			line_pmedian = r"$p($"+'Median-'+r"$%s<%s ; \,\rm{True}})=%s$"%(parlabel, parlabel.split('}')[0],round(100*QUANTILES[QUANTILES[0.5]<true_par].shape[0]/QUANTILES[0.5].shape[0],1)) + '%'
			if not include_pmedian:	median_labels = line_median
			else:					median_labels = '\n'.join([line_median,line_pmedian])
			ax[iax].plot(QUANTILES[0.5].quantile(0.5)*np.ones(2),[0,KDEmax],c='C1'	,linewidth=2,label=median_labels,linestyle=':')
			ax[iax].fill_between([QUANTILES[0.5].quantile(0.16),QUANTILES[0.5].quantile(0.84)],[0,0],[KDEmax,KDEmax],color='C1',alpha=0.2)
			print (line_median+line_pmedian)
		#Set ticks and legend
		if dress_figure:
			ax[iax].set_ylim([0,max(KDEmax,simavheight)])
			ax[iax].set_yticks([])
			ax[iax].legend(fontsize=FS,framealpha=1,loc='upper right')
			ax[iax].tick_params(labelsize=FS)

		#Set appropriate limits
		FAC  = 3
		XLIM =  FAC*np.array([sap_chain.quantile(0.16),sap_chain.quantile(0.84)])-(FAC-1)*sap_chain.quantile(0.5)#XLIMS
		with suppress(TypeError): XLIM[0] = max([XLIM[0],bounds[loop_par][0]])
		with suppress(TypeError): XLIM[1] = min([XLIM[1],bounds[loop_par][1]])
		ax[iax].set_xlim(XLIM)



def get_KEEPERS(GLOB_FITS,Nsim_keep,Rhat_threshold,loop_par,dfpar):
	"""
	Get Keepers

	Function to get indices of simulations to keep,
	loops through all sims for all parameters values
	identifies individual simulations where Rhat>Rhat_threshold

	Parameters
	----------
	GLOB_FITS : dict
		{parvalue:FITS} where FITS={ISIM:FIT}
	Nsim_keep : int
		No. of simulations to plot
	Rhat_threshold : float
		maximum allowed Rhat
	loop_par : str
		parname of parameter being plotted
	dfpar : str
		string name of parameter in posterior samples file

	Returns
	----------
	KEEPERS : list
		list of ISIM indices to keep
	"""
	for iim,parvalue in enumerate(GLOB_FITS):
		FITS = GLOB_FITS[parvalue]
		if iim==0:	KEEPERS = [True for _ in range(len(FITS))]
		for ISIM in range(len(KEEPERS)):
			try:
				fitsummary = FITS[ISIM]['fitsummary']
				Rhat = fitsummary.loc[dfpar]['r_hat']
				if Rhat>Rhat_threshold:
					print (f'{loop_par}={parvalue}, ISIM={ISIM}, Rhat={Rhat}')
					KEEPERS[ISIM] = False
			except KeyError:
				KEEPERS[ISIM] = False

	KEEPERS = [ISIM for ISIM in range(len(KEEPERS)) if KEEPERS[ISIM]]
	KEEPERS = KEEPERS[:Nsim_keep]
	print ('Max index kept is:',max(KEEPERS),f'; Nsims kept=={len(KEEPERS)}')
	return KEEPERS

def trim_to_KEEPERS(GLOB_FITS,KEEPERS):
	"""
	Trim to Keepers

	Function to trim GLOB_FITS to retain only those with ISIMs in KEEPERS

	Parameters
	----------
	GLOB_FITS : dict
		{parvalue:FITS} where FITS={ISIM:FIT}
	KEEPERS : list
		list of ISIM indices to keep

	Returns
	----------
	GLOB_FITS where values==FITS are trimmed to have only ISIM in KEEPERS as keys
	"""
	for key in GLOB_FITS:
		GLOB_FITS[key] = {ISIM:GLOB_FITS[key][ISIM] for ISIM in KEEPERS}
	return GLOB_FITS

def update_edit_dict_for_ppc(sbc_choices, edit_dict):
	"""
	Update edit_dictionary for Posterior Predictive Checks

	Parameters
	----------
	sbc_choices : dict
		dicionary of choices from .yaml regarding how datasets are simulated
	edit_dict : dict
		edits to above .yaml choices, suited for ppc

	Returns
	----------
	edit_dict : dict
		edited to have all important data to perform ppc
	"""
	#Firstly Update Paths
	if 'load_parameters' in edit_dict:
		for key,value in edit_dict['load_parameters'].items():
			sbc_choices['load_parameters'][key] = value

	#Candidates are ppc values, additional is the stuff that actually gets used
	candidates={} ; additional = {}
	#Open original BirdSnack fit to real data
	with open(f"{sbc_choices['load_parameters']['path_to_birdsnack_rootpath']}products/stan_fits/FITS/FIT{edit_dict['simulate_parameters']['pre_defined_hyps']['load_file']}.pkl",'rb') as f:
		FIT = pickle.load(f)
	#Get number of SNe (including censored SNe), and simulate this number if not specified
	candidates['S'] = FIT['stan_data']['S']
	#Get AVsimdist, and use this distribution if not specified
	candidates['AVsimdist'] = FIT['choices']['analysis_parameters']['AVprior']
	#Get RVsimdist, new entry may be buggy, implemented for AVRVBeta simulations
	if 'RVstyles' in FIT['choices']['analysis_parameters'] and FIT['choices']['analysis_parameters']['RVstyles']==['AVRVBeta']:
			candidates['RVsimdist'] = 'AVRVBeta'
	else:	candidates['RVsimdist'] = 'Norm'
	#Assign folder name using these values of dust hyps ##Old Code #dust_hyps = dict(zip(['tauA','muRV','sigRV'],list(FIT['df'].median(axis=0).round(4)[['tauA','mu_RV','sig_RV']].values)))
	import sys
	sys.path.append(sbc_choices['load_parameters']['path_to_birdsnack_rootpath']+'model_files/')
	from plotting_functions import get_parlabels
	pars,parnames,parlabels,bounds = get_parlabels(FIT['choices'])
	dust_hyps = dict(zip(pars,list(FIT['df'].median(axis=0).round(4)[parnames].values)))
	#print (dust_hyps)
	candidates= {**candidates,**dust_hyps}
	#Update changes
	for key,value in candidates.items():
		if key not in edit_dict['simulate_parameters']:
			if key in sbc_choices['simulate_parameters']:
				if sbc_choices['simulate_parameters'][key]=='ppc':
					additional[key] = value
			else:
				additional[key] = value

	#Upload additional to edit_dict
	edit_dict['simulate_parameters'] = {**edit_dict['simulate_parameters'],**additional}
	#Get reference band (for building full vector mu_int where ref-band has element entry=0)
	zero_index      = FIT['choices']['analysis_parameters']['zero_index']
	#Whether using model for colours or deviations
	IntrinsicModel  = FIT['choices']['analysis_parameters']['IntrinsicModel']
	edit_dict['simulate_parameters']['pre_defined_hyps'] = {**edit_dict['simulate_parameters']['pre_defined_hyps'],**{'ppc_zero_index':zero_index,'ppc_IntrinsicModel':IntrinsicModel}}
	return edit_dict
