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
		plot_sbc_panel()

Functions are:
	get_KEEPERS(GLOB_FITS,Nsim_keep,Rhat_threshold,loop_par,dfpar)
	trim_to_KEEPERS(GLOB_FITS,KEEPERS)

--------------------

Written by Sam M. Ward: smw92@cam.ac.uk

"""
import pandas as pd
import numpy as np
from contextlib import suppress

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

	def plot_sbc_panel(self):
		"""
		Plot SBC Panel

		Method to plot a single panel of hyperparameter recovery

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

		#Get Parameter Bounds
		dfparnames,none,bounds  = get_parlabels({'analysis_parameters':dict(zip(['AVprior','muRVmin','muRVmax'],[analysis_choices['AVprior'],analysis_choices['muRVmin'],analysis_choices['muRVmax']]))})
		bounds = {loop_par:bounds[dfparnames.index(dfpar)]}

		#Get SAMPS, QUANTILES
		SAMPS     = self.get_SAMPS()#pd.DataFrame(data={ISIM: FITS[ISIM]['df'][dfpar].values      for ISIM in FITS                    })
		QUANTILES = self.get_QUANTILES()#pd.DataFrame(data={   q:[FITS[ISIM]['df'][dfpar].quantile(q) for ISIM in FITS] for q in Quantiles})
		lims      = {loop_par:[SAMPS.min().min(),SAMPS.max().max()]}

		#Plot N68, N95
		N68   = abs(SAMPS.median()-true_par)/(SAMPS.quantile(0.84)-SAMPS.quantile(0.16)) < 1 ; N95   = abs(SAMPS.median()-true_par)/(SAMPS.quantile(0.975)-SAMPS.quantile(0.025)) < 1
		N68   = SAMPS.transpose()[N68.values].shape[0] 										 ; N95   = SAMPS.transpose()[N95.values].shape[0]
		line68 = r'$N_{68} = %s$'%round(100*(N68/SAMPS.shape[1]),1)+'%' 						 ;
		line95 = r'$N_{95} = %s$'%round(100*(N95/SAMPS.shape[1]),1)+'%'
		ax[iax].annotate(line68,xy=(0.95,0.425),xycoords='axes fraction',fontsize=FS,ha='right')
		ax[iax].annotate(line95,xy=(0.95,0.375),xycoords='axes fraction',fontsize=FS,ha='right')

		#Plot per-Sim faint posteriors
		KDEmax = 0
		for ISIM in SAMPS:
			samps = PARAMETER(SAMPS[ISIM],dfpar,parlabel,lims[loop_par],bounds[loop_par],FITS[ISIM]['fitsummary'].loc[dfpar]['r_hat'],iax,{})
			samps.get_xgrid_KDE()
			KDEmax = max(KDEmax,np.amax(samps.KDE))
			ax[iax].plot(samps.xgrid,samps.KDE,alpha=0.08,color='C0')

		#Plot True Parameter Value, and Annotate
		ax[iax].plot(true_par*np.ones(2),[0,KDEmax],c='black',linewidth=2,linestyle='--')
		ax[iax].annotate(r'True $\mu_{R_V}=%s$'%true_par,xy=(0.95-(0.95-0.0225)*('tau' in loop_par and iax>0),0.5+0.02),xycoords='axes fraction',fontsize=FS,ha='left' if ('tau' in loop_par and iax>0) else 'right')

		#Plot and Annotate Medians
		if self.quantilemode:	line_median  = r'Median-$%s=%s^{+%s}_{-%s}$'%(parlabel,QUANTILES[0.5].quantile(0.5).round(2),round(QUANTILES[0.5].quantile(0.84)-QUANTILES[0.5].quantile(0.5),2),round(QUANTILES[0.5].quantile(0.5)-QUANTILES[0.5].quantile(0.16),2))
		else:					line_median  = r'Median-$%s=%s\pm%s$'%(parlabel,QUANTILES[0.5].quantile(0.5).round(2),QUANTILES[0.5].std().round(2))
		line_pmedian = r"$p($"+'Median-'+r"$%s<%s ; \,\rm{True}})=%s$"%(parlabel, parlabel.split('}')[0],round(100*QUANTILES[QUANTILES[0.5]<true_par].shape[0]/QUANTILES[0.5].shape[0],1)) + '%'

		ax[iax].plot(QUANTILES[0.5].quantile(0.5)*np.ones(2),[0,KDEmax],c='C1'	,linewidth=2,label='\n'.join([line_median,line_pmedian]),linestyle=':')
		ax[iax].fill_between([QUANTILES[0.5].quantile(0.16),QUANTILES[0.5].quantile(0.84)],[0,0],[KDEmax,KDEmax],color='C1',alpha=0.2)

		###Plot and Simulation Averaged Posterior
		samps = PARAMETER(SAMPS.stack(),dfpar,parlabel,lims[loop_par],bounds[loop_par],None,iax,{})
		sap_chain = samps.chain
		#Plot and label
		samps.Nsamps /= 10#Otherwise too slow
		samps.get_xgrid_KDE()
		line_sap_title   = 'Simulation-Averaged Posterior; '
		if self.quantilemode:	line_sap_summary = r'$%s = %s ^{+%s}_{-%s}$'%(parlabel,sap_chain.quantile(0.5).round(2),round(sap_chain.quantile(0.84)-sap_chain.quantile(0.5),2),round(sap_chain.quantile(0.5)-sap_chain.quantile(0.16),2))
		else:					line_sap_summary = r'$%s = %s \pm %s$'%(parlabel,sap_chain.quantile(0.5).round(2),sap_chain.std().round(2))
		ax[iax].plot(samps.xgrid,samps.KDE,alpha=1,color='C0',linewidth=3,label='\n'.join([line_sap_title,line_sap_summary]))
		#Fill between with quantiles
		for qlo,qhi in zip([0.16,0.025],[0.84,0.975]):
			siglo = samps.get_KDE_values(value=sap_chain.quantile(qlo), return_only_index=True)
			sigup = samps.get_KDE_values(value=sap_chain.quantile(qhi), return_only_index=True)+1
			ax[iax].fill_between(samps.xgrid[siglo:sigup],np.zeros(sigup-siglo),samps.KDE[siglo:sigup],color='C0',alpha=0.2)
		simavheight = samps.KDE[np.argmin(np.abs(sap_chain.quantile(0.5)-samps.xgrid))]
		ax[iax].plot(sap_chain.quantile(0.5)*np.ones(2),[0,simavheight],c='C0',linewidth=2)

		#Set ticks and legend
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
