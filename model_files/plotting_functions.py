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

Functions include:
	get_parlabels(choices)
	get_Lines(choices)
	finish_corner_plot(fig,ax,Lines,save,quick,show,plotpath,savekey)
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


def get_parlabels(choices):
	'''
	Get Parameter Labels

	Uses choices r.e. AV prior to get key parameters of interest and their plot labels

	Parameters
	----------
	choices : dict
		analysis choices

	Returns
	----------
	parnames,parlabels,bounds : lists
		respectively, the parameter names in MCMC chains, the parameter labels used in plots, and the prior bounds on parameters
	'''
	AVprior = choices['analysis_parameters']['AVprior']
	muRVmin = choices['analysis_parameters']['muRVmin']
	muRVmax = choices['analysis_parameters']['muRVmax']

	parnames   = ['tauA','mu_RV','sig_RV']
	parlabels  = ['$\\tau_A$','$\\mu_{R_V}$','$\\sigma_{R_V}$']
	bounds     = [ [0,None]  , [muRVmin,muRVmax], [0,None]]
	if AVprior in ['Gamma']:
		parnames.append('nu')
		parlabels.append('$\\nu_A$')
		bounds.append([0,None])

	return parnames,parlabels,bounds

def get_Lines(choices, NSNe, NCens, posterior=True):
	"""
	Get Lines

	Function to get each str analysis choice line for plots

	Parameters
	----------
	choices: dict
		contains config.yaml choices

	NSNe: float
		number of SNe (excluding any censored SNe)

	NCens: float
		number of Censored SNe

	posterior : bool (optional; default=True)
		if False, remove lines relating to Modelling

	Returns
	----------
	list of str
	"""
	############################################################
	def Kcorrmap(apply_K_corrections,mangle):
		if not apply_K_corrections:
			return 'No K-corrections'
		else:
			if not mangle:
				return 'With K-corrections; No Mangling'
			else:
				return 'With K-corrections & Mangling'

	def interp_method_map(method,bright_mode):
		if 'GP' in method:
			mfmap = {'mag':'Magnitude', 'flux':'Flux'}
			return f"{method} {mfmap[bright_mode]} Data Interp."
		elif method=='snpy':
			return "SNooPy Template Interp."

	def Tmax_map(Tmaxchoicestr,Tmax_method):
		if Tmaxchoicestr=='Tmax_snpy_fitted':
			return 'SNooPy'
		elif Tmaxchoicestr=='Tmax_GP_restframe':
			return Tmax_method

	mass_mode_map  = dict(zip(['all_masses','high_masses','low_masses'],['No Stellar-Mass Cut','High-Stellar-Mass Hosts','Low-Stellar-Mass Hosts']))

	def dust_map(dustlaw,mwmode):
		if mwmode=='Model':
			return 'With Milky Way Extinction Correction'
		elif mwmode=='Presubtract':
			dustmap = dict(zip(['fitzpatrick99','ccm89','ccm89_o94'],['Fitzpatrick99','CCM89','CCM89+O94']))
			return f"{dustmap[dustlaw]} Dust Law"
		elif mwmode is None or mwmode=='None':
			return 'No Milky Way Extinction Correction'

	def BVcut_map(BVbool,CensoredData):
		if CensoredData:
			BVbool = True
		return {True:r'$|B-V|<0.3\,\rm{mag}$ Cut', False:r'No $B-V$ Cut'}[BVbool]

	MissingDataMap = {True:'With', False:'Without'}
	def inflate_fac_mapper(inflate_fac):
		mode = inflate_fac.split('_')[0]
		mod  = inflate_fac.split('_')[-1]
		if mode=='set':
			return f"Errs.={mod}mag"
		elif mode=='add':
			return f"Errs.+={mod}mag"
		elif mode=='mult':
			return f"Errs.*={mod}"

	def colour_choice_mapper(DataTransformation):
		if DataTransformation=='mags':
			return 'Peak Magnitudes'
		else:
			if DataTransformation=='Adjacent':
				return DataTransformation+' Colours'
			elif DataTransformation in ['B-X','X-H']:
				return '$%s$'%DataTransformation + ' Colours'


	def get_SED_model_string(IntrinsicModel):
		if IntrinsicModel=='Deviations':
			return 'Model Intrinsic Deviations'
		else:
			if IntrinsicModel=='Adjacent':
				return f'Model {IntrinsicModel} Intrinsic Colours'
			elif IntrinsicModel in ['B-X','X-H']:
				return f'Model Intrinsic $%s$ Colours'%IntrinsicModel

	def get_LC_shape_string(include_LCshape):
		if include_LCshape:
			return r'Include $\Delta m_{15}(B)$ term'
		else:
			return r'No $\Delta m_{15}(B)$ term'

	def get_lambda_string(lam_choice):
		if lam_choice=='effective_thetaON':
			return r"Pre-computed $\lambda_{\rm{eff}}$"
		if lam_choice=='effective_thetaOFF':
			return r"Using $\lambda_{\rm{eff}}$'s computed w/o LC shape"
		if lam_choice=='effective_taustar':
			return r"Using $\lambda_{\rm{eff}}$'s under $\tau_A = \tau_A^*$"
		if lam_choice=='central':
			return r"Using passband central wavelengths"

	def get_AVs_prior_string(AVprior):
		str = f"$A_V^s \sim "
		if AVprior=='Exp':
			return str + "$Exp$(\\tau_A)$"
		elif AVprior=='Gamma':
			return str + "$Gamma$(\\nu_A,\\tau_A)$"
		elif AVprior=='HalfCauchy':
			return str + "$ Half-Cauchy(0,1)"
		elif AVprior == 'HalfCauchywMean':
			return str + "$ Cauchy$(\\mu_A,\\tau_A)$"
		elif AVprior == 'Weibull':
			return str + "$ Weibull$(\\nu_A,\\tau_A)$"

	def get_CensoredSNe_string(CensoredData):
		if CensoredData:
			return f"{int(NCens)} Censored SNe"
		else:
			return "No Censored SNe"

	choices = {key:choices[glob_key][key] for glob_key in choices for key in choices[glob_key] if glob_key!='rootpath'}

	Lines = [
			f"$N_{{SNe}}={NSNe}$"									,
			f"Interpolate:${choices['interpflts']}$"							,
			dust_map(choices['dustlaw'],choices['apply_EBVMW_corrections'])							,
			f"{Kcorrmap(choices['apply_K_corrections'],choices['mangle'])}"				,
			f"$T_{{max}}$ from {Tmax_map(choices['Tmaxchoicestr'],choices['Tmax_method'])}"	,
			f"{interp_method_map(choices['mags_method'],choices['bright_mode'])}"				,
			f"Extract {colour_choice_mapper(choices['DataTransformation'])}:${''.join(choices['pblist'])}$"					,
			f"Mag. Errors <{choices['magerrcut']} mag"     				    ,
			mass_mode_map[choices['mass_mode']]								,
			BVcut_map(choices['BVcut'],choices['CensoredData'])										,
			#get_SED_model_string(choices['IntrinsicModel'])						,
			get_LC_shape_string(choices['include_LCshape'])					,
			#get_lambda_string(choices['lam_choice'])					    ,
			get_CensoredSNe_string(choices['CensoredData'])					,
			get_AVs_prior_string(choices['AVprior'])						#,
			]


	if choices['mass_mode']=='all_masses':
		Lines.remove(mass_mode_map[choices['mass_mode']])
	if not choices['magerrcut']: Lines.replace(f"Mag. Errors <{choices['magerrcut']} mag","No Mag. Error Cut")
	else: 					   Lines.remove(f"Mag. Errors <{choices['magerrcut']} mag")

	if not choices['include_LCshape']:
		Lines.remove(get_LC_shape_string(choices['include_LCshape']))
	if not posterior:
		Lines.remove(BVcut_map(choices['BVcut'],choices['CensoredData']))
		Lines.remove(get_AVs_prior_string(choices['AVprior']))
	if not choices['CensoredData'] or NCens==-1:
		Lines.remove(get_CensoredSNe_string(choices['CensoredData']))
	return Lines

def finish_corner_plot(fig,ax,Lines,save,quick,show,plotpath,savekey):
	"""
	Finish Corner Plot

	Simple function to complete corner plot

	Parameters
	----------
	fig,ax : of pl.subplots()

	Lines : list of str
		each line is analysis choice

	save,quick,show : bools
		whether, to save plot, show 2D KDE or 2D scatter samples, show plot

	plotpath : str
		path/to/plot/save/location

	savekey : str
		plot filname

	End Product(s)
	----------
	plot that is saved and/or shown
	"""
	delta_y = 0.15
	pl.annotate(r"Bird-Snack",     xy=(0.5,len(ax)-1-0.5+delta_y/2),xycoords='axes fraction',fontsize=20,color='black',weight='bold',ha='center',fontname='Courier New')
	pl.annotate(Lines[0],          xy=(0.5,len(ax)-1-0.5-delta_y/2),xycoords='axes fraction',fontsize=17,color='black',weight='bold',ha='center')
	for counter,line in enumerate(Lines[1:]):
		pl.annotate(line, xy=(1,len(ax)-delta_y*(counter-1)),xycoords='axes fraction',fontsize=15,color='C0',weight='bold',ha='right')
	fig.subplots_adjust(top=0.9)
	fig.subplots_adjust(wspace=0.075, hspace=0.075)
	ensure_folders_to_file_exist(plotpath)
	ensure_folders_to_file_exist(plotpath+'Quick/')
	if save:
		if quick:
			pl.savefig(f"{plotpath}Quick/{savekey}.pdf",bbox_inches='tight')
		else:
			pl.savefig(f"{plotpath}{savekey}.pdf",		bbox_inches='tight')
	if show:
		pl.show()
