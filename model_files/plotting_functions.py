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
	get_parlabels(choices,default=False,return_parnames_only=False)
	get_Lines(choices)
	finish_corner_plot(fig,ax,Lines,save,quick,show,plotpath,savekey)
	get_mass_label(mass,choices,nocol='black',hicol='blue',locol='magenta')
	get_list_of_colours_for_corner(pblist,Style,Input=None,parstr='parnames')
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

from miscellaneous import ensure_folders_to_file_exist
import matplotlib.pyplot as pl
import matplotlib
import string
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
			name of figure to save

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
			ensure_folders_to_file_exist(self.plotpath+savename)
			pl.savefig(self.plotpath+savename,bbox_inches='tight')
		if self.show:
			pl.show()


def get_parlabels(choices,default=False,return_parnames_only=False):
	'''
	Get Parameter Labels

	Uses choices r.e. AV prior to get key parameters of interest and their plot labels

	Parameters
	----------
	choices : dict
		analysis choices

	default : bool (optional; default=False)
		if default, ignore any analysis choices and just return standard labels

	return_parnames_only : bool (optional; default=False)
		if True, return only parnames

	Returns
	----------
	pars,parnames,parlabels,bounds : lists
		respectively, parameter names for user, the parameter names in MCMC chains, the parameter labels used in plots, and the prior bounds on parameters
	'''
	def get_binned_RV_parlabels(choice):
		AVprior  = choices['analysis_parameters']['AVprior']
		muRVmin  = choices['analysis_parameters']['muRVmin']
		muRVmax  = choices['analysis_parameters']['muRVmax']
		pars = ['tauA'] ; parnames = ['tauA'] ; parlabels = ['$\\tau_A$'] ; bounds = [[0,None]]
		N_GaussRV_dists     = len([x for x in choices['analysis_parameters']['RVstyles'] if x=='Gauss'])
		N_AVRVBeta_dists    = len([x for x in choices['analysis_parameters']['RVstyles'] if x=='AVRVBeta'])
		N_AVRVSigmoid_dists = len([x for x in choices['analysis_parameters']['RVstyles'] if x=='AVRVSigmoid'])
		alphabet = list(string.ascii_lowercase)
		for n in range(1,N_GaussRV_dists+1):
			parnames.extend([f'mu_RV.{n}',f'sig_RV.{n}'])
			if N_GaussRV_dists==1:
				pars.extend(['muRV','sigRV'])
				parlabels.extend(['$\\mu_{R_V}$','$\\sigma_{R_V}$'])
			else:
				pars.extend([f'muRV.{n}',f'sigRV.{n}'])
				parlabels.extend(['$\\mu^{%s}_{R_V}$'%alphabet[n-1],'$\\sigma^{%s}_{R_V}$'%alphabet[n-1]])
			bounds.extend([[muRVmin,muRVmax],[0,None]])
		for n in range(1,N_AVRVBeta_dists+1):
			parnames.extend([f'beta0.{n}',f'beta1.{n}',f'sig_RV_beta.{n}'])
			if N_AVRVBeta_dists==1:
				pars.extend(['beta0','beta1','sigRVbeta'])
				parlabels.extend(['$\\beta_{0}$','$\\beta_{1}$','$\\sigma_{R_V}$'])
			else:
				pars.extend([f'beta0.{n}',f'beta1.{n}',f'sigRVbeta.{n}'])
				parlabels.extend(['$\\beta^{%s}_{0}$'%alphabet[n-1],'$\\beta^{%s}_{1}$'%alphabet[n-1],'$\\sigma^{%s}_{R_V}$'%alphabet[n-1]])
			bounds.extend([[1/(muRVmax-muRVmin),None],[0,None],[0,None]])
		for n in range(1,N_AVRVSigmoid_dists+1):
			parnames.extend([f'muAV_sigmoid.{n}',f'sigAV_sigmoid.{n}',f'A_sigmoid.{n}',f'B_sigmoid.{n}',f'sig_RV_sigmoid.{n}'])
			if N_AVRVSigmoid_dists==1:
				pars.extend(['muAVsig','sigAVsig','Asig','Bsig','sigRVsig'])
				parlabels.extend(['$\\mu_{A_V}$','$\\sigma_{A_V}$','$A_{\\rm{Sigmoid}}$','$B_{\\rm{Sigmoid}}$','$\\sigma_{R_V}$'])
			else:
				pars.extend([f'muAVsig.{n}',f'sigAVsig.{n}',f'Asig.{n}',f'Bsig.{n}',f'sigRVsig.{n}'])
				parlabels.extend(['$\\mu^{%s}_{A_V}$'%alphabet[n-1],'$\\sigma^{%s}_{A_V}$'%alphabet[n-1],'$\\A^{%s}_{\\rm{Sigmoid}}$'%alphabet[n-1],'$\\B^{%s}_{\\rm{Sigmoid}}$'%alphabet[n-1],'$\\sigma^{%s}_{R_V}$'%alphabet[n-1]])
			bounds.extend([[0,None],[0,None],[muRVmin,muRVmax],[0,None],[0,None]])#B is actually dep. on A, but ignore this effect

		if AVprior in ['Gamma']:
			pars.append('nuA')
			parnames.append('nu')
			parlabels.append('$\\nu_A$')
			bounds.append([0,None])
		return pars,parnames,parlabels,bounds

	#Get parnames,labels,bounds
	pars       = ['tauA','muRV','sigRV']
	parnames   = ['tauA','mu_RV','sig_RV']
	parlabels  = ['$\\tau_A$','$\\mu_{R_V}$','$\\sigma_{R_V}$']
	bounds     = [ [0,None]  , [1,5], 			[0,None]]

	proceed=True
	if 'BinnedRVFit' in choices['analysis_parameters'] and choices['analysis_parameters']['BinnedRVFit']:
		pars,parnames,parlabels,bounds = get_binned_RV_parlabels(choices)
		proceed=False

	if not default and choices!={} and proceed:
		AVprior  = choices['analysis_parameters']['AVprior']
		muRVmin  = choices['analysis_parameters']['muRVmin']
		muRVmax  = choices['analysis_parameters']['muRVmax']
		#Update bounds with muRVmin,muRVmax
		bounds     = [ [0,None]  , [muRVmin,muRVmax], [0,None]]
		#Old Files fail because options didn't yet exist
		try:		RVprior = choices['analysis_parameters']['RVprior']
		except:		RVprior = 'Norm'
		try:		skew_RV  = choices['analysis_parameters']['skew_RV']
		except:		skew_RV = False
		try:		skew_int = choices['analysis_parameters']['skew_int']
		except:		skew_int = False

		if AVprior in ['Gamma']:
			pars.append('nuA')
			parnames.append('nu')
			parlabels.append('$\\nu_A$')
			bounds.append([0,None])
		if RVprior=='StudentT':
			pars.append('nuR')
			parnames.append('nuR')
			parlabels.append('$\\nu_{R_V}$')
			bounds.append([0,None])
		if skew_RV:
			pars.append('askewRV')
			parnames.append('alpha_skew_RV')
			parlabels.append('$\\alpha^{\\rm{skew}}_{R_V}$')
			bounds.append([None,None])
		if skew_int:
			pars.append('askewint')
			parnames.append('alpha_skew_int')
			parlabels.append('$\\alpha^{\\rm{skew}}_{\\rm{int}}$')
			bounds.append([None,None])
	if return_parnames_only:
		return parnames
	else:
		return pars,parnames,parlabels,bounds


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

	def BVcut_map(BVbool,CensoredData,BVcutval):
		if CensoredData:
			BVbool = True
		return {True:r'$|B-V|<%s\,\rm{mag}$ Cut'%BVcutval, False:r'No $B-V$ Cut'}[BVbool]

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

	def RVstylemapper(RVstyles):
		alphabet = list(string.ascii_lowercase)
		Ngauss = len([x for x in RVstyles if x=='Gauss']) ; counter = 1
		strlist = []
		for RVstyle in RVstyles:
			if RVstyle=='Gauss':
				if Ngauss>1:
					x = '$\\mathcal{N}(\mu^{%s}_{R_V},\sigma^{%s}_{R_V})$'%(alphabet[counter-1],alphabet[counter-1])
					counter+=1
				else:
					x = '$\\mathcal{N}(\mu_{R_V},\sigma_{R_V})$'
			elif RVstyle=='Flat':
				x= r'$U(%s,%s)$'%(choices['flatRVsmin'],choices['flatRVsmax'])
			elif 'Fix' in RVstyle:
				x= r'$R_V^*=%s$'%RVstyle.split('_')[-1]
			elif RVstyle=='AVRVBeta':
				x='$\\mu^s_{R_V} = %s+(\\beta_0 + \\beta_1 A_V^s)^{-1}$'%choices['muRVmin']
			elif RVstyle=='AVRVSigmoid':
				x='$\\mu^s_{R_V} = A+B/(1 + exp(\\frac{A_V^s-\\mu_{A_V}}{\\sigma_{A_V}}))$'
			strlist.append(x)
		return strlist

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
			BVcut_map(choices['BVcut'],choices['CensoredData'],choices['BVcutval'])										,
			#get_SED_model_string(choices['IntrinsicModel'])						,
			get_LC_shape_string(choices['include_LCshape'])					,
			#get_lambda_string(choices['lam_choice'])					    ,
			get_CensoredSNe_string(choices['CensoredData'])					,
			get_AVs_prior_string(choices['AVprior'])						#,
			]

	if 'BinnedRVFit' in choices and choices['BinnedRVFit']:
		if len(choices['BVbinboundaries'])==0:
			Lines.append(r'$R_V$ distribution: %s'%', '.join(RVstylemapper(choices['RVstyles'])))
		else:
			if len(choices['BVbinboundaries'])==1:
				Lines.append(r'%s SNe in $B-V$ Bins using boundary at %s mag'%(list(choices['N_in_each_bin']), float(round(choices['BVbinboundaries'][0],1))))
				pass
			else:
				Lines.append(r'%s SNe in $B-V$ Bins using boundaries at %s mag'%(list(choices['N_in_each_bin']), [float(round(x,1)) for x in choices['BVbinboundaries']]))
			Lines.append(r'$R_V$ distributions: [%s]'%', '.join(RVstylemapper(choices['RVstyles'])))

	if choices['mass_mode']=='all_masses':
		Lines.remove(mass_mode_map[choices['mass_mode']])
	if not choices['magerrcut']: Lines.replace(f"Mag. Errors <{choices['magerrcut']} mag","No Mag. Error Cut")
	else: 					   Lines.remove(f"Mag. Errors <{choices['magerrcut']} mag")

	if not choices['include_LCshape']:
		Lines.remove(get_LC_shape_string(choices['include_LCshape']))
	if not posterior:
		Lines.remove(BVcut_map(choices['BVcut'],choices['CensoredData'],choices['BVcutval']))
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


def get_mass_label(mass,choices,nocol='black',hicol='blue',locol='magenta'):
	"""
	Get Mass Label

	Parameters
	----------
	mass : float
		logM/Modot host galaxy stellar mass

	nocol,hicol,locol : strs (optional; defaults are black,blue,magenta)
		colour strings for whether sn is low/high mass or unknown

	Returns
	----------
	str, bool: colour for data point, if Mass is high
	"""
	logMcut = choices['additional_cut_parameters']['logMcut']
	if mass is None:#Unknown mass
		return nocol, None
	else:
		if mass>=logMcut:#High Mass
			return hicol, True
		else:			  #Low Mass
			return locol, False



def get_list_of_colours_for_corner(pblist,Style,Input=None,parstr='parnames'):
	"""
	Get list of colours for corner

	Function to return the colours plotted in the colour-colour corner plot based on input Style

	Parameters
	----------
	pblist: list
		list of passbands

	Style: str
		dictates which set of colours will be plotted, e.g. 'FixedBlue','Adjacent' or 'Input' is B-?, BV->Vr etc., or ['BV','VH',etc.]

	Input: list or None (optional; default is None)
		if list, a list of manually inputted colours to inspect

	parstr: str (optional; default='parnames')
		arbitrary keyword

	Returns
	----------
	grid_sets: dict
		key is the plot name, value is list of colours in that plot (and also the plot axis values e.g. '$B-V$ (mag)')

	grid_names: list of str
		the list of plot names

	"""

	if Style == 'FixedBlue':#For example BV,Br,Bi,BJ,BH and then another set is Vr,Vi,VJ,VH, all through to JH
		dstr   = {f"{Style}_{pb}":[] for pb in pblist[:-2]} ; dname  = {f"{Style}_{pb}_{parstr}":[] for pb in pblist[:-2]}
		grid_sets = {**dstr,**dname}
		for ifblue,fblue in enumerate(pblist[:-2]):
			for ifred,fred in enumerate(pblist[ifblue+1:]):
				grid_sets[f"{Style}_{fblue}"].append(fblue+fred)
				grid_sets[f"{Style}_{fblue}_{parstr}"].append(r'$%s-%s$ (mag)'%(fblue,fred))

	if Style == 'Adjacent':
		grid_sets = {f"{Style}":[],f"{Style}_{parstr}":[]}
		for if1,pb1 in enumerate(pblist[:-1]):
			pb2 = pblist[if1+1]
			grid_sets[f"{Style}"].append(pb1+pb2)
			grid_sets[f"{Style}_{parstr}"].append(r'$%s-%s$ (mag)'%(pb1,pb2))

	if Style == 'Input':
		grid_sets = {f"{Style}":[],f"{Style}_{parstr}":[]}
		for col in Input:
			pb1,pb2 = col[0],col[1]
			grid_sets[f"{Style}"].append(pb1+pb2)
			grid_sets[f"{Style}_{parstr}"].append(r'$%s-%s$ (mag)'%(pb1,pb2))

	grid_names = [s for s in grid_sets if parstr not in s]
	return grid_sets, grid_names
