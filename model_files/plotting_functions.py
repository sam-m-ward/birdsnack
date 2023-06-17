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

def get_Lines(choices, NSNe, NCens):
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

	choices = {key:choices[glob_key][key] for glob_key in choices for key in choices[glob_key]}

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
			get_SED_model_string(choices['IntrinsicModel'])						,
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

#########################
def kde(x_data, x_target, y_data=None, y_target=None, x_bounds=None, y_bounds=None, smoothing=1.0):
	"""
	Kernel Density Estimate

	Returns 1D or 2D KDE

	Parameters
	----------
	x_data: array
		x-samples

	x_target: array
		grid of x to evaluate KDE at

	y_data: array
		y-samples

	y_target: array
		grid of y to evaluate KDE at

	x_bounds: None or [lo,up] (optional; default = None)
		define prior lower and upper bounds on the x-samples

	x_bounds: None or [lo,up] (optional; default = None)
		define prior lower and upper bounds on the y-samples

	smoothing: float (optional; default=1.0)
		level of KDE smoothing

	Returns
	----------
	KDE
	"""
	if y_data is None:
		n = len(x_data)
		d = 1
	else:
		if len(x_data) == len(y_data):
			n = len(x_data)
			d = 2
		else:
			raise ValueError("Data vectors should be same length.")
	b = smoothing*n**(-1./(d+4)) #Scott Rule x Smoothing Factor
	if d==1:
		h = np.std(x_data)*b
	else:
		h = np.cov([x_data, y_data])*b**2

	x = x_target[:,None] - x_data[None,:]
	if d==2:
		y = y_target[:,None] - y_data[None,:]
		KH = scipy.stats.multivariate_normal.pdf(np.stack([x,y], axis=-1), cov=h)
		if x_bounds is not None:
			if x_bounds[0] is not None:
				x_minus = 2*x_bounds[0] - x_data
				KH += scipy.stats.multivariate_normal.pdf(np.stack([x_target[:,None] - x_minus[None,:], y], axis=-1), cov=h)
				if y_bounds is not None:
					if y_bounds[0] is not None:
						y_minus = 2*y_bounds[0] - y_data
						KH += scipy.stats.multivariate_normal.pdf(np.stack([x_target[:,None] - x_minus[None,:], y_target[:,None] - y_minus[None,:]], axis=-1), cov=h)
					if y_bounds[1] is not None:
						y_plus = 2*y_bounds[1] - y_data
						KH += scipy.stats.multivariate_normal.pdf(np.stack([x_target[:,None] - x_minus[None,:], y_target[:,None] - y_plus[None,:]], axis=-1), cov=h)
			if x_bounds[1] is not None:
				x_plus = 2*x_bounds[1] - x_data
				KH += scipy.stats.multivariate_normal.pdf(np.stack([x_target[:,None] - x_plus[None,:], y], axis=-1), cov=h)
				if y_bounds is not None:
					if y_bounds[0] is not None:
						y_minus = 2*y_bounds[0] - y_data
						KH += scipy.stats.multivariate_normal.pdf(np.stack([x_target[:,None] - x_plus[None,:], y_target[:,None] - y_minus[None,:]], axis=-1), cov=h)
					if y_bounds[1] is not None:
						y_plus = 2*y_bounds[1] - y_data
						KH += scipy.stats.multivariate_normal.pdf(np.stack([x_target[:,None] - x_plus[None,:], y_target[:,None] - y_plus[None,:]], axis=-1), cov=h)
		if y_bounds is not None:
			if y_bounds[0] is not None:
				y_minus = 2*y_bounds[0] - y_data
				KH += scipy.stats.multivariate_normal.pdf(np.stack([x, y_target[:,None] - y_minus[None,:]], axis=-1), cov=h)
			if y_bounds[1] is not None:
				y_plus = 2*y_bounds[1] - y_data
				KH += scipy.stats.multivariate_normal.pdf(np.stack([x, y_target[:,None] - y_plus[None,:]], axis=-1), cov=h)
		f = np.sum(KH, axis=1)/n
	else:
		Kh = scipy.stats.norm.pdf(x, scale=h)
		if x_bounds is not None:
			if x_bounds[0] is not None:
				x_minus = 2*x_bounds[0] - x_data
				Kh += scipy.stats.norm.pdf(x_target[:,None] - x_minus[None,:], scale=h)
			if x_bounds[1] is not None:
				x_plus = 2*x_bounds[1] - x_data
				Kh += scipy.stats.norm.pdf(x_target[:,None] - x_plus[None,:], scale=h)
		f = np.sum(Kh, axis=1)/n

	return f

def corner(chains, names, colour="#1f77b4", quick=False,lims=None, bounds=None, title=None, fig_ax=None, show_summary=True, figsize=None, smoothing=2.0, warn_tolerance=0.05,FS=15):
	"""
	corner(chains, get_param_labels(), quick=quick,lims=lims, bounds=bounds,colour=colour,fig_ax=[fig,ax],show_summary=False,FS=FS+1)
	Corner

	Returns 2D contours

	Parameters
	----------
	chains: list
		list of samples

	names: list
		list of axis labels

	colour: str (optional; default="#1f77b4")
		colour for contours

	quick: bool (optional; default=True)
		if True, don't plot 2D contours

	lims: list of tuple (optional; default=None)
		each tuple defines panel limits

	bounds: list of tuple (optional; default=None)
		each tuple defines prior lower upper bounds on each parameter

	title: str (optional; default=None)
		title of plot

	fig_ax: list (optional; default=None)
		[fig,ax]

	show_summary: bool (optional; default = True)

	figsize: (optional; default=None)
		sets size of figure if creating new figure inside function

	smoothing: float (optional; default=2.0)
		level of KDE smoothing

	warn_tolerance: float (optional; default=0.05)
		width in percentage (i.e. 0.05==5%) which contours need to be correct for before displaying warning message

	FS: float (optional; default=15)
		fontsize

	Returns
	----------
	fig,ax
	"""
	if len(chains) != len(names):
		raise ValueError("First dimension of input list/array should equal number of parameter names.")
	d = len(names)
	n = len(chains[0])
	if figsize is None:
		figsize = (3*d, 3*d)
	if lims is None:
		lims = [[None,None] for _ in range(d)]
	if bounds is None:
		bounds = [[None,None] for _ in range(d)]

	for p in range(d):
		if lims[p][0] is None:
			if bounds[p][0] is None:
				lims[p][0] = np.average(chains[p])-3*np.std(chains[p])
			else:
				lims[p][0] = max([np.average(chains[p])-3*np.std(chains[p]),bounds[p][0]])
		if lims[p][1] is None:
			if bounds[p][1] is None:
				lims[p][1] = np.average(chains[p])+3*np.std(chains[p])
			else:
				lims[p][1] = min([np.average(chains[p])+3*np.std(chains[p]),bounds[p][1]])


	if fig_ax is None:
		fig, ax = pl.subplots(d, d, figsize=figsize, sharex="col", sharey=False)
	else:
		fig = fig_ax[0]
		ax = fig_ax[1]
	for a in ax[np.triu_indices(d, k=1)]:
		a.axis("off")
	for row in range(d):
		pyrange = np.linspace(lims[row][0] - (bounds[row][0] is not None)*(lims[row][1]-lims[row][0]), lims[row][1] + (bounds[row][1] is not None)*(lims[row][1]-lims[row][0]), 100*int(1 + (bounds[row][0] is not None) + (bounds[row][1] is not None)))
		ax[row,row].plot(pyrange, kde(chains[row], pyrange, x_bounds=bounds[row], smoothing=smoothing), color=colour)
		if row == d-1:
			if lims is not None:
				ax[row,row].set_xlim(*lims[row])
		if show_summary:
			ax[row,row].set_title(names[row] + " = {:.3f} $\pm$ {:.3f}".format(np.median(chains[row]), np.std(chains[row])), fontsize=11)
		ax[row,row].set_yticklabels("")
		ax[row,row].set_yticks([])
		for col in range(row):
			ax[row,col].get_shared_y_axes().remove(ax[row,row])
			if not quick:
				pxrange = np.linspace(lims[col][0] - (bounds[col][0] is not None)*(lims[col][1]-lims[col][0]), lims[col][1] + (bounds[col][1] is not None)*(lims[col][1]-lims[col][0]), 100*int(1 + (bounds[col][0] is not None) + (bounds[col][1] is not None)))
				pxgrid, pygrid = np.meshgrid(pxrange, pyrange)
				try:
					cons = ax[row,col].contour(pxgrid, pygrid, np.reshape(kde(chains[col], pxgrid.flatten(), chains[row], pygrid.flatten(), x_bounds=bounds[col], y_bounds=bounds[row], smoothing=smoothing), pxgrid.shape), levels=25, colors=colour, alpha=0.1)
				except Exception as e:
					print (e)
				fracs = []
				for c, con in enumerate(cons.collections):
					paths = con.get_paths()
					if len(paths) == 1:
						fracs.append(sum(paths[0].contains_points(np.vstack([chains[col], chains[row]]).T))/n)
					elif len(paths) == 0:
						fracs.append(np.inf)
					else:
						fracs.append(sum([sum(path.contains_points(np.vstack([chains[col], chains[row]]).T)) for path in paths])/n)
				c68 = np.fabs(np.array(fracs) - 0.68).argmin()
				c95 = np.fabs(np.array(fracs) - 0.95).argmin()
				if not 0.68 - warn_tolerance < fracs[c68] < 0.68 + warn_tolerance:
					print("WARNING: Fraction of samples contained in estimated ({}, {}) 68 percent credible interval is {:.3f}, plotted contour may be suspect!".format(names[col], names[row],fracs[c68]))
				if not 0.95 - warn_tolerance < fracs[c95] < 0.95 + warn_tolerance:
					print("WARNING: Fraction of samples contained in estimated ({}, {}) 95 percent credible interval is {:.3f}, plotted contour may be suspect!".format(names[col], names[row], fracs[c95]))
				try:
					ax[row,col].contour(pxgrid, pygrid, np.reshape(kde(chains[col], pxgrid.flatten(), chains[row], pygrid.flatten(), x_bounds=bounds[col], y_bounds=bounds[row], smoothing=smoothing), pxgrid.shape), levels=[cons.levels[c95], cons.levels[c68]], colors=colour)
				except Exception as e:
					print (e)
			elif quick:
				ax[row,col].scatter(chains[col],chains[row],alpha=0.05,color=colour,s=50)
			if col == 0:
				ax[row,col].set_ylabel(names[row],fontsize=FS)
				if lims is not None:
					ax[row,col].set_ylim(*lims[row])
			else:
				ax[row,col].set_yticklabels("")
				ax[row,col].set_yticks([])
				ax[row,col].set_ylim(ax[row,0].get_ylim())
			ax[row,col].tick_params(labelsize=FS)
			if row == d-1:
				ax[row,col].set_xlabel(names[col],fontsize=FS)
				if lims is not None:
					ax[row,col].set_xlim(*lims[col])
	ax[d-1,d-1].set_xlabel(names[d-1],fontsize=FS)
	ax[d-1,d-1].tick_params(labelsize=FS)
	if title is not None:
		fig.suptitle(title)
	fig.subplots_adjust(top=0.9)
	fig.subplots_adjust(wspace=0.075, hspace=0.075)

	return fig, ax


def find_conf_interval_samples_count(CHAIN,x,conf,meanormode,m_loc,m_loc_index,activated):
	"""
	Find Intervals using Sample Counts

	Function that takes posterior samples, orders them, and extracts the Xth quantiles

	Parameters
	----------
	CHAIN: list or array
		list of samples

	x: list or array
		x-grid for KDE

	conf: float
		e.g. 68% is conf=0.68

	meanormode: str
		either 'mean' or 'mode', to find location on KDE

	m_loc: float
		location on x-grid of meanormode

	m_loc_index:int
		index on x-grid of meanormode

	activated : bool
		if True, then count samples assuming posterior peaks at prior boundary

	Returns
	----------
	lb,ub,sigma_lower,sigma_upper,siglo,sigup
	the lower bound, upper bound, difference from m_loc of each, and index of each
	"""
	lbub,sigloup,sigmalowupp = [],[],[]
	chain = copy.deepcopy(CHAIN)
	chain.sort()
	#############################
	if meanormode=='mean':
		confs        = [(1-conf)/2,(1+conf)/2]
	if meanormode=='mode':#e.g. peaks at boundary on LHS
		confs        = [0,conf]
		if m_loc_index==len(x)-1:#e.g. peaks at boundary on RHS
			confs    = [1,1-conf]
	#############################
	conf_indices = [int((len(chain)-1)*c) for c in confs]

	if m_loc_index<0.5*(len(x)-1):
		left_biased=True ; right_biased = False
	else:
		left_biased=False ; right_biased = True

	for ii,ic in enumerate(conf_indices):
		curve = abs(x-chain[ic])
		index = (curve.tolist()).index(np.amin(np.asarray(curve)))
		if meanormode=='mode' and ii==0 and (m_loc_index==0 or (activated and left_biased)):
			index=0
		if meanormode=='mode' and ii==0 and (m_loc_index==len(x)-1 or (activated and right_biased)):
			index=len(x)-1
		lbub.append(x[index])
		sigloup.append(index)
		sigmalowupp.append((m_loc-x[index])*(-1)**(ii))

	return lbub[0],lbub[1],sigmalowupp[0],sigmalowupp[1],sigloup[0],sigloup[1]

def get_conf_interval(chain,conf,x,KDE,meanormode,activated=False):
	"""
	Find Intervals using Sample Counts

	Function that takes posterior samples, orders them, and extracts the Xth quantiles

	Parameters
	----------
	chain: list or array
		list of samples

	conf: float
		e.g. 68% is conf=0.68

	x: list or array
		x-grid for KDE

	KDE: list or array
		kernel density estimate

	meanormode: str
		either 'mean' or 'mode', to find location on KDE

	activated : bool (optional; default=False)
		if True, then count samples assuming posterior peaks at prior boundary

	Returns
	----------
	lb,ub,sigma_lower,sigma_upper,siglo,sigup,KDE[m_loc_index],m_loc,m_loc_index
	the lower bound, upper bound, difference from m_loc of each, and index of each, KDE at mloc, meanormode, index of meanormode
	"""
	######################################################################
	if meanormode=='mode':
		mm          = np.amax(KDE)
		m_loc_index = (KDE.tolist()).index(mm)
		m_loc       = x[m_loc_index]
	######################################################################
	if meanormode=='mean':
		m_loc       = np.median(chain)#np.average(chain)#sum(x*KDE)/sum(KDE)
		curve       = abs(x-m_loc)
		m_loc_index = (curve.tolist()).index(np.amin(np.asarray(curve)))
		m_loc       = x[m_loc_index]
	#####################################################################
	lb,ub,sigma_lower,sigma_upper,siglo,sigup = find_conf_interval_samples_count(chain,x,conf,meanormode,m_loc,m_loc_index,activated)
	return lb,ub,sigma_lower,sigma_upper,siglo,sigup,KDE[m_loc_index],m_loc,m_loc_index




def plot_marginal_posterior(xtest,KDE,chain,parname,parlabel,bound,Rhat,conf,row,colour,alph,FS,paperstyle,ax,plot=True):
	"""
	Plot Marginal Posterior

	Function to plot the 1D marginal panels in the posterior

	Parameters
	----------
	xtest: list or array
		x-grid for KDE

	KDE: list or array
		kernel density estimate

	chain: list or array
		list of samples

	parname: str
		samples name of parameter

	parlabel: str
		plot axis label of parameter

	bound: list of tuple
		prior lower and upper bound of parameter

	Rhat: float
		the Rhat value

	conf: float
		e.g. 68% is conf=0.68

	row: int
		index of parameter in plot

	colour: str
		colour of KDE

	alph: float
		alpha of plot items

	FS: float
		fontsize

	paperstyle: bool
		if True, don't plot Rhat and make plot look nicer

	ax: ax of figure

	plot: bool (optional; default=True)

	Returns
	----------
	xtest,KDE,lb,ub,sigma_lower,sigma_upper,siglo,sigup,Km,m_loc,m_loc_index,Summary_Str

	Summary_Str is the posterior summary statistic that is displayed in the plot
	"""
	y0       = len(ax)-0.85
	delta_y  = 0.15

	if bound[0] != None:#Having got a full mirror reflected kde, slice down to half kde bounded by support
		for i,xi in enumerate(xtest):
			if xi>=bound[0]:
				lobound = i
				break
		xtest, KDE = xtest[lobound:], KDE[lobound:]
	if bound[1] != None:#Having got a full mirror reflected kde, slice down to half kde bounded by support
		upbound=len(xtest)
		for i,xi in enumerate(xtest):
			if xi>=bound[1]:
				upbound = i
				break
		xtest, KDE = xtest[:upbound], KDE[:upbound]
	#####################
	mmean,sstd = np.median(chain),np.std(chain)
	if not paperstyle:
		ax[row,row].set_title(r'$\hat{R}$('+parlabel+f') = {Rhat:.3}')

	#KDE doesnt peak at prior boundary
	condition1 = (KDE.tolist()).index(np.amax(KDE))!=0 and (KDE.tolist()).index(np.amax(KDE))!=len(KDE)-1
	lb,ub,sigma_lower,sigma_upper,siglo,sigup,Km,m_loc,m_loc_index = get_conf_interval(chain,conf,xtest,KDE,'mode')
	#KDE is not near flat topped at prior boundary
	hh = np.exp(-1/8)#Gaussian height at sigma/2
	condition2 = not (KDE[0]>=hh*Km or KDE[-1]>=hh*Km)# and siglo!=0
	if condition1 and condition2:
		lb,ub,sigma_lower,sigma_upper,siglo,sigup,Km,m_loc,m_loc_index = get_conf_interval(chain,conf,xtest,KDE,'mean')
		print (parname,":",round(m_loc,3),"+",round(sigma_upper,3),"-",round(sigma_lower,3))

		ax[row,row].plot([m_loc,m_loc],[0,KDE[m_loc_index]],c=colour)
		ax[row,row].fill_between(xtest[siglo:sigup],np.zeros(sigup-siglo),KDE[siglo:sigup],color=colour,alpha=alph)

		if not paperstyle:
			pl.annotate("%s ="%parlabel,            xy=(0.3  ,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='right')#String broken up for alignment
			pl.annotate("{:.3f}".format(mmean),     xy=(0.65 ,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='right')
			pl.annotate("$\pm$ {:.3f}".format(sstd),xy=(0.665,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='left' )
		elif paperstyle:
			ax[row,row].set_title(parlabel + " = {:.2f} $\pm$ {:.2f}".format(mmean, sstd), fontsize=FS)
			if parname=='none':
				Summary_Str = "${:.2f} \pm {:.2f}$".format(m_loc, sstd)
			else:
				#Choose 90% here therefore each bound has 95% on the other side
				x1,x2,sigma_lower95,sigma_upper95,siglo95,sigup95,x7,x8,x9 = get_conf_interval(chain,0.9,xtest,KDE,'mean')
				pdchain = pd.DataFrame(data={'x':chain})
				X = "{:.2f},{:.2f},{:.2f},{:.2f}".format(pdchain.quantile(0.5).values[0],pdchain.std().values[0],pdchain.quantile(0.95).values[0],pdchain.quantile(0.05).values[0])
				Y = X.split(',')
				Summary_Str = f"${Y[0]}\\pm {Y[1]}^{str('{')}\\,\\,{Y[2]}{str('}')}_{str('{')}\\,\\,{Y[3]}{str('}')}$"
	else:
		storeinfo = [[],[]] ; Summary_Str = ''
		print ('peaks at index:',(KDE.tolist()).index(np.amax(KDE)),'which is par value:',xtest[(KDE.tolist()).index(np.amax(KDE))])
		for ic,CONF in enumerate([0.68,0.95]):
			lb,ub,sigma_lower,sigma_upper,siglo,sigup,Km,m_loc,m_loc_index = get_conf_interval(chain,CONF,xtest,KDE,'mode',activated=True)
			print (parname,":",round(m_loc,3),"+",round(sigma_upper,3),"-",round(sigma_lower,3),siglo)
			ax[row,row].plot([m_loc+sigma_upper,m_loc+sigma_upper],[0,KDE[sigup]],c=colour)
			if siglo==len(xtest)-1:
				siglist     = [siglo,sigup]
				newsiglist  = siglist.copy()
				siglo,sigup = newsiglist[::-1][:]
				storeinfo[0].append('>')
			elif siglo==0:
				storeinfo[0].append('<')

			pdchain = pd.DataFrame(data={'x':chain})
			storeinfo[1].append(pdchain.quantile(CONF).values[0])
			print (storeinfo)
			ax[row,row].fill_between(xtest[siglo:sigup],np.zeros(sigup-siglo),KDE[siglo:sigup],color=colour,alpha=alph*(1-0.5*ic)*0.5)#half factor because gets doubled over
			###########################For RHS Boundary
			if sigup==len(xtest)-1:
				special_fac = 100
				ax[row,row].annotate(str(int(CONF*100))+str("%"),xy=(m_loc+sigma_upper+(xtest[1]-xtest[0])*special_fac,KDE[siglo]+0.08*(KDE[-1]-KDE[0])), color=colour,fontsize=FS+1,weight='bold',ha='right')
				if not paperstyle:
					pl.annotate("%s >"%parlabel,                              xy=(0.3  ,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='right')
			###########################For LHS Boundary
			elif siglo==0:
				ax[row,row].annotate(str(int(CONF*100))+str("%"),xy=(m_loc+sigma_upper,KDE[sigup]), color=colour,fontsize=FS+1,weight='bold')#String broken up for alignment
				if not paperstyle:
					pl.annotate("%s <"%parlabel,                              xy=(0.3  ,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='right')
			###########################
			if not paperstyle:
				if ic==0:
					pl.annotate("{:.3f}".format(m_loc+sigma_upper),  xy=(0.65 ,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='right')
				if ic==1:
					pl.annotate("({:.3f})".format(m_loc+sigma_upper),xy=(0.735,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='left')

		if paperstyle:
			ax[row,row].set_title(parlabel + f" {storeinfo[0][0]} {storeinfo[1][0]:.2f} ({storeinfo[1][1]:.2f})", fontsize=FS)
			Summary_Str = f"${storeinfo[0][0]} {storeinfo[1][0]:.2f} ({storeinfo[1][1]:.2f})$"
	if not paperstyle: Summary_Str=''
	return xtest,KDE,lb,ub,sigma_lower,sigma_upper,siglo,sigup,Km,m_loc,m_loc_index,Summary_Str


def PLOTTER(parnames,parlabels,bounds,samples,Lines,RVdist,savekey,Rhats=None,paperstyle=True,quick=False,show=True,save=False,pathappender=''):
	"""
	PLOTTER

	Function to plot the posterior samples

	Parameters
	----------
	parnames:list
		each item is sample parameter name

	parlabels:list
		each item is parameter label for plotting

	bounds:list of lists, each of a tuple
		each list is a tuple with the parameters' lower and upper prior bound

	samples:dict
		each value is posterior samples for that parameter

	Lines: list
		each item is an analysis choice that will be displayed in plot line by line

	RVdist: str
		either 'globRV' or 'popdist'

	savekey: str
		string used for figure savename

	Rhats: dict (optional; default is None)
		each value is the Rhat from sampling

	paperstyle: bool (optional; default=True)
		if True, don't plot Rhat and make plot look nicer

	quick: bool (optional; default=False)
		if True, don't plot the 2D contours, so as to speed up plotting process

	show: bool (optional; default=True)
		if True, show the plot via pl.show()

	save: bool (optional; default=False)
		if True, save the figure

	pathappender: str
		append str to start of savepath, e.g. '' or '../'

	Returns
	----------
	Summary_Strs: dict
		each value is the parameters' Summary_Str: the posterior summary statistic that is displayed in the plot
	"""
	pardict          = {key:value for key,value in zip(parnames,parlabels)}
	lims             = [[None,None] for _ in range(len(pardict))]

	sfigx,sfigy = 3.1*len(pardict),2.7*len(pardict)
	fig,ax = pl.subplots(len(pardict),len(pardict),figsize=(sfigx,sfigy),sharex='col',sharey=False)

	def get_chains():          return [samples[k] for k in pardict]
	def get_param_labels():    return [pardict[k] for k in pardict]

	chains = get_chains()
	for p in range(len(chains)):
		if lims[p][0] is None:
			if bounds[p][0] is None:
				lims[p][0] = np.average(chains[p])-3*np.std(chains[p])
			else:
				lims[p][0] = max([np.average(chains[p])-3*np.std(chains[p]),bounds[p][0]])
		if lims[p][1] is None:
			if bounds[p][1] is None:
				lims[p][1] = np.average(chains[p])+3*np.std(chains[p])
			else:
				lims[p][1] = min([np.average(chains[p])+3*np.std(chains[p]),bounds[p][1]])
	################################
	colour    = 'C0'
	conf      = 0.68
	alph      = 0.2
	smoothing = 2
	FS        = 14

	print ("####################")
	Summary_Strs = {}
	for row,chain in enumerate(chains):#Get Conf Intervals Here
		parname  =  parnames[row]
		parlabel = parlabels[row]
		###################################Apply this to extend KDE so we can check it integrates to 1, temporary change
		lims[row][0] += -5*np.std(chain)
		lims[row][1] +=  5*np.std(chain)
		###################################
		Nsteps  = len(chain)*2
		pyrange = np.linspace(lims[row][0] - (bounds[row][0] is not None)*(lims[row][1]-lims[row][0]), lims[row][1] + (bounds[row][1] is not None)*(lims[row][1]-lims[row][0]), Nsteps*int(1 + (bounds[row][0] is not None) + (bounds[row][1] is not None)))
		################################################
		xtest,KDE = pyrange, kde(chain, pyrange, x_bounds=bounds[row], smoothing=smoothing)
		xtest,KDE,lb,ub,sigma_lower,sigma_upper,siglo,sigup,Km,m_loc,m_loc_index, Summary_Str = plot_marginal_posterior(xtest,KDE,chain,parname,parlabel,bounds[row],Rhats[parname],conf,row,colour,alph,FS,paperstyle,ax,plot=bool(not quick))
		ax[row,row].set_ylim([0,None])
		######################################Revert back
		lims[row][0] +=  5*np.std(chain)
		lims[row][1] += -5*np.std(chain)
		######################################
		Summary_Strs[parname] = Summary_Str
		print ("####################")
	####################################
	fig,ax = corner(chains, get_param_labels(), quick=quick,lims=lims, bounds=bounds,colour=colour,fig_ax=[fig,ax],show_summary=False,FS=FS+1)
	print ("Corner samples/contours plotted successfully")
	for col in range(len(pardict)):
		increments = ax[len(pardict)-1,col].get_xticks()
		step_order = -int(np.floor(np.log10(np.average(increments[1:]-increments[:-1]))))
		ax[len(pardict)-1,col].set_xticklabels(np.round(ax[len(pardict)-1,col].get_xticks(),step_order))
		#if col==1:
		#	ax[len(pardict)-1,col].set_xticks([1.5,2.0,2.5,3.0])
		#	ax[len(pardict)-1,col].set_xticklabels([1.5,2.0,2.5,3.0])
	######################################
	if RVdist=='popdist':
		delta_y = 0.15
	elif RVdist=='globRV':
		delta_y= 0.11
	print (Lines[0])
	if RVdist=='popdist':
		pl.annotate(r"Bird-Snack",      xy=(0.5,len(ax)-1-0.5+delta_y/2),xycoords='axes fraction',fontsize=FS+6,color='black',weight='bold',ha='center',fontname='Courier New')#'monospace')
		pl.annotate(Lines[0],          xy=(0.5,len(ax)-1-0.5-delta_y/2),xycoords='axes fraction',fontsize=FS+3,color='black',weight='bold',ha='center')
	for counter,line in enumerate(Lines[1:]):
		pl.annotate(line, xy=(1,len(ax)-delta_y*(counter-1+1*(RVdist=='globRV'))),xycoords='axes fraction',fontsize=FS+1,color=colour,weight='bold',ha='right')#bbox=dict(facecolor='none', edgecolor='black'))
	########################################################################################################
	fig.subplots_adjust(top=0.9)
	fig.subplots_adjust(wspace=0.075, hspace=0.075)
	if save:
		if quick:
			pl.savefig(pathappender+f"plots/RVGP_Posteriors/Quick/{savekey}.pdf",bbox_inches='tight')
		else:
			pl.savefig(pathappender+f"plots/RVGP_Posteriors/{savekey}.pdf",bbox_inches='tight')
	if show:
		pl.show()
	return Summary_Strs
#########################
