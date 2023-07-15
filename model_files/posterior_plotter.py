"""
Posterior Plotter Module

Module containing POSTERIOR_PLOTTER class
Used to take in posterior samples and plot up corner plot


Contains:
--------------------
Functions:
	kde()

POSTERIOR_PLOTTER class:
	inputs: samples, parnames, parlabels, bounds, Rhats, choices, smoothing=2

	Methods are:
		update_lims(Nsig=3)
		corner(fig_ax, Ngrid=100, colour="C0", warn_tolerance=0.05, FS=15)
		corner_plot()
PARAMETER class:
	inputs: chain,parname,parlabel,lim,bound,Rhat,row,choices,smoothing=2,XGRID=None

	Methods are:
		get_xgrid(fac=0.1)
		slice_reflection_KDE()
		get_xgrid_KDE()
		get_KDE_values(location=None,value=None, return_only_index=False)
		plot_marginal_posterior(ax,colour='C0',alph=0.2,FS=14)
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""
import copy,scipy
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd

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

	y_data: array (optional; default=None)
		y-samples

	y_target: array (optional; default=None)
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


class POSTERIOR_PLOTTER:

	def update_lims(self, Nsig = 3):
		"""
		Update Limits

		This method updates limits to extend Nsig sigma either way of sample average and within bounds
		Limits used for KDE and visualisation

		Parameters
		----------
		Nsig : float (optional; default = 5)
			limits set to Nsig sigma either way of sample average, but ensures it doesn't cross parameter bounds

		End Product(s)
		----------
		self.lims updated
		"""
		chains,lims,bounds = self.chains,self.lims,self.bounds
		for p in range(self.Npar):
			if lims[p][0] is None:
				if bounds[p][0] is None:
					lims[p][0] = np.average(chains[p])-Nsig*np.std(chains[p])
				else:
					lims[p][0] = max([np.average(chains[p])-Nsig*np.std(chains[p]),bounds[p][0]])
			if lims[p][1] is None:
				if bounds[p][1] is None:
					lims[p][1] = np.average(chains[p])+Nsig*np.std(chains[p])
				else:
					lims[p][1] = min([np.average(chains[p])+Nsig*np.std(chains[p]),bounds[p][1]])

		self.lims = lims

	def __init__(self, samples, parnames, parlabels, bounds, Rhats, choices, smoothing=2):
		"""
		Initialisation

		Parameters
		----------
		samples  : dict
			{parname : array of samples}

		parnames : list
			[parname1, parname2...]

		parlabels : list
			[script_parname1, script_parname2...]

		bounds : list of list
			[[low_bound_1, up_bound_1],[low_bound_2, up_bound_2],...]

		Rhats : dict
			{parname : Rhat value}

		choices : dict
			plotting_parameters choices

		smoothing : float
			smoothing of KDEs for 1D and 2D marginals
		"""
		self.samples   = samples
		self.parnames  = parnames
		self.parlabels = parlabels
		self.bounds    = bounds
		self.Rhats     = Rhats
		self.choices   = choices
		self.smoothing = smoothing

		self.chains    = [self.samples[par] for par in self.parnames]
		self.Npar      = len(self.parnames)
		self.pardict   = {key:value for key,value in zip(self.parnames,self.parlabels)}
		self.lims      = [[None,None] for _ in range(self.Npar)]
		#Set limits to extend Nsig sigma either side of sample average, but ensures doesn't cross parameter boundaries
		self.update_lims()


	def corner(self, fig_ax, Ngrid=100, colour="C0", warn_tolerance=0.05, FS=15):
		"""
		Corner Method

		Returns 2D contours

		Parameters
		----------
		fig_ax: list (optional; default=None)
			[fig,ax]

		Ngrid : int (optional; default=100)
			if using KDE for 2D contour, sets no. of grid points, 100 is intensive, 10 is fairly fast but still noticeably jagged

		colour: str (optional; default="C0")
			colour for contours

		warn_tolerance: float (optional; default=0.05)
			width in percentage (i.e. 0.05==5%) which contours need to be correct for before displaying warning message

		FS: float (optional; default=15)
			fontsize

		End Products(s)
		----------
		fig,ax
		"""
		chains    = self.chains
		names     = self.parlabels
		lims      = self.lims
		bounds    = self.bounds
		quick     = self.choices['quick']
		smoothing = self.smoothing

		if len(chains) != len(names):
			raise ValueError("First dimension of input list/array should equal number of parameter names.")
		d = len(names) ; n = len(chains[0])

		fig,ax = fig_ax[:]
		for a in ax[np.triu_indices(d, k=1)]:
			a.axis("off")

		for row in range(d):
			pyrange = np.linspace(lims[row][0] - (bounds[row][0] is not None)*(lims[row][1]-lims[row][0]), lims[row][1] + (bounds[row][1] is not None)*(lims[row][1]-lims[row][0]), Ngrid*int(1 + (bounds[row][0] is not None) + (bounds[row][1] is not None)))
			ax[row,row].set_yticklabels("")
			ax[row,row].set_yticks([])
			for col in range(row):
				ax[row,col].get_shared_y_axes().remove(ax[row,row])
				if not quick:
					#PLOT THE 2D CONTOURS
					print ('Corner 2D KDE on row,col indices',row,col)
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
					#PLOT THE SAMPLES INSTEAD
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
		fig.subplots_adjust(top=0.9)
		fig.subplots_adjust(wspace=0.075, hspace=0.075)
		self.fig = fig
		self.ax  = ax

	def corner_plot(self):
		"""
		Plot posterior samples

		Function to plot the posterior samples

		End Product(s)
		----------
		Summary_Strs: dict
			each value is the parameters' Summary_Str: the posterior summary statistic that is displayed in the plot
		"""

		bounds  = self.bounds
		chains  = self.chains
		pardict = self.pardict

		#Create figure
		sfigx,sfigy = 3.1*len(pardict),2.7*len(pardict)
		fig,ax = pl.subplots(len(pardict),len(pardict),figsize=(sfigx,sfigy),sharex='col',sharey=False)

		print ("###"*5)
		Summary_Strs = {}
		for row,parname in enumerate(self.pardict):#Get Conf Intervals Here
			parameter = PARAMETER(self.chains[row],parname,self.pardict[parname],self.lims[row],self.bounds[row],self.Rhats[parname],row,self.choices,self.smoothing)
			parameter.get_xgrid_KDE()
			#Plot 1D Marginal KDEs
			Summary_Strs[parname] = parameter.plot_marginal_posterior(ax)
			print ("###"*5)
		#Plot 2D marginals
		self.corner([fig,ax])
		print ("Corner samples/contours plotted successfully")
		return Summary_Strs

class PARAMETER:

	def get_xgrid(self,fac=0.1):
		"""
		Get xgrid

		Get grid of points for KDE

		Parameters
		----------
		fac : float (optional; default=2)
			factor determines number of grid points: Ngrid = fac*Nsamps

		Returns
		----------
		xgrid : array
			grid points for KDE
		"""
		if self.XGRID is None:
			lim   = self.lim
			bound = self.bound
			Ngrid = int(self.Nsamps*fac)
			#If Bounds is not None, make xgrid larger, to allow for reflection KDE (i.e. KDE extends beyond prior boundary, then is reflected and doubled over within boundary)
			xgrid = np.linspace(lim[0] - (bound[0] is not None) * (lim[1]-lim[0]),
								lim[1] + (bound[1] is not None) * (lim[1]-lim[0]),
								Ngrid*int(1 + (bound[0] is not None) + (bound[1] is not None))
								)
			self.xgrid = xgrid
		else:
			self.xgrid = np.linspace(self.XGRID[0],self.XGRID[1],self.XGRID[2])

	def slice_reflection_KDE(self):
		"""
		Slice Reflection KDE

		Reflection KDEs go outside parameter bounds; this method trims to reside within bounds

		End Product(s)
		----------
		self.xgrid,KDE trimmed to reside within parameters bounds
		"""
		xgrid, KDE, bound = self.xgrid, self.KDE, self.bound
		#Having got a full mirror reflected kde, slice down to half kde bounded by lower support
		if bound[0] != None:
			lobound = 0
			for i,xi in enumerate(xgrid):
				if xi>=bound[0]:
					lobound = i
					break
			xgrid, KDE = xgrid[lobound:], KDE[lobound:]
		#Having got a full mirror reflected kde, slice down to half kde bounded by upper support
		if bound[1] != None:
			upbound = len(xgrid)
			for i,xi in enumerate(xgrid):
				if xi>=bound[1]:
					upbound = i
					break
			xgrid, KDE = xgrid[:upbound], KDE[:upbound]

		self.xgrid = xgrid
		self.KDE   = KDE

	def __init__(self,chain,parname,parlabel,lim,bound,Rhat,row,choices,smoothing=2,XGRID=None):
		"""
		See POSTERIOR_PLOTTER class docstring for input descriptions
		"""
		self.chain     = chain
		self.parname   = parname
		self.parlabel  = parlabel
		self.lim       = lim
		self.bound     = bound
		self.Rhat      = Rhat
		self.row       = row
		self.choices   = choices
		self.smoothing = smoothing #Smoothing for 1D marginal KDE
		self.XGRID     = XGRID

		self.Nsamps      = len(self.chain)
		self.samp_median = np.median(self.chain)
		self.samp_std    = np.std(self.chain)

		self.dfchain = pd.DataFrame(data={'par':chain}).sort_values(ascending=True,by='par')


	def get_xgrid_KDE(self):
		"""
		Get xgrid KDE

		Gets trimmed grid of points for KDE that resides within parameter boundaries
		Gets reflection KDE taking into account prior boundaries

		End Product(s)
		----------
		self.xgrid, KDE
		"""
		#Get grid for KDE
		self.get_xgrid()
		#Get full reflection KDE
		self.KDE = kde(self.chain, self.xgrid, x_bounds=self.bound, smoothing=self.smoothing)
		#Trim xgrid,KDE to reside within boundaries
		self.slice_reflection_KDE()

	def get_KDE_values(self,location=None,value=None, return_only_index=False):
		"""
		Get KDE values

		Method to get array location/values given some parameter value

		Parameters
		----------
		location : str (optional; default=None)
			if 'mode', use mode of KDE

		value : float (optionl; default=None)
			if float, get closest grid location to this parameter value
			checks are performed to ensure this parameter value resides within parameter boundaries

		return_only_index : bool (optional; default=False)
			bool dicatates if only grid index is returned

		Returns
		----------
		grid_index, xgrid_loc, KDE_height
		if return_only_index: return grid_index only
		"""
		#Error check
		if location is not None:
			if location not in ['mode']:
				raise Exception("To get KDE feature at a location, set location to valid value such as location='mode'")
		if value is not None:
			if self.bound[0] is not None:
				assert(self.bound[0]<=value)
			if self.bound[1] is not None:
				assert(value<=self.bound[1])
		if location is not None and value is not None:
			raise Exception('Choose either to look at KDE location OR numerical value, not both')

		#Get grid location
		if location=='mode':
			grid_index = np.argmax(self.KDE)
		elif value is not None:
			grid_index = np.argmin(np.abs(self.xgrid-value))
		#Return x,y values
		KDE_height = self.KDE[grid_index]
		xgrid_loc  = self.xgrid[grid_index]
		if return_only_index:
			return grid_index
		else:
			return grid_index, xgrid_loc, KDE_height

	def plot_marginal_posterior(self,ax,colour='C0',alph=0.2,FS=14):
		"""
		Plot Marginal Posterior

		Function to plot the 1D marginal panels in the posterior

		Parameters
		----------
		ax : ax of figure
			ax[row,row] = marginal panels

		colour : str (optional; default='C0')
			colour of KDE

		alph : float (optional; default=0.2)
			alpha of plot items (e.g KDE fill_between)

		FS : float
			fontsize
		Returns
		----------
		Summary_Str : str
		   Posterior summary for table
		"""

		paperstyle = self.choices['paperstyle']
		chain      = self.chain
		row        = self.row
		parname    = self.parname
		parlabel   = self.parlabel
		KDE        = self.KDE
		xgrid      = self.xgrid
		#Initialisations
		y0 = len(ax)-0.85 ; delta_y = 0.15 ; Summary_Str = ''
		if not paperstyle: ax[row,row].set_title(r'$\hat{R}$('+self.parlabel+f') = {self.Rhat:.3}')

		#Plot KDE
		ax[row,row].plot(xgrid, KDE, color=colour)
		ax[-1,row].set_xlim(*self.lim)

		#KDE doesnt peak at prior boundary
		condition1 = np.argmax(KDE)!=0 and np.argmax(KDE)!=len(KDE)-1
		#KDE is not near flat topped at prior boundary
		hh = np.exp(-1/8)#Gaussian height at sigma/2
		imode, xmode, KDEmode = self.get_KDE_values(location='mode')
		condition2 = not (KDE[0]>=hh*KDEmode or KDE[-1]>=hh*KDEmode)

		#If typical Gaussian-like posterior, plot median, 16 and 84th percentiles
		print (f"5%, 50%, 68%, 95% quantiles: {round(self.dfchain.par.quantile(0.05),2)}, {round(self.dfchain.par.quantile(0.5),2)}, {round(self.dfchain.par.quantile(0.68),2)},{round(self.dfchain.par.quantile(0.95),2)}")
		if condition1 and condition2:
			print (f"{parname}: {round(self.samp_median,2)} +/- {round(self.samp_std,2)}; 16th and 84th intervals: -{round(self.dfchain.par.quantile(0.5)-self.dfchain.par.quantile(0.16),2)}+{round(self.dfchain.par.quantile(0.84)-self.dfchain.par.quantile(0.5),2)}")
			#Plot median and 16th 84th quantiles
			i_med, x_med, KDE_med = self.get_KDE_values(value=self.samp_median)
			ax[row,row].plot(np.ones(2)*x_med,[0,KDE_med],c=colour)
			i_16 = self.get_KDE_values(value=self.dfchain.par.quantile(0.16),return_only_index=True)
			i_84 = self.get_KDE_values(value=self.dfchain.par.quantile(0.84),return_only_index=True)+1
			ax[row,row].fill_between(xgrid[i_16:i_84],np.zeros(i_84-i_16),KDE[i_16:i_84],color=colour,alpha=alph)
			if not paperstyle:
				pl.annotate("%s ="%parlabel,                     xy=(0.3  ,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='right')#String broken up for alignment
				pl.annotate("{:.3f}".format(self.samp_median),   xy=(0.65 ,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='right')
				pl.annotate("$\pm$ {:.3f}".format(self.samp_std),xy=(0.665,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='left' )
			elif paperstyle:
				ax[row,row].set_title(parlabel + " = {:.2f} $\pm$ {:.2f}".format(self.samp_median, self.samp_std), fontsize=FS)
				summary = ["{:.2f}".format(x) for x in [self.samp_median, self.samp_std,self.dfchain.par.quantile(0.95),self.dfchain.par.quantile(0.05)]]
				Summary_Str = f"${summary[0]}\\pm {summary[1]}^{str('{')}\\,\\,{summary[2]}{str('}')}_{str('{')}\\,\\,{summary[3]}{str('}')}$"
		#Otherwise, posterior peaks at/near prior boundary, choose to summarise posterior using quantiles
		else:
			storeinfo = {}
			for ic,conf in enumerate([0.68,0.95]):
				if imode>0.5*(len(xgrid)-1):#If peaks at RHS
					CONF = copy.deepcopy(1-conf)
					lg = '>'
					irange = [None,len(xgrid)]
				else:
					CONF = copy.deepcopy(conf)
					irange = [0,None]
					lg = '<'

				storeinfo[conf] = self.dfchain.par.quantile(CONF)

				i_conf, x_conf, KDE_conf = self.get_KDE_values(value=self.dfchain.par.quantile(CONF))
				irange = [i if i is not None else i_conf for i in irange]

				ax[row,row].plot(np.ones(2)*x_conf,[0,KDE_conf],c=colour)
				ax[row,row].fill_between(xgrid[irange[0]:irange[1]],np.zeros(irange[1]-irange[0]),KDE[irange[0]:irange[1]],color=colour,alpha=alph*(1-0.5*ic)*0.5)#half factor because gets doubled over
				#For RHS Boundary
				if irange[-1]==len(xgrid): ax[row,row].annotate(str(int(round(CONF*100,0)))+str("%"),xy=(x_conf+(xgrid[-1]-xgrid[0])*0.1,KDE_conf+0.08*KDEmode),color=colour,fontsize=FS+1,weight='bold',ha='right')
				#For LHS Boundary
				elif irange[0]==0:         ax[row,row].annotate(str(int(round(CONF*100,0)))+str("%"),xy=(x_conf,KDE_conf),color=colour,fontsize=FS+1,weight='bold')
				if not paperstyle:
					pl.annotate("%s %s"%(parlabel,lg),xy=(0.3  ,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='right')
					if ic==0: pl.annotate("{:.3f}".format(x_conf),  xy=(0.65 ,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='right')
					if ic==1: pl.annotate("({:.3f})".format(x_conf),xy=(0.735,y0-delta_y*(row+1)),xycoords='axes fraction',fontsize=FS,color=colour,ha='left')

			if paperstyle:
				ax[row,row].set_title(parlabel + f" {lg} {storeinfo[0.68]:.2f} ({storeinfo[0.95]:.2f})", fontsize=FS)
				Summary_Str = f"${lg} {storeinfo[0.68]:.2f} ({storeinfo[0.95]:.2f})$"
				print (f"{parname} {lg} {storeinfo[0.68]:.2f} ({storeinfo[0.95]:.2f})")

		ax[row,row].set_ylim([0,None])
		return Summary_Str
