"""
LC Object Module

Module containing LCObj class and additional functions
Used to take snana lc, and perform various transformations including getting Tmax with GP, and interpolating with GP


Contains:
--------------------
LCObj class:
	inputs: lc, choices

	Methods are:
		get_Tmax(return_samps=False)
		test_data_availability(ref_band='B', tref=0, Nlist=[2,10,2,40], tcol='phase', local=False)
		get_phase(returner=False)
		get_GP_interpolation()
		extract_interpolations()
		update_DF_M(DF_M, sn)
		plot_lc(plotter, mjd_or_phase='phase', bright_mode=None)
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

from snana_functions import *
from GP_functions import *
from plotting_functions import *
import matplotlib.pyplot as pl


class LCObj:

	def __init__(self,lc,choices):
		self.lc      = lc
		self.choices = choices

	def get_Tmax(self,return_samps=False):
		"""
		Get Tmax Method

		Method uses GPs (1D or 2D) to interpolate data and estimate Tmax
		Uses model maximum brightness as initial estimate, then cuts data to reside within phase range,
		And finally, draws Ngpdraws samples, and records sample mean and std of these.

		Parameters
		----------
		return_samps : bool (optional; default=False)
			if True, return the Tmax samples

		End Product(s)
		----------
		lc.meta has Tmax and Tmax_std

		Returns
		----------
		if return_samps: tmax_grid,f_samps,tmaxs
			respectively, the GP Tgrid, the flux interpolations, and the tmax estimates over Ngpdraws
		"""
		def trim_lcarrays_on_phase(bright_mode,GPfit,zhelio,phasemin,phasemax,data):
			"""
			Trim LC on Phase

			Simple function to trim arrays to reside within phase range

			Parameters
			----------
			bright_mode: str
				mag or flux data

			GPfit : GPfitobj
				for GP interpolations

			zhelio : flt
				heliocentric redshift for computing phase

			phasemin,phasemax : floats
				minimum and maximum phase
			data : list of arrays
				[mjd,bright,brighterr]

			Returns
			----------
			phase trimmed mjd, bright, brighterr
			"""
			mjd,bright,brighterr = data[:]
			#Get intial Tmax guess
			if bright_mode=='mag':      initial_Tmax_guess = GPfit.x[np.argmin(GPfit.y)]
			elif bright_mode=='flux':   initial_Tmax_guess = GPfit.x[np.argmax(GPfit.y)]
			phase        = (mjd-initial_Tmax_guess)/(1+zhelio)
			#Trim on Phase
			drop_indices = [i for i,pp in enumerate(phase) if pp<phasemin or pp>phasemax]
			mjd = np.delete(mjd,drop_indices) ; bright = np.delete(bright,drop_indices) ; brighterr = np.delete(brighterr,drop_indices)
			return mjd, bright, brighterr

		Tmaxchoicestr = self.choices['Tmaxchoicestr']
		if self.lc.meta[Tmaxchoicestr] is None or self.lc.meta[f'{Tmaxchoicestr}_std'] is None:
			#Get bright_mode==either interpolate flux or magnitude data
			bright_mode = self.choices['bright_mode']

			#If 1DGP, use reference band to interpolate, use initial guess of Tmax, then trim on phase range
			if self.choices['Tmax_method']=='1DGP':
				for flt in [self.choices['ref_band']]:
					#Get Ref-band Light Curve
					lcf = get_lcf(self.lc,flt)
					#Get Time in MJD, LC and LCerr
					mjd,bright,brighterr = get_time_lc_arrays(lcf,mjd_or_phase='mjd',flux_or_mag=bright_mode)
					#If LC has more than 1 data point
					if len(bright)>1:
						#Fit Ref-band Light Curve with 1D Squared Exponential GP
						GPfit     = GP_1D_squaredexp(mjd,bright,brighterr,tau_guess=10)
						if GPfit is not None:
							mjd, bright, brighterr = trim_lcarrays_on_phase(bright_mode,GPfit,self.lc.meta['REDSHIFT_HELIO'],self.choices['phasemin'],self.choices['phasemax'],[mjd,bright,brighterr])
							GPfit     = GP_1D_squaredexp(mjd,bright,brighterr,tau_guess=10)
					else:
						print (f"{self.lc.meta['SNID']} length of {bright_mode} vector in band {flt} is {len(bright)}")
						GPfit = None
			#If 2DGP, interpolate all bands, then use initial guess of Tmax in reference band to trim on phase range
			elif self.choices['Tmax_method']=='2DGP':
				lc  = self.lc
				#Get Time in MJD, LC and LCerr
				mjd,bright,brighterr = get_time_lc_arrays(lc,mjd_or_phase='mjd',flux_or_mag=bright_mode)
				#Get lam to flt mapping
				lam_to_flt  = dict(zip(lc.meta['lams'],lc.meta['flts']))
				flt_to_lam  = dict(zip(lc.meta['flts'],lc.meta['lams']))
				wavelengths = np.array([flt_to_lam[flt] for flt in lc['flt']])
				lambdaC     = np.asarray(list(lam_to_flt.keys()))
				#Get 2DGP Interpolation
				gFIT   = GP_2D_Matern(mjd,bright,brighterr,lambdaC,wavelengths)#,x_pred = x_pred)
				#Extract the reference band component for Tmax estimation
				lamref = flt_to_lam[self.choices['ref_band']]
				GPfit  = gFIT[lamref]
				if GPfit is not None:#Re-fit data within phase range
					mjd, bright, brighterr = trim_lcarrays_on_phase(bright_mode,GPfit,lc.meta['REDSHIFT_HELIO'],self.choices['phasemin'],self.choices['phasemax'],[mjd,bright,brighterr])
					gFIT  = GP_2D_Matern(mjd,bright,brighterr,lambdaC,wavelengths)
					GPfit = gFIT[lamref]

			#With trimmed data, re-sample GP Ngpdraws times
			if GPfit is not None:
				tmax_resolution, tmax_window, Ngpdraws = self.choices['tmax_resolution'], self.choices['tmax_window'], self.choices['Ngpdraws']
				if self.choices['Tmax_method']=='2DGP':
					self.choices['tmax_resolution'] *= 10
				#Using tmax window and resolution, identify maximum brightness by drawing GP samples
				if bright_mode=='flux': tmax_grid = np.linspace(GPfit.x[np.argmax(GPfit.y)]-(tmax_window/2),GPfit.x[np.argmax(GPfit.y)]+(tmax_window/2),int(tmax_window/tmax_resolution))
				elif bright_mode=='mag':tmax_grid = np.linspace(GPfit.x[np.argmin(GPfit.y)]-(tmax_window/2),GPfit.x[np.argmin(GPfit.y)]+(tmax_window/2),int(tmax_window/tmax_resolution))
				if self.choices['Tmax_method']=='1DGP':
					f_samps     = GPfit.gp.sample_conditional(bright, tmax_grid, size=Ngpdraws)
				elif self.choices['Tmax_method']=='2DGP':
					print (f"{lc.meta['SNID']}; 2DGP samples to get Tmax take longer, therefore: tmax_resolution *= 10; Ndraws /= 10")
					new_tpred   = np.vstack([np.hstack((tmax_grid for _ in range(len(lambdaC)))),np.array([lam for lam in lambdaC for _ in range(len(tmax_grid))])]).T
					f_samps_all = GPfit.gp.sample_conditional(bright, new_tpred, size=int(Ngpdraws/10))
					il          = list(lambdaC).index(lamref) ; Ngridt = copy.deepcopy(len(tmax_grid))
					f_samps     = f_samps_all[:,il*Ngridt:(il+1)*Ngridt]
				#Samples of tmax
				tmaxs     = np.array([tmax_grid[i] for i in np.argmax(f_samps-2*f_samps*(bright_mode=='mag'),axis=1)])
				#Tmax estimate and uncertainty is sample average and std. of tmax samples
				Tmax = np.average(tmaxs) ; Tmax_std = np.std(tmaxs)
			else:
				Tmax, Tmax_std = None, None
				tmax_grid,f_samps,tmaxs = None, None, None

			print (Tmax, Tmax_std)
			self.lc.meta[Tmaxchoicestr] = Tmax
			self.lc.meta[f'{Tmaxchoicestr}_std'] = Tmax_std
		else:
			print ('Tmax already estimated')
			pass

		if return_samps: return tmax_grid,f_samps,tmaxs


	def test_data_availability(self, ref_band = 'B', tref = 0, Nlist = [2,10,2,40], tcol='phase', local=False):
		"""
		Test data availability method

		Given some reference band and phase, check if data is available either side

		Parameters
		----------
		ref_band : str (optional; default='B')
			reference band to check data availability in

		tref : float (optional; default=0)
			reference time

		Nlist : list of floats, len=4 (optional; default=[2,10,2,40])
			Dictates required number of data points before and after reference time, and phase gaps either side from reference time

		tcol : str (optional; default='phase')
			time column to check data availability in

		local : bool (optional; default=False)
			if True, phaselim is with respect to reference phase, otherwise, use absolute phase w.r.t t=0

		Returns
		----------
		trim : bool
			if True, data is not available as required, so trim SN
		"""
		if ref_band is None or tref is None or Nlist is None:
			raise Exception('Must specify all of ref_band, tref, Nlist to run .test_data_availability() method')

		#Unpack Nlist
		Nbeforemax,Nbeforegap,Naftermax,Naftergap = Nlist[:]

		#Trim lc to phase range, get reference band, check data availability
		trim = get_data_availability_bool(
					get_lcf(self.get_phase(returner=True),ref_band),
						Nbeforemax,Nbeforegap,Naftermax,Naftergap,tref=tref,tcol=tcol,local=local
						)
		return trim


	def get_phase(self, returner=False):
		"""
		Get Phase Method

		Computes phase column using Tmax, cuts lc on phaselims, then updates metadata flts

		Parameters
		----------
		returner : bool (optional; default=False)
			if True, return lc

		End Products(s)
		----------
		lc: :py:class:`astropy.table.Table`
			light curve object with phase column cut on phaselims (and returned if returner=True)
		"""
		#Create phase column and cut data on phasemin and phasemax
		self.lc            = create_phase_column(self.lc,self.lc.meta[self.choices['Tmaxchoicestr']],phasemin=self.choices['phasemin'],phasemax=self.choices['phasemax'])
		#Given data has been cut, there is potential that some bands are removed completely
		self.lc            = set_lcmeta_ordered_flts(self.lc)
		return self.lc

	def get_GP_interpolation(self):
		"""
		Get GP Interpolation

		Function to fit GP to light curve data

		End Product(s)
		----------
		FIT: dict
			FIT[flt] = GPfit object; FIT['lc']=lc; FIT['flts']=flts, contains both data and GP interpolation
		"""
		#Interpolate in phase by default, but if Tmax is not estimated, interpolate mjd
		mjd_or_phase = 'phase' if self.lc.meta[self.choices['Tmaxchoicestr']] is not None else 'mjd'

		def get_tpred_gp(mjd_or_phase,phase):
			"""
			Get Predicted Times for GP

			Simple function to get grid of times at which GP interpolation will be estimated

			Parameters
			----------
			mjd_or_phase : str
			 	either 'phase' or 'mjd'; choose time column

			phase : array
				array of times

			Returns
			----------
			tpred : array
				gridded array of times for GP interpolation
			"""
			if mjd_or_phase=='phase':
				tpred = np.arange(self.choices['phasemin'],self.choices['phasemax']+self.choices['dp'],self.choices['dp'])
			elif mjd_or_phase=='mjd':
				tpred = np.arange(np.amin(phase),np.amax(phase)+self.choices['dp'],self.choices['dp'])
			return tpred

		FIT = {'lc':self.lc} ; flts = self.lc.meta['flts']
		print (f"{self.lc.meta['SNID']}: Beginning GP Interpolation; flts are:{flts}")
		if self.choices['mags_method']=='1DGP':#For 1DGP Interpolation
			for iif,flt in enumerate(flts):#For each filter
				phase,bright,brighterr = get_time_lc_arrays(get_lcf(self.lc,flt),mjd_or_phase=mjd_or_phase,flux_or_mag=self.choices['bright_mode'])#Get data
				tpred = get_tpred_gp(mjd_or_phase,phase)#Get GP time grid
				if len(bright)>1:
					try:#Interpolation
						GPfit     = GP_1D_squaredexp(phase,bright,brighterr,x_pred=tpred,tau_guess=10)
						if bright_mode=='mag' and GPfit is None:#If fit fails try with larger errors
							try:
								print (f"Try Inflating errors by 10% for SN{lc.meta['SNID']} {flt}-band")
								GPfit     = GP_1D_squaredexp(phase,bright,brighterr*1.1,x_pred = tpred,tau_guess=10)
								if GPfit is None:	print (f"Inflating errors by 10% did not work for SN{lc.meta['SNID']} {flt}-band")
								else:				print (f"Inflating errors success for SN{lc.meta['SNID']} {flt}-band")
							except Exception as e:
								print (e)
								FIT[flt]  = None
						#If we interpolate in flux space, we want to transform to mags by transforming each gp flux draw,
						#NOT finding mean flux then transforming to mags (Jensens Inequality)
						if self.choices['bright_mode']=='flux':#If interpolating flux
							f_samps      = GPfit.gp.sample_conditional(bright, GPfit.x, size=self.choices['Ngpdraws'])
							m_samps      = -2.5*np.log10(f_samps)+27.5
							m_mean_curve = np.average(m_samps,axis=0) ; m_std_curve  = np.std(m_samps,axis=0)

							GPfit.y    = m_mean_curve
							GPfit.yerr = m_std_curve
							GPfit.df   = pd.DataFrame({'x':GPfit.x,'y':GPfit.y,'yerr':GPfit.yerr,'f':np.average(f_samps,axis=0),'ferr':np.std(f_samps,axis=0)})

						FIT[flt]  = copy.deepcopy(GPfit)
					except ValueError:
						print (f"{self.lc.meta['SNID']} GP Interpolation Failed for Band: {flt}")
						FIT[flt]  = None
				else:
					print (f"{self.lc.meta['SNID']} GP Interpolation not attempted because {len(bright)} points in Band: {flt}")
					FIT[flt] = None

		elif self.choices['mags_method']=='2DGP':
			lc  = self.lc
			phase,bright,brighterr = get_time_lc_arrays(lc,mjd_or_phase=mjd_or_phase,flux_or_mag=self.choices['bright_mode'])#Get data
			tpred = get_tpred_gp(mjd_or_phase,phase)#Get GP time grid
			#Get lam to flt mapping
			lam_to_flt  = dict(zip(lc.meta['lams'],lc.meta['flts']))
			flt_to_lam  = dict(zip(lc.meta['flts'],lc.meta['lams']))
			wavelengths = np.array([flt_to_lam[flt] for flt in lc['flt']])
			lambdaC     = np.asarray(list(lam_to_flt.keys()))
			#Get 2DGP Interpolation
			gFIT   = GP_2D_Matern(phase,bright,brighterr,lambdaC,wavelengths,x_pred=tpred)

			#If we interpolate in flux space, we want to transform to mags by transforming each gp flux draw,
			#NOT finding mean flux then transforming to mags (Jensens Inequality)
			if self.choices['bright_mode']=='flux':
				Ngridt      = copy.deepcopy(len(tpred))
				new_tpred   = np.vstack([np.hstack((tpred for _ in range(len(lambdaC)))),np.array([ll for ll in lambdaC for _ in range(len(tpred))])]).T
				f_samps_all = gFIT[wavelengths[0]].gp.sample_conditional(bright, new_tpred, size=self.choices['Ngpdraws'])
				for ll in gFIT:
					GPfit        = copy.deepcopy(gFIT[ll])
					il           = list(lambdaC).index(ll)
					f_samps      = f_samps_all[:,il*Ngridt:(il+1)*Ngridt]
					m_samps      = -2.5*np.log10(f_samps)+27.5
					m_mean_curve = np.average(m_samps,axis=0) ; m_std_curve  = np.std(m_samps,axis=0)

					GPfit.y      = m_mean_curve
					GPfit.yerr   = m_std_curve
					GPfit.df     = pd.DataFrame({'x':GPfit.x,'y':GPfit.y,'yerr':GPfit.yerr,'f':np.average(f_samps,axis=0),'ferr':np.std(f_samps,axis=0)})

					gFIT[ll]     = copy.deepcopy(GPfit)

			for ll in gFIT:
				FIT[lam_to_flt[ll]] = gFIT[ll]
		#Set attribute
		self.FIT = FIT

	def extract_interpolations(self):
		"""
		Extract Interpolations

		Method to extract FIT dictionary interpolations in passbands in pblist at times in tilist

		End Product(s)
		----------
		df_m, df_m_extra : pandas dfs of interpolations
			respectively, GP evaluated in pblist at tilist, and evaluated at Extra_Features e.g. {15:['B']}
		"""
		#Load interpolations, pblist and tilist
		FIT  = self.FIT ; pblist = self.choices['pblist'] ; tilist = self.choices['tilist']
		#Intialise empty matrix of interpolations and measurement errors
		mags = np.empty((len(tilist),len(pblist),2))
		for ipb,pb in enumerate(pblist):
			for iti,ti in enumerate(tilist):
				if FIT[pb] is not None:
					gp_index        = np.where(FIT[pb].x==ti)[0][0]
					y,yerr          = FIT[pb].y[gp_index],FIT[pb].yerr[gp_index]
					mags[iti,ipb,0] = y
					mags[iti,ipb,1] = yerr
				else:
					if iti==0: print (f"{FIT['lc'].meta['SNID']} GP fit is None in Band: {pb}")
					mags[iti,ipb,0] = np.nan
					mags[iti,ipb,1] = np.nan
		#Transform to pandas df
		df_m = pd.DataFrame(np.concatenate((mags[:,:,0],mags[:,:,1]),axis=1),index=tilist,columns=pblist+[pb+self.choices['errstr'] for pb in pblist])

		#Do the same for df_extra
		#Slightly different because each chosen phase may have different set of columns==passbands, therefore 1 df for each ti in Extra_Features
		df_extra = {}
		tilist_extra = list(self.choices['Extra_Features'].values())
		for ti,pbmini in self.choices['Extra_Features'].items():
			mags = np.empty((len(pbmini),2))
			for ipb,pb in enumerate(pbmini):
				if FIT[pb] is not None:
					gp_index        = np.where(FIT[pb].x==ti)[0][0]
					y,yerr          = FIT[pb].y[gp_index],FIT[pb].yerr[gp_index]
					mags[ipb,0] = y
					mags[ipb,1] = yerr
				else:
					if iti==0: print (f"{FIT['lc'].meta['SNID']} GP fit is None in Band: {pb}")
					mags[ipb,0] = np.nan
					mags[ipb,1] = np.nan

				df_extra[ti] = pd.DataFrame(np.hstack((mags[:,0],mags[:,1])).reshape(1,len(pbmini*2)),index=[ti],columns=pbmini+[pb+self.choices['errstr'] for pb in pbmini])#index=[ti]

		#Set attributes
		self.df_m       = df_m
		self.df_m_extra = df_extra

	def update_DF_M(self, DF_M, sn):
		"""
		Update DF_M

		Takes in sn-specific interpolations dataframes, df_m and df_m_extra, and appends to sample level DF_M

		Parameters
		----------
		DF_M : dict of pandas df
			sample level collection of GP interpolations

		sn : str
			name of sn being appended to sample object

		Returns
		----------
		DF_M with sn interpolations appended
		"""
		for ti in self.choices['tilist']:
			DF_M[ti].loc[sn] = self.df_m.loc[ti]
		for ti,pbmini in self.choices['Extra_Features'].items():
			DF_M['extra'][ti].loc[sn] = self.df_m_extra[ti].loc[ti]
		return DF_M

	def plot_lc(self, plotter, mjd_or_phase='phase', bright_mode=None):
		"""
		Plot GP Interpolation

		Function to plot phase light curves (flux or mag) including both data and GP interpolation

		Parameters
		----------
		plotter : class object
			for plotting, contains e.g. fontsize, choice to show/save etc.

		mjd_or_phase : str (optional; default='phase')
			choice of time data to plot

		bright_mode : str or None (optional; default=None)
			choice of LC data to plot, if None, uses self.choices['bright_mode']

		End product(s)
		----------
		Saves and/or shows plot
		"""
		#Data/Model to plot
		FIT  = self.FIT
		lc   = FIT['lc']
		flts = lc.meta['flts']

		#Get bright_mode
		if bright_mode is None: bright_mode = self.choices['bright_mode']
		legend = {'flux':True,'mag':False}[bright_mode]

		#Add this functionality later
		reasons = None

		#Initialise plot
		pl.figure(figsize=(8,6))
		pl.title(f"{self.choices['mags_method']} Interpolation: SN {lc.meta['SNID']}"+(reasons is not None)*f"\n{reasons}",fontsize=plotter.FS)

		#Plot each band
		for iif,flt in enumerate(flts):
			color = f'C{iif}'; const = 1.5*(len(flts)-1-iif)
			#Get Data
			phase,bright,brighterr = get_time_lc_arrays(get_lcf(lc,flt),mjd_or_phase=mjd_or_phase,flux_or_mag=bright_mode)
			#Plot Data
			pl.errorbar(phase, bright+const*(bright_mode=='mag'), yerr=brighterr, marker=".", markersize=10, linestyle="", color=color, label=legend*flt)
			#Plot GP interpolation
			if len(bright)>1 and self.choices['mags_method']=='1DGP' or self.choices['mags_method']=='2DGP':
				GPfit = copy.deepcopy(FIT[flt])
				if bright_mode=='flux':	GPfit.y = GPfit.df['f'].values ; GPfit.yerr = GPfit.df['ferr'].values
				if GPfit is not None:
					if not legend: 	#Plot passband at final phase
						pl.annotate(flt,xy=(GPfit.x[-1]+0.5,GPfit.y[-1]+const*(bright_mode=='mag')),weight='bold',fontsize=15)
					if bright_mode=='mag':#Offset each LC for mag plot
						GPfit.y += const
					pl.fill_between(GPfit.x, GPfit.y - GPfit.yerr, GPfit.y + GPfit.yerr,color="k", alpha=0.2)
					pl.plot(GPfit.x, GPfit.y, "k", lw=1.5, alpha=0.5)
				else:
					print (f'No successful GP fit to plot in Band: {flt}')
			elif len(bright)==1 and self.choices['mags_method']=='1DGP':
				if not legend:
					pl.annotate(flt,xy=(phase+0.5, bright+const*(bright_mode=='mag')),weight='bold',fontsize=15)
			else:
				print (f'Only 1 data point for Band: {flt}; therefore, no GP plotted')

		#Finishing touches
		xlabel   = {'phase':'rest-frame phase (days)', 'mjd':'MJD'}[mjd_or_phase]
		ylabel   = bright_mode+' + const.'*(bright_mode=='mag')
		savename = f"{plotter.plotpath}GPinterp_plots/SN{lc.meta['SNID']}_{self.choices['mags_method']}.pdf"
		plotter.finish_plot(xlabel,ylabel,savename=savename,invert=bright_mode=='mag',legend=legend)
