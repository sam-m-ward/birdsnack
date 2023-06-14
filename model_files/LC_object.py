"""
LC Object Module

Module containing LCObj class and additional functions
Used to take snana lc, and perform various transformations including getting Tmax with GP, and interpolating with GP


Contains:
--------------------
function1(inputs):
	description

LCObj class:
	inputs: snana lc

	Methods are:
		method1
		method2
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

from snana_functions import *
from GP_functions import *


class LCObj:

	def __init__(self,lc,choices):
		self.lc      = lc
		self.choices = choices

	def get_Tmax(self,Tmaxchoicestr='Tmax_GP_restframe',return_samps=False):
		"""
		Get Tmax Method

		Method uses GPs (1D or 2D) to interpolate data and estimate Tmax
		Uses model maximum brightness as initial estimate, then cuts data to reside within phase range,
		And finally, draws Ngpdraws samples, and records sample mean and std of these.

		Parameters
		----------
		Tmaxchoicestr : str (optional; default='Tmax_GP_restframe')
			string used to identify GP Tmax estimate in lc.meta

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
		def trim_lc_on_phase(bright_mode,GPfit,zhelio,phasemin,phasemax,data):
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
							mjd, bright, brighterr = trim_lc_on_phase(bright_mode,GPfit,self.lc.meta['REDSHIFT_HELIO'],self.choices['phasemin'],self.choices['phasemax'],[mjd,bright,brighterr])
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
				lambdaC     = np.asarray(lam_to_flt.keys())
				#Get 2DGP Interpolation
				gFIT   = GP_2D_Matern(mjd,bright,brighterr,lambdaC,wavelengths)#,x_pred = x_pred)
				#Extract the reference band component for Tmax estimation
				lamref = flt_to_lam[self.choices['ref_band']]
				GPfit  = gFIT[lamref]
				if GPfit is not None:#Re-fit data within phase range
					mjd, bright, brighterr = trim_lc_on_phase(bright_mode,GPfit,lc.meta['REDSHIFT_HELIO'],self.choices['phasemin'],self.choices['phasemax'],[mjd,bright,brighterr])
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
					new_xpred   = np.vstack([np.hstack((tmax_grid for _ in range(len(lambdaC)))),np.array([lam for lam in lambdaC for _ in range(len(tmax_grid))])]).T
					f_samps_all = GPfit.gp.sample_conditional(bright, new_xpred, size=int(Ngpdraws/10))
					il          = list(lambdaC).index(lamref) ; Ngridt = copy.deepcopy(len(tmax_grid))
					f_samps     = f_samps_all[:,il*Ngridt:(il+1)*Ngridt]
				#Samples of tmax
				tmaxs     = np.array([tmax_grid[i] for i in np.argmax(f_samps-2*f_samps*(bright_mode=='mag'),axis=1)])
				#Tmax estimate and uncertainty is sample average and std. of tmax samples
				Tmax = np.average(tmaxs) ; Tmax_std = np.std(tmaxs)
			else:
				Tmax, Tmax_std = None, None
				tmax_grid,f_samps,tmaxs = None, None, None

			#print (Tmax, Tmax_std)
			self.lc.meta[Tmaxchoicestr] = Tmax
			self.lc.meta[f'{Tmaxchoicestr}_std'] = Tmax_std
		else:
			#print ('Tmax already estimated')
			pass

		if return_samps: return tmax_grid,f_samps,tmaxs





def get_phase_method(lc,phaselims,Tmaxchoicestr,interpflts,fupperbool):
	"""
	Get Phase Method

	Called upon by RVGP class to ensure phase column is computed, and lc is cut on phaselims

	Parameters
	----------
	lc: :py:class:`astropy.table.Table`
		light curve object

	phaselims: list
		phasemin,phasemax,dp = phaselims[:]

	Tmaxchoicestr: str
			lc.meta entry name for Tmax, e.g. 'Tmax_GP_restframe' or 'Tmax_snpy_fitted'

	interpflts: str
		string of filters we allow e.g. 2DGP or SNooPy to use when interpolating or fitting, e.g. 'BVrRiIYJH'
		either 'all' or string of filts e.g. 'BVriYJH'

	fupperbool: bool
		if True, then lower and upper case belong to same family

	Returns
	----------
	lc: :py:class:`astropy.table.Table`
		light curve object with phase column cut on phaselims
	"""
	phasemin,phasemax,dp = phaselims[:]
	#Create phase column and cut data on phasemin and phasemax
	lc            = create_phase_column(lc,lc.meta[Tmaxchoicestr],phasemin=phasemin,phasemax=phasemax)
	#Given data has been cut, there is potential that some bands are removed completely, e.g. Y_WIRC
	lc            = update_lcmetaflts(lc,interpflts,fupperbool)
	return lc
