# -*- coding: UTF-8 -*-
"""
BayeSN SED Model. Defines a class which allows one to fit or simulate from the
BayeSN Optical+NIR SED model.
"""

import os
import copy
import time
import tempfile
import numpy as np
import scipy.sparse
from astropy.table import Table, unique
from astropy.cosmology import FlatLambdaCDM
import extinction
import cmdstanpy

from . import io
from . import passband
from . import spline_utils
from . import utils

class SEDmodel(object):
	"""
	BayeSN-SED Model

	Class which imports a BayeSN model, and allows one to fit or simulate
	Type Ia supernovae based on this model.

	Parameters
	----------
	model : str, optional
		Can be either a pre-defined BayeSN model name (see table below), or
		a path to directory containing a set of .txt files from which a
		valid model can be constructed. Currently implemented default models
		are listed below - default is M20. See README in `BayeSNmodel/model_files`
		for more info.
		  ``'M20'`` | Mandel+20 BayeSN model (arXiv:2008.07538). Covers
					|   rest wavelength range of 3000-18500A (BVRIYJH). No
					|   treatment of host mass effects. Global RV assumed.
					|   Trained on low-z Avelino+19 (ApJ, 887, 106)
					|   compilation of CfA, CSP and others.
		  ``'T21'`` | Thorp+21 No-Split BayeSN model (arXiv:2102:05678). Covers
					|   rest wavelength range of 3500-9500A (griz). No
					|   treatment of host mass effects. Global RV assumed.
					|   Trained on Foundation DR1 (Foley+18, Jones+19).
	fiducial_cosmology :  dict, optional
		Dictionary containg kwargs ``{H0, Om0}`` for initialising a
		:py:class:`astropy.cosmology.FlatLambdaCDM` instance. Defaults to
		Riess+16 (ApJ, 826, 56) cosmology ``{H0:73.24, "Om0":0.28}``.
	compile : bool, optional
		Decides whether to precompile Stan code necessary for fitting.
		If False, Stan code can be compiled later with a call to
		`compile_stan_model`. If True, compiles Stan model (takes around
		30 seconds)
	fix_tmax : bool, optional
		If precompiling Stan model, decides whether to load model where tmax
		is fixed (defaults to True).

	Attributes
	----------
	params : dict
		Dictionary containing BayeSN model parameters.
	zpt : float
		SNANA zero point which is to be assumed
	t_k : :py:class:`numpy.array`
		Array of time knots which the model is defined at
	l_k : :py:class:`numpy.array`
		Array of wavelength knots which the model is defined at
	knots : :py:class:`numpy.array`
		Cartesian product of t_k and l_k
	dl_int : float
		Wavelength pacing of the Hsiao template (defaults to 10A)
	hsiao : dict
		Dictionary containing specification of the Hsiao template from
		Hsiao+07 (ApJ, 663, 1187)
	passbands : dict
		Dictionary of all passbands with available specification. Indexed
		by filter name.
	cosmo : :py:class:`astropy.cosmology.FlatLambdaCDM`
		:py:class:`astropy.cosmology.FlatLambdaCDM` instance defining the
		fiducial cosmology which the model was trained using.
	stan_model : None or :py:class:`cmdstanpy.CmdStanModel`
		Compiled Stan model for fitting the BayeSN photometric distance
		model to a supernova.
	fixed_tmax : bool
		Indicates whether the currently compiled Stan model assumes a fixed
		time of maximum.

	Returns
	-------
	out : :py:class:`BayeSNmodel.bayesn_model.SEDmodel` instance

	Notes
	-----
		Attributes should not be messed with unless you know what you're
		doing. The same is probably true for any method which begins with
		`load_`.
	"""
	def __init__(self, model="M20", fiducial_cosmology={"H0":73.24, "Om0":0.28}, compile=False, fix_tmax=True):
		path_to_try = io.get_pkgfile("model_files/")
		if os.path.exists(path_to_try + model + "_model/"):
			print("Loading " + model + " model...")
			self.load_model_parameters(path_to_try + model + "_model/")
			print("Loaded successfully!")
		else:
			print("Non standard model name - trying as path...")
			if os.path.exists(model):
				self.load_model_parameters(model)
				print("Loaded succesfully!")
			else:
				raise IOError("Could not find model at: " + model)

		self.load_hsiao_template()
		self.load_passbands()

		self.cosmo = FlatLambdaCDM(**fiducial_cosmology)

		if compile is True:
			self.compile_stan_model(fix_tmax)
		else:
			self.stan_model = None
		self.fixed_tmax = fix_tmax

	def load_model_parameters(self, path):
		"""
		Loads parameters of a BayeSN model.

		Goes to the provided path and loads in the files which define a
		BayeSN model. Tampers with the SEDmodel internals, so you should only
		use this if you're sure it's the right thing to do. You should be
		especially wary of providing a path to something other than one
		of the provided models.

		Parameters
		----------
		path : str
			Path to a directory containing the necessary files for constructing
			a BayeSN model (see notes).

		Notes
		-----
		A valid BayeSN model can be specified by 6 .txt files:
		  l_knots           | A single column of wavelength knots
		  tau_knots         | A single colun of phase knots
		  W0                | A grid with len(l_knots) rows and len(tau_knots)
							|   columns which defines the W0 component.
		  W1                | A grid with len(l_knots) rows and len(tau_knots)
							|   columns which defines the W1 component.
		  L_Sigma_epsilon   | A square grid with side of length
							|   len(tau_knots)*len(l_knots) containing the
							|   Cholesky decomposition of the Sigma_epsilon matrix.
		  M0_sigma0_RV_tauA | A list of the global M0, sigma0, RV, and tauA, may
							|   constain multiple versions in the case of a
							|   partial-split style model. In this case, a fifth
							|   column should list the split mass.

		"""
		self.params = {}
		for param in ["l_knots", "tau_knots", "W0", "W1", "M0_sigma0_RV_tauA", "L_Sigma_epsilon"]:
			try:
				p = np.genfromtxt(path + param + ".txt")
			except OSError as e:
				raise Exception("Missing {}.txt file in provided model.".format(param)) from e
			else:
				if param == "M0_sigma0_RV_tauA":
					if p.shape == (4,):
						for n, name in enumerate(param.split("_")):
							self.params[name] = p[n]
						self.M_split = False
					elif p.shape == (2,5):
						for n, name in enumerate(param.split("_")):
							self.params[name + "_HM"] = p[0,n]
							self.params[name + "_LM"] = p[1,n]
						if p[0,-1] == p[1,-1]:
							self.M_split = p[0,-1]
						else:
							raise ValueError("Ambiguous split mass given in M0_sigma0_RV_tauA.txt file!")
					else:
						raise ValueError("Unexpected format of M0_sigma0_RV_tauA.txt file!")
				elif param == "tau_knots":
					self.t_k = p
				elif param == "l_knots":
					self.l_k = p
				else:
					self.params[param] = p
					if param == "L_Sigma_epsilon":
						self.params["Sigma_epsilon"] = np.matmul(p, p.T)

		self.zpt = 27.5

		self.knots = spline_utils.cartesian_prod(self.t_k, self.l_k)

	def load_hsiao_template(self, dl_int=10):
		"""
		Loads the Hsiao template from the internal HDF5 file.

		Stores the template, evaluated at the chosen wavelength spacing,
		as an attribute of `SEDmodel`.

		Parameters
		----------
		dl_int : float, optional
			Wavelength spacing in Angstroms to load the template at.
			Defaults to 10A.
		"""
		t_hsiao, l_hsiao, f_hsiao = io.read_model_grid()
		l_int = np.arange(min(l_hsiao),max(l_hsiao)+dl_int,dl_int,dtype=float)
		f_hsiao = f_hsiao[:,np.in1d(l_hsiao,l_int)]
		l_hsiao = l_int

		self.dl = dl_int
		self.hsiao = {"t": t_hsiao, "l": l_hsiao, "f": f_hsiao}

	def load_passbands(self, pbfile=None, cutoff_factor=0.002):
		"""
		Loads the passbands we have available.

		These loaded by default from `SNmodel_pb_obsmode_map.txt`.

		A different set of passbands can be loaded by recalling this
		function with a different ``pbfile`` argument.

		Stores the loaded passbands internally to the `SEDmodel` class
		instance.

		Parameters
		----------
		pbfile : str, optional
			Full path to a file specifying the passbands to be used.
			Propagates through to :py:func:`BayeSNmodel.passband.get_pbmodel`.
			If not provided, defaults to `SNmodel_pb_obsmode_map.txt`.
		cutoff_factor : float, optional
			Fraction of max throughput to truncate at. Passband is truncated
			beyond the last wavelength where the throughput is > `cutoff_factor`
			times the maximum. Default is 0.002 (0.2% of max throughput).

		See Also
		--------
		:py:func:`BayeSNmodel.passband.get_pbmodel`
		"""
		if pbfile is None:
			pb_file = io.get_pkgfile("SNmodel_pb_obsmode_map.txt")
		else:
			pb_file = pbfile
		pb_list = np.genfromtxt(pb_file, usecols=0, dtype=str, skip_header=1)

		self.passbands = passband.get_pbmodel(pb_list, pbfile=pbfile, cutoff_factor=cutoff_factor)

	def compile_stan_model(self, fix_tmax=True):
		"""
		Compiles the Stan code needed for fitting.

		Loads a version of the BayeSN model configuered
		for phototmetric distance estimation for a single supernova.
		Can load a version which treats the time of maximum as known and
		fixed, or which fits for time of maximum.

		Takes some time.

		Parameters
		----------
		fix_tmax : bool, optional
			If True, compiles a photometric distance model which
			assumes the time of maximum is known and fixed. Else,
			loads a version which enables one to fit for tmax.
		"""
		path_to_try = io.get_pkgfile("stan_files/")
		if fix_tmax is True:
			self.stan_model = cmdstanpy.CmdStanModel(stan_file=path_to_try + "photometric_distance_model.stan")
		else:
			self.stan_model = cmdstanpy.CmdStanModel(stan_file=path_to_try + "photometric_distance_model_free_tmax.stan")
		self.fixed_tmax = fix_tmax

	def fit_supernova(self, lc, z=None, z_cmb=None, ebv_mw=None, logM=None, tmax=None, fix_tmax=True, phase_shift=5, snr_cut=3, n_chains=4, n_workers=None, n_warmup=250, n_sampling=250, show_progress=True, show_console=False, cmdstan_output_dir=None, fit_mode="hmc"):
		"""
		Fits the light curves of a single supernova.

		Carries out necessary precomputation, and runs the Stan code required
		to fit the BayeSN photometric distance model to a provided supernova.

		This can take some time depending on the number of observations
		available.

		Parameters
		----------
		lc : :py:class:`astropy.table.Table`
			Astropy table containing light curve information. This must contain
			`phase` (and/or `mjd`), `fluxcal`, `fluxcalerr` and `flt` columns.
			It is recommended that this table be initialised from an SNANA format
			light curve file using :py:func:`BayeSNmodel.io.read_snana_lcfile`, since
			this should retrieve correctly formatted columns, and well populated
			metadata.
		z : float or None, optional
			Heliocentric redshift used for redshifting SED. Doesn't need to
			be specified explicitly if contained in `lc.meta` as ``'REDSHIFT_HELIO'``.
		z_cmb : float or None, optional
			CMB redshift, used to compute prior mean distance modulus. Needn't
			be given if sufficient information given in `lc.meta` as
			``'REDSHIFT_CMB'``, ``'REDSHIFT_FINAL'``, or ``'MU'``.
		ebv_mw : float or None, optional
			Milky Way reddening. Needn't be given if included in `lc.meta` as
			``'MWEBV'``.
		logM : float or None, optional
			log_10 of host galaxy stellar mass (only needed when using a model with
			a host galaxy mass split).
		tmax : float or None, optional
			Time of B-band maximum in MJD. If `fix_tmax` is True, this will be
			fixed when carrying out light curve fits. If `fix_tmax` is False,
			this will be taken as the centre of a uniform prior on tmax,
			tmax ~ U(`tmax` - `phase_shift*(1+z)`, `tmax` + `phase_shift*(1+z)`).
			Needn't be given if included in `lc.meta` as ``'PEAKMJD'`` or
			``'SEARCH_PEAKMJD'``.
		fix_tmax : bool, optional
			If True, assumes tmax is known, and fixes this to `tmax`, or the
			value given by ``'SEARCH_PEAKMJD'`` or ``'PEAKMJD'``
		phase_shift : float, optional
			Maximum positive or negative shift in rest-frame time of maximum
			permitted if this is being fit for. Defaults to 5 rest-frame days.
		snr_cut : float, optional
			Minimum S/N cut to apply to light curves before fitting. Default=3.
		n_chains : int, optional
			Number of independent HMC chains to run. Default is 4.
		n_workers : int or None, optional
			Number of cores to parallelise the MCMC over. Will only paralellise
			up to the chain level (i.e. `n_workers` > `n_chains` has no effect).
			Default is `n_workers == n_chains`, i.e. one chain per core.
		n_warmup : int, optional
			Number of HMC warmup iterations to use. Default is 250 per chain.
		n_sampling : int, optional
			Number of post-warmup HMC iterations. Default is 250 per chain.
		show_progress : bool, optional
			Display a progress bar for the HMC (need tqdm). Default is True
		show_console : bool, optional
			Display cmdstan console output. Default is False
		cmdstan_output_dir : str or None, optional
			Specify a permanent directory where cmdstan's outputs will be
			written. If not provided, they will be written to a tempfile
			then dumped once the run has finished. If the provided path is
			invalid, the default behaviour of writing to a tempfile will
			be reverted to.
		fit_mode : str, optional
			Choose the fitting mode. Default is ``'hmc'``, which samples
			from the posterior using Stan's Hamiltonian Monte Carlo. Can
			also select ``'opt'``, which produces a MAP estimate using
			Stan's L-BFGS optimizer.


		Returns
		-------
		fit : :py:class:`cmdstanpy.CmdStanMCMC` or :py:class:`cmdstanpy.CmdStanMLE`
			If ``fit_mode='hmc'``, :py:class:`cmdstanpy.CmdStanMCMC`
			instance containing complete output of HMC run.
			If ``fit_mode='opt'``, :py:class:`cmdstanpy.CmdStanMLE`
			instance returned by :py:func:`cmdstanpy.CmdStanModel.optimize`.
		results : dict
			If ``fit_mode='hmc'``, dictionary containing post-warmup HMC
			samples, with all chains stacked.
			If ``fit_mode='opt'``, dictionary containing MAP estimates of parameters.
			Indexed by name of parameter ``{"lp__", "mu", "delM", "AV",
			"theta", "epsilon",	["tmax"], ["dtmax"]}``.
		summary : :py:class:`pandas.DataFrame`
			Dataframe of summary statistics

		Notes
		-----
		In the output samples dict, each entry will be a :py:class:`numpy.array`
		whose first dimension corresponds to the HMC iterations. Subsequent
		dimensions are inheritted by parameter dimensions. The parameter
		names which index the dictionary are:
		  ``'lp__'``     | Log posterior + const.
		  ``'mu'``       | Distance modulus
		  ``'delM'``     | Gray offset, delta M_s
		  ``'AV'``       | Host galaxy dust extinction
		  ``'theta'``    | Scalar coefficient of W1 SED component
		  ``'epsilon'``  | Matrix of residual perturbations
		  ``'tmax'``     | Time of maximum in MJD (only if `fix_tmax` is False)
		  ``'dtmax'``    | Shift in tmax compared to input
		"""
		#Check the provided fitting mode is valid
		if fit_mode not in ["hmc", "opt"]:
			raise ValueError("Only available fitting modes are 'hmc' and 'opt'! Please request one of these!")
		else:
			if fit_mode == "opt":
				n_chains = 1

		#Check light curve in recognised format
		lc = copy.deepcopy(lc)
		if not isinstance(lc, Table):
			raise TypeError("Light curve must be provided as Astropy table with phase (or mjd), fluxcal, fluxcalerr, flt columns.")

		#Check for heliocentric redshift
		if z is None:
			try:
				z = lc.meta["REDSHIFT_HELIO"]
			except KeyError as e:
				raise Exception("Redshift not provided or found in light curve metadata.") from e
		#Check for MW reddening
		if ebv_mw is None:
			try:
				ebv_mw = lc.meta["MWEBV"]
			except KeyError as e:
				raise Exception("Milky Way reddening not provided or found in light curve metadata.") from e
		#Check for hostmass (if necessary)
		if self.M_split is not False and logM is None:
			try:
				logM = lc.meta["HOSTGAL_LOGMASS"]
			except KeyError as e:
				raise Exception("Host galaxy mass not provided or found in light curve metadata") from e
		#Check for tmax
		if tmax is None:
			try:
				tmax = lc.meta["SEARCH_PEAKMJD"]
			except KeyError:
				try:
					tmax = lc.meta["PEAKMJD"]
				except KeyError as e:
					raise Exception("Time of maximum not provided or found in light curve metadata.") from e
		#Check for CMB redshift or external distance estimate
		#External distance estimate (provided as "MU" in metadata) takes precedence
		#over redshift metadata. Something passed direct to the function has ultimate
		#precedence
		if z_cmb is None:
			try:
				mu = lc.meta["MU"]
			except KeyError:
				try:
					z_cmb = lc.meta["REDSHIFT_CMB"]
				except KeyError:
					try:
						z_cmb = lc.meta["REDSHIFT_FINAL"]
					except KeyError as e:
						raise Exception("CMB redshift/distance estimate not provided or found in light curve metadata.") from e
				mu = self.cosmo.distmod(z_cmb).value
		else:
			mu = self.cosmo.distmod(z_cmb).value

		#Population params which might be mass split
		if self.M_split is not False:
			if logM >= self.M_split:
				RV = self.params["RV_HM"]
				tauA = self.params["tauA_HM"]
				M0 = self.params["M0_HM"]
				sigma0 = self.params["sigma0_HM"]
			else:
				RV = self.params["RV_LM"]
				tauA = self.params["tauA_LM"]
				M0 = self.params["M0_LM"]
				sigma0 = self.params["sigma0_LM"]
		else:
			RV = self.params["RV"]
			tauA = self.params["tauA"]
			M0 = self.params["M0"]
			sigma0 = self.params["sigma0"]
		#Turn of sigma0 if optimizing to prevent undesireable behaviour
		if fit_mode == "opt":
			sigma0 = 1e-6

		#Computes phase
		if "phase" not in lc.keys():
			if "mjd" in lc.keys():
				lc["phase"] = (lc["mjd"] - tmax)/(1 + z)
			else:
				raise ValueError("Light curve table doesn't have phase or mjd column")
		if fix_tmax is False and "mjd" not in lc.keys():
			lc["mjd"] = lc["phase"]*(1 + z) + tmax

		#Phase and S/N cut
		mask = (lc["phase"] > min(self.t_k))*(lc["phase"] < max(self.t_k))*(lc["fluxcal"]/lc["fluxcalerr"] > snr_cut)
		lc = lc[mask]

		#Matrix needed for wavelength spline
		KD_l = spline_utils.invKD_irr(self.l_k)

		#Dictionary where we'll store precomputed things that are
		#only dependent on passband (and not observation time)
		flt_dep_pieces = {}

		#Cycles through unique filters
		for flt in unique(lc, keys="flt")["flt"]:
			if flt not in self.passbands.keys():
				raise ValueError("Unrecognised passband: " + flt)

			#Interpolates transmision curve to appropriate wavelengths
			#l_r is the rest frame wavelength vector this filter will probe
			l_r, tp = utils.interpolate_passband(self.hsiao["l"], self.passbands[flt][0]["wave"], self.passbands[flt][0]["norm_throughput"], z)

			if np.any(l_r < min(self.l_k)) or np.any(l_r > max(self.l_k)):
				print("WARNING: " + flt + " filter extends beyond model coverage - discarding.")
				lc = lc[lc["flt"] != flt]
				continue

			#Observer frame wavelengths probed by filter
			l_o = l_r*(1+z)
			#MW reddening curve in this filter (in flux units)
			R_MW = utils.extinction_mw(l_o, ebv_mw)

			#Rectangle rule for this passband (also includes MW reddening)
			#Forms row of large H matrix passed to Stan
			H_row = np.array([tp*l_o*self.dl*R_MW])

			#Coefficients for interpolating SED to wavelengths in this filter
			J_l_block = spline_utils.spline_coeffs_irr(l_r, self.l_k, KD_l, allow_extrap=False)

			#Adds all these pieces to the dict
			flt_dep_pieces[flt] = {"l_o": l_o, "l_r": l_r, "H": H_row, "J_l": J_l_block}

			#If fitting for tmax, also need to save some parts of the Hsiao template
			#In this case, we have to interpolate this within Stan
			if fix_tmax is False:
				f_hsiao_block = self.hsiao["f"][:, np.in1d(self.hsiao["l"], l_r)]
				flt_dep_pieces[flt]["f_hsiao"] = f_hsiao_block

			#A(lambda)/AV host extinction curve in this filter (in magnitudes)
			R_host_block = extinction.fitzpatrick99(l_r, 1.0, RV)
			flt_dep_pieces[flt]["Alam_AV"] = R_host_block

		#Matrix needed for time spline
		KD_t = spline_utils.invKD_irr(self.t_k)
		#Coefficients for interpolating SED in time
		#If we're fitting for tmax, this has to be computed in Stan
		if fix_tmax is True:
			J_t = spline_utils.spline_coeffs_irr(lc["phase"], self.t_k, KD_t).T

		#Parts of matrix for passband integral + Milky Way reddening
		H = []
		#Parts of matrix for wavelength interpolation
		J_l = []

		#Empty vector for applying host reddening
		R_host = np.empty((0,))

		#Empty container for Hsiao template
		if fix_tmax is True:
			S0 = np.empty((0,))
		else:
			f_indices = np.ones((1,), dtype=int)
			f_hsiao = []

		#Cycle through observations, doing all precomputation,
		#and preparing different matrices
		for t, flt in lc.iterrows("phase", "flt"):
			if flt not in flt_dep_pieces.keys():
				continue
			#Hsiao template (or parts to be interpolated in time)
			if fix_tmax is True:
				S0 = np.append(S0, utils.interpolate_hsiao(t, flt_dep_pieces[flt]["l_r"], self.hsiao["t"], self.hsiao["l"], self.hsiao["f"]))
			else:
				f_indices = np.append(f_indices, len(flt_dep_pieces[flt]["l_r"]))
				f_hsiao.append(flt_dep_pieces[flt]["f_hsiao"])

			#Host reddening for this observation
			R_host = np.append(R_host, flt_dep_pieces[flt]["Alam_AV"])

			#Passband integral + MW reddening for this observation
			H.append(flt_dep_pieces[flt]["H"])
			#Wavelength interpolation for this observation
			J_l.append(flt_dep_pieces[flt]["J_l"])

		#Sparse matrix for passband integral
		H = scipy.sparse.block_diag(H, format="csr")
		#Sparse matrix for wavelength interpolation
		J_l = scipy.sparse.block_diag(J_l, format = "csr")
		if fix_tmax is False:
			f_hsiao = np.hstack(f_hsiao)

		#Data to be given to Stan
		stan_data = {"Nknots": len(self.t_k)*len(self.l_k), "Ntknots": len(self.t_k), "Nlknots": len(self.l_k), "tk": self.knots[:,0], "lk": self.knots[:,1], "Jl_nnz": J_l.nnz, "Jl_w": J_l.data, "Jl_v": J_l.indices+1, "Jl_u": J_l.indptr+1, "Alam_AV": R_host, "H_nnz": H.nnz, "H_w": H.data, "H_v": H.indices+1, "H_u": H.indptr+1, "W0": self.params["W0"], "W1": self.params["W1"], "L_Sigma": self.params["L_Sigma_epsilon"], "tauA": tauA, "M0": M0, "sigma0": sigma0, "Nobs": len(lc), "Dint": H.shape[1], "fobs": lc["fluxcal"].data, "fobserr": lc["fluxcalerr"].data, "muhat": mu, "muhaterr": 100, "ZPT": self.zpt}

		#Extra bits which depend on if tmax is being fitted
		if fix_tmax is True:
			stan_data.update({"tobs": lc["phase"].data, "S0": S0, "Jt": J_t})
		else:
			stan_data.update({"tobs": lc["mjd"].data, "zhel": z, "tmax_guess": tmax, "phase_shift": phase_shift, "tk_unique": self.t_k, "KDinv":  KD_t, "Nhsiao": len(self.hsiao["t"]), "nint": np.cumsum(f_indices), "thsiao": self.hsiao["t"], "fhsiao": f_hsiao})

		#Initial values for model parameters
		stan_init = [{"Ds": mu + np.random.normal(0, np.sqrt(sigma0**2 + (fit_mode=="opt")*0.1**2)), "AV": np.random.exponential(tauA), "theta": np.random.normal(0,1), "epsilon_tform": np.random.normal(0, 1, self.params["L_Sigma_epsilon"].shape[0])} for _ in range(n_chains)]
		par_list = ("lp__", "mu", "delM", "AV", "theta")
		#Tmax initialisation
		if fix_tmax is False:
			par_list = par_list + ("tmax", "dtmax")
			tmax_shift = phase_shift*(1+z)
			for ini in stan_init:
				ini.update({"dtmax": np.random.uniform(-tmax_shift, tmax_shift)})

		#Compile, if not done so already
		if self.stan_model is None or self.fixed_tmax != fix_tmax:
			self.compile_stan_model(fix_tmax)
			self.fixed_tmax = fix_tmax

		#Check output_dir is valid, if given:
		if cmdstan_output_dir is not None and os.path.exists(cmdstan_output_dir):
			print("cmdstan outputs being written to: " + cmdstan_output_dir)
			output_dir = cmdstan_output_dir
		else:
			if cmdstan_output_dir is not None:
				print("WARNING: The provided cmdstan_output_dir ({}) is not real - run will be logged in tempfiles!".format(cmdstan_output_dir))
			output_dir = None

		with tempfile.TemporaryDirectory() as init_dir:
			#Write initialisation to temporary files
			print("initial values written to: " + init_dir)
			init_filenames = []
			for i in range(n_chains):
				init_filenames.append("{}/{}-init-{}.json".format(init_dir, self.stan_model.name, i))
				cmdstanpy.utils.write_stan_json(init_filenames[-1], stan_init[i])

			#Run Stan
			start = time.time()
			if fit_mode == "hmc":
				fit = self.stan_model.sample(data=stan_data, iter_sampling=n_sampling, iter_warmup=n_warmup, chains=n_chains, parallel_chains=n_workers, inits=init_filenames, refresh=1+9*(show_console), show_progress=show_progress, show_console=show_console, output_dir=output_dir)
			elif fit_mode == "opt":
				fit = self.stan_model.optimize(data=stan_data, inits=init_filenames[0], output_dir=output_dir, show_console=show_console)
			end = time.time()
			print("elapsed time: {:.2f} seconds".format(end - start))

		if fit_mode == "hmc":
			#Print diagnostics
			print(fit.diagnose())
			#Print summary table
			summary = fit.summary(sig_figs=6).round(2)
			summary = summary.loc[[*par_list]].append(summary.loc["epsilon_free[1]":"epsilon_free[{}]".format(self.params["L_Sigma_epsilon"].shape[0])])
			print("\n" + summary.to_string(index_names=False) + "\n")

			#Extract stacked post-warmup chains
			model_params = fit.stan_variables()
			sampler_params = fit.method_variables()
			samples = {par: model_params[par] for par in par_list + ("epsilon_free", "epsilon") if par != "lp__"}
			samples["lp__"] = sampler_params["lp__"]
		elif fit_mode == "opt":
			#Print summary table
			summary = fit.optimized_params_pd.round(2)
			summary = summary[[*par_list, *["epsilon_free[{:d}]".format(1+i) for i in range(self.params["L_Sigma_epsilon"].shape[0])]]].loc[0]
			print("\n" + summary.to_string() + "\n")

			#Extract result
			model_params = fit.stan_variables()
			samples = {par: model_params[par] for par in par_list + ("epsilon_free", "epsilon") if par != "lp__"}
			samples["lp__"] = fit.optimized_params_dict["lp__"]

		return fit, samples, summary

	def sample_AV(self, n=None, logM=None):
		"""
		Draw host extinction value(s) from the prior.

		Samples one or more host extinction values from the
		exponential distribution implied by the loaded model's
		tau_A value.

		Parameters
		----------
		n : int or None, optional
			Number of values to draw. Default is one.
		logM : float, :py:class:`numpy.array`, or None, optional
			log_10 of host galaxy stellar mass (only needed when using a
			model with a host galaxy mass split).

		Returns
		-------
		AV : float or :py:class:`numpy.array`
			AV value(s) drawn from the prior
		"""
		if self.M_split is not False:
			if logM is None:
				raise ValueError("This model contains a mass split - host galaxy mass must be provided!")
			else:
				tauA = self.params["tauA_HM"]*(logM >= self.M_split) + self.params["tauA_LM"]*(logM < self.M_split)
				return np.random.exponential(tauA, n)
		else:
			return np.random.exponential(self.params["tauA"], n)

	def sample_theta(self, n=None):
		"""
		Draw theta value(s) from the prior.

		Samples one or more values of the W1 coefficient, theta,
		from the assumed N(0,1) prior.

		The distribution here may not be reflective of the
		theta distribution of the training set.

		Parameters
		----------
		n : int or None, optional
			Number of values to draw. Default is one.

		Returns
		-------
		theta : float or :py:class:`numpy.array`
			theta value(s) drawn from the prior
		"""
		return np.random.normal(0, 1, n)

	def sample_del_M(self, n=None, logM=None):
		"""
		Draw gray offset value(s) from the prior.

		Samples one or more values of the delta M_s gray offset
		from the Gaussian implied by the loaded model's
		sigma_0 value.

		Parameters
		----------
		n : int or None, optional
			Number of values to draw. Default is one.
		logM : float, :py:class:`numpy.array`, or None, optional
			log_10 of host galaxy stellar mass (only needed when using a
			model with a host galaxy mass split).

		Returns
		-------
		del_M : float or :py:class:`numpy.array`
			delta M_s value(s) drawn from the prior
		"""
		if self.M_split is not False:
			if logM is None:
				raise ValueError("This model contains a mass split - host galaxy mass must be provided!")
			else:
				sigma0 = self.params["sigma0_HM"]*(logM >= self.M_split) + self.params["sigma0_LM"]*(logM < self.M_split)
				return np.random.normal(0, sigma0, n)
		else:
			return np.random.normal(0, self.params["sigma0"], n)

	def sample_epsilon(self, n=None, inc_zeros=True):
		"""
		Draw residual perturbation realisations(s) from the prior.

		Samples one or more epsilon (residual SED perturbation) matrices
		from the multivariate Gaussian implied by the loaded model's
		Sigma_epsilon matrix.

		Parameters
		----------
		n : int or None, optional
			Number of values to draw. Default is one.
		inc_zeros : bool, optional
			Toggles whether to fill in elements which are pinned to
			zero. Defaults to True - this should be left alone
			unless you have a specific reason for turning this off.

		Returns
		-------
		epsilons : :py:class:`numpy.array`
			epsilon_s value(s) drawn from the prior. These are reshaped
			so that each indivdual realisation matches the W0 and W1
			matrices in shape (number of rows = number of wavelength knots,
			number of cols = number of time knots). If `n != None`, first
			dimension of returned array will index realisations.
		"""
		mean = np.zeros((self.params["Sigma_epsilon"].shape[0],))
		e = np.random.multivariate_normal(mean, self.params["Sigma_epsilon"], n)
		eps = e.reshape((-1, len(self.l_k)-2, len(self.t_k)), order="F")
		if inc_zeros is True:
			eps = np.insert(eps, [0, eps.shape[-2]], np.zeros((eps.shape[-1],)), axis=-2)

		return np.squeeze(eps)

	def simulate_light_curve(self, flt, t=None, z=0, mu=0, ebv_mw=0, logM=None, del_M=None, AV=None, RV=None, theta=None, epsilon=None, mag=False, no_warn=False):
		"""
		Simulate a light curve from the BayeSN model.

		Generates a light curve in a single passband at times, `t`,
		given a set of latent parameters.

		Parameters
		----------
		flt : str
			String denoting which filter to simulate in. Must be
			a passband recognised by BayeSN.
		t : :py:class:`numpy.array` or None, optional
			Vector of rest frame phases to simulate at. If not provided,
			a daily cadence within the model range will be used.
		z : float, optional
			Heliocentric redshift to simulate light curve at. Defaults to
			0 (i.e. simulating in rest frame).
		mu : float, optional
			Distance modulus. Defaults to zero (simulating absolute mag
			light curve)
		ebv_mw : float, optional
			Milky Way reddening to apply. Defaults to zero.
		logM : float or None, optional
			log_10 of host galaxy stellar mass (only needed when using a model with
			a host galaxy mass split.
		del_M : float or None, optional
			Gray offset value. If None, generates a random value from
			the prior.
		AV : float or None, optional
			Host galaxy dust extinction. If None, generates a random
			value from the prior.
		theta : float or None, optional
			Coefficient for the W1 SED component. If None, generates from
			the N(0,1) prior.
		epsilon : :py:class:`numpy.array` or None, optional
			Intrinsic SED perturbation matrix. Should be same shape as
			W0 and W1 (number of rows = number of wavelength knots,
			number of cols = number of time knots). If None, generates a
			random realisation from the prior.
		mag : bool, optional
			If True, the returned light curve will be in magnitudes. Else,
			fluxes are given (defaults to False).
		no_warn : bool, optional
			If True, warnings about phase extrapolation will not be issued.
			Defaults to False.

		Returns
		-------
		t : :py:class:`numpy.array`
			Phases at which light curve has been simulated
		y : :py:class:`numpy.array`
			If `mag=False`, gives fluxes at `t`. Else, gives magnitudes at `t`.
		"""
		if flt not in self.passbands.keys():
			raise ValueError("Unrecognised passband: " + flt)
		#Interpolate filter transmission
		l_r, tp = utils.interpolate_passband(self.hsiao["l"], self.passbands[flt][0]["wave"], self.passbands[flt][0]["norm_throughput"], z)
		if np.any(l_r < min(self.l_k)) or np.any(l_r > max(self.l_k)):
			raise ValueError(flt + " filter extends beyond model coverage.")

		#Define daily sampled phase vector, if alternative not provided
		if t is None:
			t = np.linspace(min(self.t_k), max(self.t_k), int(max(self.t_k) - min(self.t_k)) + 1)
		elif no_warn is False and (max(t) > max(self.t_k) or min(t) < min(self.t_k)):
			print("WARNING: Points simulated at phases outside of [{:.1f}, {:.1f}] will be extrapolated - treat with caution.".format(min(self.t_k), max(self.t_k)))

		if self.M_split is not False and logM is None:
			raise ValueError("This model contains a mass split - host galaxy mass must be provided!")

		#Draw random values of model params if not provided
		if del_M is None:
			del_M = self.sample_del_M(logM=logM)
		if AV is None:
			AV = self.sample_AV(logM=logM)
		if RV is None:
			RV = self.params["RV"]
		if theta is None:
			theta = self.sample_theta()
		if epsilon is None:
			epsilon = self.sample_epsilon()

		#Spline interpolation matrices
		KD_t = spline_utils.invKD_irr(self.t_k)
		KD_l = spline_utils.invKD_irr(self.l_k)
		J_t = spline_utils.spline_coeffs_irr(t, self.t_k, KD_t).T
		J_l = spline_utils.spline_coeffs_irr(l_r, self.l_k, KD_l)

		#Observer frame wavelengths in this filter
		l_o = l_r*(1+z)
		#Milky Way reddening (in flux units)
		R_MW = utils.extinction_mw(l_o, ebv_mw)
		#Passband integral
		H = tp*l_o*self.dl*R_MW

		#Host galaxy reddening vector (in magnitudes)
		if self.M_split is not False:
			if logM >= self.M_split:
				R_host = extinction.fitzpatrick99(l_r, AV, self.params["RV_HM"])
				M0 = self.params["M0_HM"]
			else:
				R_host = extinction.fitzpatrick99(l_r, AV, self.params["RV_LM"])
				M0 = self.params["M0_LM"]
		else:
			R_host = extinction.fitzpatrick99(l_r, AV, RV)
			M0 = self.params["M0"]

		#Intrinsic SED spline knots
		W = self.params["W0"] + theta*self.params["W1"] + epsilon
		#Interpolate
		JWJ = np.linalg.multi_dot([J_l, W, J_t])

		#Hsiao template
		S0 = np.zeros(JWJ.shape)
		for i in range(len(t)):
			S0[:,i] = utils.interpolate_hsiao(t[i], l_r, self.hsiao["t"], self.hsiao["l"], self.hsiao["f"])

		gamma = np.log(10)/2.5

		#Host-extinguished rest frame SED
		S_tilde = S0*np.exp(-gamma*(JWJ + R_host[:,None]))
		#Light curve
		if mag is False:
			y = np.exp(gamma*(self.zpt - M0 - mu - del_M))*np.dot(H, S_tilde)
		else:
			y = mu + M0 + del_M - 2.5*np.log10(np.dot(H, S_tilde))

		return t, y

	def simulate_spectrum(self, t, l_o=None, z=0, mu=0, ebv_mw=0, logM=None, del_M=None, AV=None, theta=None, epsilon=None, log=False):
		"""
		Simulate a slice of the (log) SED from the BayeSN model.

		Generates a the SED or log SED a single time, `t`, for a range of
		wavelengths given a set of latent parameters.

		Parameters
		----------
		t : float
			Rest frame phase to simulate at.
		l_o : :py:class:`numpy.array` or None, optional
			Vector of observer frame wavelengths to simulate at. If not provided,
			10A spacing within the model range will be used.
		z : float, optional
			Heliocentric redshift to simulate SED at. Defaults to
			0 (i.e. simulating in rest frame).
		mu : float, optional
			Distance modulus. Defaults to zero (simulating absolute SED)
		ebv_mw : float, optional
			Milky Way reddening to apply. Defaults to zero.
		logM : float or None, optional
			log_10 of host galaxy stellar mass (only needed when using a model with
			a host galaxy mass split.
		del_M : float or None, optional
			Gray offset value. If None, generates a random value from
			the prior.
		AV : float or None, optional
			Host galaxy dust extinction. If None, generates a random
			value from the prior.
		theta : float or None, optional
			Coefficient for the W1 SED component. If None, generates from
			the N(0,1) prior.
		epsilon : :py:class:`numpy.array` or None, optional
			Intrinsic SED perturbation matrix. Should be same shape as
			W0 and W1 (number of rows = number of wavelength knots,
			number of cols = number of time knots). If None, generates a
			random realisation from the prior.
		log : bool, optional
			If True, the log SED is returned. Else, SED is returned in
			flux units (defaults to False).

		Returns
		-------
		l_o : :py:class:`numpy.array`
			Observer frame wavelengths at which SED has been simulated
		y : :py:class:`numpy.array`
			If `log=False`, gives SED at `l_o`. Else, gives log SED at `l_o`.

		Notes
		-----
		If `log=False`, returned quantity, `F_obs`, is the observable flux
		density (M20, eq. 5). This can be integrated through a passband per
		M20 eq. 4 to yield an apparent magnitude, If `log=True`, returned
		quantity is `2.5*log10` of the host-extinguished, Milky-Way reddened
		apparent SED. This equals `2.5*log10[F_obs*(1+z)]`. With `z=0`,
		`mu=0`, `ebv_mw=0`, this is `2.5*log10(S_s)`, where S_s is the
		host-dust-extinguished SED (c.f. M20, eq. 12) in the rest frame.
		"""

		#Generic wavelength vector if one isn't provided
		if l_o is None:
			l_r = np.linspace(min(self.l_k), max(self.l_k), int((max(self.l_k) - min(self.l_k))/self.dl) + self.dl)
			l_o = l_r*(1+z)
		else:
			l_r = l_o/(1+z)
			if max(l_r) > max(self.l_k) or min(l_r) < min(self.l_k):
				raise ValueError("Rest wavelength extends beyond model coverage.")

		if self.M_split is not False and logM is None:
			raise ValueError("This model contains a mass split - host galaxy mass must be provided!")

		#Draw random values of model params if not provided
		if del_M is None:
			del_M = self.sample_del_M(logM=logM)
		if AV is None:
			AV = self.sample_AV(logM=logM)
		if theta is None:
			theta = self.sample_theta()
		if epsilon is None:
			epsilon = self.sample_epsilon()

		#Spline interpolation matrices
		KD_t = spline_utils.invKD_irr(self.t_k)
		KD_l = spline_utils.invKD_irr(self.l_k)
		J_t = spline_utils.spline_coeffs_irr([t], self.t_k, KD_t).T
		J_l = spline_utils.spline_coeffs_irr(l_r, self.l_k, KD_l)

		#Milky Way reddening (in flux units)
		R_MW = utils.extinction_mw(l_o, ebv_mw)

		#Host galaxy reddening vector (in magnitudes)
		if self.M_split is not False:
			if logM >= self.M_split:
				R_host = extinction.fitzpatrick99(l_r, AV, self.params["RV_HM"])
				M0 = self.params["M0_HM"]
			else:
				R_host = extinction.fitzpatrick99(l_r, AV, self.params["RV_LM"])
				M0 = self.params["M0_LM"]
		else:
			R_host = extinction.fitzpatrick99(l_r, AV, self.params["RV"])
			M0 = self.params["M0"]

		#Intrinsic SED spline knots
		W = self.params["W0"] + theta*self.params["W1"] + epsilon
		#Interpolate
		JWJ = np.linalg.multi_dot([J_l, W, J_t]).squeeze()

		#Hsiao template
		S0 = utils.interpolate_hsiao(t, l_r, self.hsiao["t"], self.hsiao["l"], self.hsiao["f"])

		gamma = np.log(10)/2.5

		if log is False:
			#Host-extinguished rest frame SED
			S_tilde = S0*np.exp(-gamma*(JWJ + R_host))
			#SED slice in flux units
			y = np.exp(-gamma*(M0 + mu + del_M))*R_MW*S_tilde/(1+z)
		else:
			#log of SED slice
			y = 2.5*np.log10(R_MW*S0) - mu - M0 - del_M - JWJ - R_host

		return l_o, y

	def pppvalue(self, lc, samples, discrep="chi2", z=None, ebv_mw=None, logM=None, tmax=None, snr_cut=3):
		"""
		Estimates posterior predictive p-value

		Estimates the posterior predictive p-value (Meng94, Gelman+96)
		for a light curve fit, using chi-square as the default
		discrepancy statistic.

		Parameters
		----------
		lc : :py:class:`astropy.table.Table`
			Astropy table containing light curve information. This must contain
			`phase` (and/or `mjd`), `fluxcal`, `fluxcalerr` and `flt` columns.
			It is recommended that this table be initialised from an SNANA format
			light curve file using :py:func:`BayeSNmodel.io.read_snana_lcfile`, since
			this should retrieve correctly formatted columns, and well populated
			metadata.
		samples : dict
			Dictionary of MCMC samples as returned by :py:func:`fit_supernova'.
		discrep : string, optional
			Discrepancy statistic to use when computing p-value. Defaults to
			chi-square (``'chi2'``). Other options are mean (``'resid_mean'``)
			and stddev (``'resid_std'``) of residuals, or mean squared error (``'mse'``).
		z : float or None, optional
			Heliocentric redshift used for redshifting SED. Doesn't need to
			be specified explicitly if contained in `lc.meta` as ``'REDSHIFT_HELIO'``.
		ebv_mw : float or None, optional
			Milky Way reddening. Needn't be given if included in `lc.meta` as
			``'MWEBV'``.
		logM : float or None, optional
			log_10 of host galaxy stellar mass (only needed when using a model with
			a host galaxy mass split).
		tmax : float or None, optional
			Time of B-band maximum in MJD. If tmax was fixed in fits, this will be
			used as the time of maximum. If tmax was fit for, it will be assumed that
			this is the initial guess which was used. Needn't be given if included in `lc.meta` as ``'PEAKMJD'`` or	``'SEARCH_PEAKMJD'``.
		snr_cut : float, optional
			Minimum S/N cut applied to light curves before fitting. Default=3.

		Returns
		-------
		pppv : float
			Posterior predictive p-value
		D_y : :py:class:`numpy.array`
			Vector of realised discrepancies
		D_y_rep : py:class:`numpy.array`
			Vector of predictive discrepancies
		"""

		#Check requested discrepancy stat is valid
		if discrep not in ["chi2", "resid_mean", "resid_std", "mse"]:
			raise ValueError("Requested discrepancy statistic is not recognised. Please chose from chi2, resid_mean or resid_std.")

		#Check for heliocentric redshift
		if z is None:
			try:
				z = lc.meta["REDSHIFT_HELIO"]
			except KeyError as e:
				raise Exception("Redshift not provided or found in light curve metadata.") from e
		#Check for MW reddening
		if ebv_mw is None:
			try:
				ebv_mw = lc.meta["MWEBV"]
			except KeyError as e:
				raise Exception("Milky Way reddening not provided or found in light curve metadata.") from e
		#Check for hostmass (if necessary)
		if self.M_split is not False and logM is None:
			try:
				logM = lc.meta["HOSTGAL_LOGMASS"]
			except KeyError as e:
				raise Exception("Host galaxy mass not provided or found in light curve metadata") from e
		#Check for tmax
		if tmax is None:
			try:
				tmax = lc.meta["SEARCH_PEAKMJD"]
			except KeyError:
				try:
					tmax = lc.meta["PEAKMJD"]
				except KeyError as e:
					raise Exception("Time of maximum not provided or found in light curve metadata.") from e

		#See if tmax is free
		if "tmax" in samples.keys():
			fix_tmax = False
		else:
			fix_tmax = True

		#Computes phase
		lc["phase"] = (lc["mjd"] - tmax)/(1 + z)

		#Phase and S/N cut
		mask = (lc["phase"] > min(self.t_k))*(lc["phase"] < max(self.t_k))*(lc["fluxcal"]/lc["fluxcalerr"] > snr_cut)
		lc = lc[mask]

		D_y_rep = np.empty((0,))
		D_y = np.empty((0,))

		#Cycle through MCMC chains
		for i in range(len(samples["mu"])):

			d_rep_i = np.empty((0,))
			d_i = np.empty((0,))

			if fix_tmax is True:
				phase = lc["phase"]
			else:
				phase = (lc["mjd"] - samples["tmax"][i])/(1 + z)

			#Cycle through filters
			for f, flt in enumerate(unique(lc, keys="flt")["flt"]):
				flt_mask = lc["flt"] == flt

				#Simulate model light curve for this MCMC sample
				_, y_model_i_f = self.simulate_light_curve(flt, t=phase[flt_mask], z=z, mu=samples["mu"][i], ebv_mw=ebv_mw, logM=logM, del_M=samples["delM"][i], AV=samples["AV"][i], theta=samples["theta"][i], epsilon=samples["epsilon"][i], mag=False, no_warn=True)
				#Generate repeated observation
				y_rep_i_f = np.random.normal(y_model_i_f, lc["fluxcalerr"][flt_mask])

				if discrep == "chi2":
					d_rep_i = np.append(d_rep_i, (y_rep_i_f - y_model_i_f)**2/lc["fluxcalerr"][flt_mask]**2)
					d_i = np.append(d_i, (lc["fluxcal"][flt_mask] - y_model_i_f)**2/lc["fluxcalerr"][flt_mask]**2)
				elif discrep in ["resid_mean", "resid_std", "mse"]:
					d_rep_i = np.append(d_rep_i, y_rep_i_f - y_model_i_f)
					d_i = np.append(d_i, lc["fluxcal"][flt_mask] - y_model_i_f)

			if discrep == "chi2":
				D_y_rep = np.append(D_y_rep, np.sum(d_rep_i))
				D_y = np.append(D_y, np.sum(d_i))
			elif discrep == "resid_mean":
				D_y_rep = np.append(D_y_rep, np.mean(d_rep_i))
				D_y = np.append(D_y, np.mean(d_i))
			elif discrep == "resid_std":
				D_y_rep = np.append(D_y_rep, np.std(d_rep_i))
				D_y = np.append(D_y, np.std(d_i))
			elif discrep == "mse":
				D_y_rep = np.append(D_y_rep, np.mean(d_rep_i**2))
				D_y = np.append(D_y, np.mean(d_i**2))

		pppv = np.sum(D_y_rep >= D_y)/len(D_y)

		return pppv, D_y, D_y_rep
