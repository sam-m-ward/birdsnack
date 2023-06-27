"""
SBC Module==Simulation-Based Calibration

Module contains SBC_CLASS for simulating data with known hyperparameters then testing posterior recovery

Contains
--------------------
SBC_CLASS class
	inputs : choices, edit_dict={}

	Methods are:
		get_simulations_folder()
		simulate_truths()
		get_truths(index=None)
		get_leffs_df(TRUTHS_DICT=None)
		fit_truths(TRUTHS_DICT=None)
		get_fits(TRUTHS_DICT=None)

SIMULATOR class
	inputs : choices

	Methods are:
		trim_Sigma(Sigma)
		get_dust_pop()
		get_dust_ind()
		get_int_pop()
		get_int_ind()
		get_latent_parameters()
		get_observations()

Functions include:
	get_Lint_prior_samples(Nc=6,savefolder='products/Lint_sims/',generator_file='generator.stan',n_sampling=10000,n_chains=1,n_warmup=10)
	def get_DF_M_from_truths(truths,pblist,errstr,tref)

--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""
from snpy import fset
import matplotlib#Reset matplotlib params
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
import numpy as np
import copy, extinction, os, pickle, sys
from contextlib import suppress
from glob import glob
import pandas as pd
from scipy.stats import truncnorm, gamma

class SBC_CLASS:
	def get_simulations_folder(self):
		"""
		Get Simulations Folder

		Creates folder name based on .yaml choices

		End Product(s)
		----------
		simfolder,savepath
		"""
		#Get Miscellaneous function(s)
		sys.path.append(self.path_to_birdsnack_rootpath+'model_files/')
		from miscellaneous import ensure_folders_to_file_exist

		#NSNe
		folder     = f"S{self.S}"
		#Simulator
		folder    += f"_simulator{self.simulator}"
		if self.simulator=='BayeSN':#BayeSN Parameters
			bypars = self.bayesn_parameters
			folder  += f"_bymodel{bypars['bymodel']}"
			folder  += f"_theta{bypars['thetamode']}"
		#Intrinsic Variations
		folder    += f"_epsilon{self.epsmode}"
		folder    += f"_FPC0m{self.FPC0m}"
		#Dust Hyperparameters
		folder    += f"_tauA{self.tauA}_muRV{self.muRV}_sigRV{self.sigRV}"
		#AV Population Distribution
		folder    += f"_AVsimdist{self.AVsimdist}"
		if self.AVsimdist!='Exp':#AVsimdist (hyper)parameters
			if self.AVsimdist in ['Gamma']:
				folder  += f"_nuA{str(self.AV_parameters['nuA'])}"
		#Intrinsic Hyps
		folder    += f"_PredefinedIntrinsicHyps{self.PredefinedIntrinsicHyps}"
		if self.PredefinedIntrinsicHyps:
			folder += f"_loadedfrom{self.pre_defined_hyps['load_file']}"
		#Final Addition to make it a folder
		folder += '/'

		self.simfolder   = folder
		self.simsavepath = self.simspath + self.simfolder
		#Make sure path to where these truth sims will be saved exists
		ensure_folders_to_file_exist(self.simsavepath)
		#Also, make sure path to folder where random draws of correlation matrices are saved exists
		ensure_folders_to_file_exist(f"{self.productpath}Lint_sims/")

	def __init__(self, choices, edit_dict={}):
		##Load up sbc.yaml
		#Ignore global keys and update with local key,value pairs only
		for glob_key in choices:
			for key,value in choices[glob_key].items():
				setattr(self, key, value)
		#Update key,value edit_dict pairs
		for glob_key in edit_dict:
			for key,value in edit_dict[glob_key].items():
				setattr(self, key, value)

		#Initialise SBC Paths
		self.rootpath     = self.path_to_rootpath
		self.plotpath     = f"{self.rootpath}plots/"
		self.productpath  = f"{self.rootpath}products/"
		self.BayeSNpath   = f"{self.rootpath}BAYESN_GITHUB/"
		self.simspath     = f"{self.productpath}sbc_sims/"

		#Get folder where truth simulations will be saved
		self.get_simulations_folder()

		#Initialise BayeSN items
		if self.simulator=='BayeSN':
			bypars = self.bayesn_parameters
			self.bymodel      = bypars['bymodel']
			self.thetamode    = bypars['thetamode']
			os.environ["PYSYN_CDBS"] = f"{self.BayeSNpath}pysynphot_downloads/grp/hst/cdbs/"
			sys.path.append(self.BayeSNpath+'bayesn-public/')
			from BayeSNmodel import bayesn_model
			b = bayesn_model.SEDmodel(model=self.bymodel)
			self.b = b

		self.dict = self.__dict__
		self.TRUTHS_DICT = None

	def simulate_truths(self):
		"""
		Simulate Truths

		Method to simulate SNe and store true (hyper)parameters
		Limiting factor on computation time is .get_latent_parameters() when using simulator='BayeSN'

		End Product(s)
		----------
		Sim_truths.pkl files
		"""
		for _ in range(self.Nsims):
			if not os.path.exists(self.simsavepath+f'Sim{_}_truths.pkl'):
				if _%10==0:	print (_,'/',self.Nsims)
				truths = SIMULATOR(self.dict)
				with open(self.simsavepath+f'Sim{_}_truths.pkl','wb') as f:
					pickle.dump(truths,f)
		print ('Done')
		print ('####'*10)


	def get_truths(self,index=None):
		"""
		Get Truths

		Used to load up set of simulations

		Parameters
		----------
		index : float (optional; default is None)
			if not None, load up Sim{index}
			otherwise load up all the files

		End Product(s)/Return
		----------
		TRUTHS_DICT : dict
			key,value are index,SIMULATOR class object
		"""
		files = glob(self.simsavepath+'*truths.pkl')
		TRUTHS_DICT = {}
		for _,file in enumerate(files):
			if index is None:
				with open(file,'rb') as f:
					truths = pickle.load(f)
				TRUTHS_DICT[_] = truths
			else:
				if file == sbc.simsavepath+f'Sim{index}_truths.pkl':
					with open(file,'rb') as f:
						truths = pickle.load(f)
					TRUTHS_DICT[index] = truths
					break
		self.TRUTHS_DICT = TRUTHS_DICT
		return self.TRUTHS_DICT

	def get_leffs_df(self, TRUTHS_DICT = None):
		"""
		Get Effective Wavelengths DataFrame

		Method to load up simulated SN datasets,
		then compute median of leffs for 1 dataset
		then store df of leffs across datasets

		Parameters
		----------
		TRUTHS_DICT : dict (optional; default=None)
			key,value are index,SIMULATOR class object

		End Product(s)/Returns
		----------
		lameff_df : df
			each column is passband
			each row is leff for that simulation
		"""
		#Get Truths
		if TRUTHS_DICT is None:
			if self.TRUTHS_DICT is None:
				TRUTHS_DICT = self.get_truths()
			else:
				TRUTHS_DICT = self.TRUTHS_DICT
		#Initialise pandas df
		df = pd.DataFrame(columns=[s[0] for s in self.flts]) ; ISIM = -1
		#Loop through simulated datasets
		for key,truths in TRUTHS_DICT.items():
			ISIM += 1
			if ISIM%10==0:	print (f'{ISIM} / {len(TRUTHS_DICT)}')
			#Get lameff entries, 1 for each SN
			l_effs_list = []
			for s in range(truths.S):
				eps_hat = (truths.mexts[s]-truths.mints[s])/truths.AVs[s]
				l_effs = []
				for il,eps in enumerate(eps_hat):
					lams    = np.arange(truths.l_eff_rest[il]-300,truths.l_eff_rest[il]+100+1,0.1)
					eps_eff = extinction.fitzpatrick99(lams,1,truths.RVs[s])
					l_effs.append(lams[np.argmin(np.abs(eps_eff-eps))])
				l_effs_list.append(np.asarray(l_effs))
			#With all SNe, get median, and append this to global row
			l_effs_list = np.asarray(l_effs_list)
			append_dict = dict(zip(df.columns,np.median(l_effs_list,axis=0)))
			df          = df.append(append_dict,ignore_index=True)
		#Assign attribute and return
		self.lameff_df = df
		return self.lameff_df

	def fit_truths(self, TRUTHS_DICT=None):
		"""
		Fit Truths Method

		Parameters
		----------
		TRUTHS_DICT : dict (optional; default=None)
			key,value are index,SIMULATOR class object

		End Product(s)
		----------
		FIT.pkl files of HBM fits to simulated data
		"""
		#Get Truths
		if TRUTHS_DICT is None:
			if self.TRUTHS_DICT is None:
				TRUTHS_DICT = self.get_truths()
			else:
				TRUTHS_DICT = self.TRUTHS_DICT

		#Initialise Bird-Snack Model
		sys.path.append(f"{self.path_to_birdsnack_rootpath}model_files/")
		from birdsnack_model import BIRDSNACK
		bs      = BIRDSNACK(configname=f"{self.birdsnack_yaml}.yaml")
		pblist  = bs.choices['preproc_parameters']['pblist']
		errstr  = bs.choices['preproc_parameters']['errstr']
		tref    = bs.choices['preproc_parameters']['tilist'][bs.choices['preproc_parameters']['tref_index']]

		#Loop through truths
		for ISIM,truths in TRUTHS_DICT.items():
			save_filename = f"{self.simsavepath}FIT{bs.choices['analysis_parameters']['HBM_savekey']}_Sim{ISIM}.pkl"
			if not os.path.exists(save_filename):
				print ("###"*10)
				print (f"Performing fit to ISIM={ISIM};")
				print (f"save_filename is {save_filename}")
				#Get DF_M for BirdSnack
				DF_M = get_DF_M_from_truths(truths,pblist,errstr,tref)
				bs.DF_M = DF_M
				#Fit HBM
				bs.fit_stan_model(save=False,Rhat_threshold=1.05)
				#Extract FIT, thin df to x3 dust hyperparameters
				FIT = bs.FIT
				FIT['df'] = FIT['df'][['mu_RV','sig_RV','tauA']]
				#Save FIT
				with open(save_filename,'wb') as f:
					pickle.dump(FIT,f)

	def get_fits(self,TRUTHS_DICT=None):
		"""
		Fit FITS Method

		Parameters
		----------
		TRUTHS_DICT : dict (optional; default=None)
			key,value are index,SIMULATOR class object

		End Product(s)
		----------
		FITS : dict
			{ISIM:FIT} of SBC fits
		"""
		#Get Truths
		if TRUTHS_DICT is None:
			if self.TRUTHS_DICT is None:
				TRUTHS_DICT = self.get_truths()
			else:
				TRUTHS_DICT = self.TRUTHS_DICT

		#Initialise Bird-Snack Model
		sys.path.append(f"{self.path_to_birdsnack_rootpath}model_files/")
		from birdsnack_model import BIRDSNACK
		bs = BIRDSNACK(configname=f"{self.birdsnack_yaml}.yaml")

		#Load FITS
		FITS = {}
		for ISIM,truths in TRUTHS_DICT.items():
			save_filename = f"{self.simsavepath}FIT{bs.choices['analysis_parameters']['HBM_savekey']}_Sim{ISIM}.pkl"
			with open(save_filename,'rb') as f:
				FIT = pickle.load(f)
			FITS[ISIM] = FIT
		self.bs   = bs
		self.FITS = FITS
		return self.FITS



class SIMULATOR():

	def trim_Sigma(self,Sigma):
		"""
		Trim Sigma

		Simple method to trim covariance matrix down to appropriate passband

		Returns
		----------
		Sigma : array
			trimmed version of Sigma==Covariance Matrix
		"""
		Sigma = np.delete(Sigma, self.kill_indices, 0)
		Sigma = np.delete(Sigma, self.kill_indices, 1)
		return Sigma

	def __init__(self,  choices,
						tauAmin=0.1,tauAmax=1,muRVmin=1,muRVmax=5,sigRVmin=0.1,sigRVmax=2,
						sigma_int_min=0,sigma_int_max=0.25,RVmin=1,RVmax=np.inf,
						mumin=30,mumax=60,Emin = 0.005, Emax = np.inf, Eloc = 0.02, Escal = 0.02,
						**kwargs
						):

		#Set key,value pairs in choices
		for key,value in choices.items():
			setattr(self, key, value)

		#Number of passbands
		self.Nm         = len(self.flts)
		#For use on perfect simulator only when adding on AV*xi(lambda,RV)
		self.l_eff_rest = np.array([fset[flt[0]].ave_wave for flt in self.flts])

		#Dust Population Hyperparameters Simulation Limits
		self.tauAmin  = tauAmin
		self.tauAmax  = tauAmax
		self.muRVmin  = muRVmin
		self.muRVmax  = muRVmax
		self.sigRVmin = sigRVmin
		self.sigRVmax = sigRVmax

		#Intrinsic SED Hyperparameter Simulation Limits
		self.sigma_int_min = sigma_int_min
		self.sigma_int_max = sigma_int_max

		#Individual Dust Parameter Simulation Limits
		self.RVmin   = RVmin
		self.RVmax   = RVmax

		#Other Simulation Limits
		self.mumin = mumin
		self.mumax = mumax

		self.Emin  = Emin
		self.Emax  = Emax
		self.Eloc  = Eloc
		self.Escal = Escal

		#BayeSN parameters
		if self.simulator=='BayeSN':
			self.Model_Passband_Dict =  { 'M20':dict(zip(['B','V','r','i','Y','J','H'],list(np.arange(7)))),
										  'W22':dict(zip(['B','g','V','r','i','z','Y','J','H'],list(np.arange(9))))
										}
			self.keep_indices = [self.Model_Passband_Dict[self.bymodel][flt[0]] for flt in self.flts]
			self.kill_indices = [i for i in range(len(self.Model_Passband_Dict[self.bymodel])) if i not in self.keep_indices]

		#Update with additional/overwritten kwargs
		for key,value in kwargs.items():
			setattr(self, key, value)
		#Get Dust Hyperparameters (if not already defined)
		self.get_dust_pop()
		#Simulate each SN's dust parameters (AVs,RVs)
		self.get_dust_ind()
		#Get Intrinsic Hyperparameters (if not already defined)
		self.get_int_pop()
		#Simulate each SN's intrinsic parameters
		self.get_int_ind()
		#Get Latent Params Before Observation
		self.get_latent_parameters()
		#Get Simulated Data using fake measurement errors
		self.get_observations()
		#We save the entire SIMULATOR class, but don't save BayeSN model, so instead replace this with a string
		with suppress(AttributeError):
			self.b = 'not_saved_due_to_space_complexity'

	def get_dust_pop(self):
		"""
		Get Dust Population Hyperparameters

		Method to get random draw of dust hyperparameters if they are not already defined

		End Products
		----------
		tauA,muRV,sigRV
		"""
		if self.tauA is None:
			self.tauA  = np.random.uniform(self.tauAmin,self.tauAmax)
		if self.muRV is None:
			self.muRV  = np.random.uniform(self.muRVmin,self.muRVmax)
		if self.sigRV is None:
			self.sigRV = np.random.uniform(self.sigRVmin,self.sigRVmax)

	def get_dust_ind(self):
		"""
		Get Dust Individual

		Method to get draws of dust parameters from population distributions and hyperparameters

		End Products
		----------
		AVs,RVs,dustlaws
		"""
		AVs,RVs = [],[]
		for s in range(self.S):
			#Get AVs
			if self.AVsimdist=='Exp': AV = np.random.exponential(self.tauA)
			elif self.AVsimdist=='Gamma': AV = gamma.rvs(a=self.nuA,loc=0,scale=self.tauA)
			AVs.append(AV)
			#Get RVs
			RV = truncnorm.rvs( (self.RVmin-self.muRV)/self.sigRV, (self.RVmax-self.muRV)/self.sigRV, loc=self.muRV,scale=self.sigRV)
			RVs.append(RV)
		#Assign Attributes
		self.AVs      = np.asarray(AVs)
		self.RVs      = np.asarray(RVs)
		self.dustlaws = [extinction.fitzpatrick99(self.l_eff_rest,1,RV) for RV in RVs]#Only used if simulator=='perfect'

	def get_int_pop(self):
		"""
		Get Intrinsic Population

		Method to get intrinsic population hyperpameters if not already defined

		End Product(s)
		----------
		L_mint, sigma_mints, FPC0m
		"""
		if self.PredefinedIntrinsicHyps:
			with open(f"{self.path_to_birdsnack_rootpath}products/stan_fits/FITS/FIT{self.pre_defined_hyps['load_file']}.pkl",'rb') as f:
				FIT = pickle.load(f)
			row    = FIT['df'].median(axis=0)
			FPC0m  = row[[col for col in row.keys() if 'FPC0m' in col and 'simp' not in col and 'mbar' not in col]].values
			L_mint = row[[col for col in row.keys() if 'L_mint' in col and 'eta' not in col and 'Cens' not in col]].values.reshape(6,6).T
		else:
			#Get Epsilon Hyperparameters
			if self.epsmode=='random':
				L_mint_eta    = get_Lint_prior_samples(Nc=self.Nm,savefolder=f"{self.productpath}Lint_sims/",generator_file=f"{self.rootpath}model_files/stan_files/generator.stan")
				L_mint_eta    = L_mint_eta[np.random.choice(np.arange(len(L_mint_eta)))]
				sigma_mints   = np.random.uniform(self.sigma_int_min,self.sigma_int_max,self.Nm)
				L_mint        = np.diag(sigma_mints)@L_mint_eta@np.diag(sigma_mints)
			elif self.epsmode=='bymodel':
				Lint = np.loadtxt(f"{self.BayeSNpath}bayesn-model-files/BAYESN.{self.bymodel}/L_Sigma_epsilon.txt")
				Sigma_mint = Lint@Lint.T
				#indices here trim to peak time, method removes e.g. g,z,Y filters
				if self.bymodel == 'W22':
					Sigma_mint = self.trim_Sigma(Sigma_mint[9:18,9:18])
				if self.bymodel == 'M20':
					Sigma_mint = self.trim_Sigma(Sigma_mint[7:14,7:14])
				sigma_mints = np.diag(Sigma_mint)**0.5
				L_mint      = np.linalg.cholesky(Sigma_mint)
			elif self.epsmode==0:
				sigma_mints = np.zeros(self.Nm)
				L_mint      = np.zeros((self.Nm,self.Nm))
			#Get FPC0m hyperparameter
			if self.FPC0m=='random':
				FPC0m = np.random.normal(0,1,self.Nm)
			elif self.FPC0m=='bymodel':
				FPC0m = np.zeros(self.Nm)
				for i,flt in enumerate(self.flts):
					t, m     = self.b.simulate_light_curve(flt, AV=0, t=[0], del_M=0, mag=True, theta=0, epsilon=0)
					FPC0m[i] = m[0]

		#Assign hyperparameters
		self.sigma_mints = sigma_mints
		self.L_mint      = L_mint
		self.FPC0m = FPC0m


	def get_int_ind(self):
		"""
		Get Intrinsic Individual

		Method to get draws of intrinsic parameters from population distributions and population hyperparameters

		End Product(s)
		----------
		dMis, mus, epsilons, thetas
		"""
		#Get distances and draws from covariance matrix
		self.mus  = np.random.uniform(self.mumin,self.mumax,self.S)
		self.dMis = [self.L_mint@np.random.normal(0,1,self.Nm) for s in range(self.S)]
		#Get epsilon draws and thetas if BayeSN is used as simulator
		if self.simulator == 'BayeSN':
			epsilons = [] ; thetas = []
			for s in range(self.S):
				####################
				#Epsilons
				if self.epsmode=='random':
					epsilon = self.b.sample_epsilon()
					counter = 0 ; keep_indices = np.asarray(copy.deepcopy(self.keep_indices))+1#Add on 1 because first wavelength knot is tied at zero
					for i in range(epsilon.shape[0]):
						if i in keep_indices:
							epsilon[i,1] = self.dMis[s][counter]#1 is time knot of 0days (-10,0,10,20,30,40)
							counter += 1
				elif self.epsmode=='bymodel':
					epsilon = self.b.sample_epsilon()
				elif self.epsmode==0:
					epsilon = 0
				epsilons.append(epsilon)
				####################
				#Thetas
				if self.thetamode=='random':
					theta = np.random.normal(0,1)
					while theta<-1.5 or theta>2:
						theta = np.random.normal(0,1)
				elif self.thetamode == 0:
					theta = 0
				thetas.append(theta)
				####################
			self.epsilons = epsilons
			self.thetas   = thetas

	def get_latent_parameters(self):
		"""
		Get Latent Parameters

		Combines extrinsic, intrinsic and achromatic parameters together to get latent apparent magnitude parameters

		End Product(s)
		----------
		mints,mexts,dm15Bs
		"""
		mints = [] ; mexts = [] ; dm15Bs = []
		#Perfect simulator is simple linear combination
		if self.simulator == 'perfect':
			for s in range(self.S):
				mint = self.mus[s] + self.FPC0m + self.dMis[s]
				mints.append(mint)
				mext = mint + self.AVs[s]*self.dustlaws[s]
				mexts.append(mext)
				dm15Bs.append(1.05)
		#BayeSN simulation
		elif self.simulator == 'BayeSN':
			for s in range(self.S):
				mi = [] ; me = []
				for flt in self.flts:
					t, mext = self.b.simulate_light_curve(flt,RV=self.RVs[s],AV=self.AVs[s], t=[0], del_M=0, mag=True, theta=self.thetas[s], epsilon=self.epsilons[s], mu=self.mus[s])
					me.append(mext[0])
					t, mint = self.b.simulate_light_curve(flt,RV=self.RVs[s],AV=0,           t=[0], del_M=0, mag=True, theta=self.thetas[s], epsilon=self.epsilons[s], mu=self.mus[s])
					mi.append(mint[0])
					if flt=='B_CSP':
						t, mext15 = self.b.simulate_light_curve(flt,RV=self.RVs[s],AV=self.AVs[s], t=[15], del_M=0, mag=True, theta=self.thetas[s], epsilon=self.epsilons[s], mu=self.mus[s])
						dm15Bs.append(mext15[0]-mext[0])
				mexts.append(np.asarray(me))
				mints.append(np.asarray(mi))
		#Assign attributes
		self.mints  = mints
		self.mexts  = mexts
		self.dm15Bs = dm15Bs

	def get_observations(self):
		"""
		Get Observations

		Takes in latent parameters, and gets observed 'data' using measurement errors

		End Product(s)
		----------
		errors,mobs,dm15Bobs
		"""
		errors = [] ; mobs = [] ; dm15Bobs = []
		Emin,Emax,Eloc,Escal = self.Emin,self.Emax,self.Eloc,self.Escal
		for s in range(self.S):
			#Get random draw of error dispersion from trunctaed normal
			error = truncnorm.rvs((Emin-Eloc)/Escal,(Emax-Eloc)/Escal,loc=Eloc,scale=Escal,size=self.Nm)
			#Get Gaussian draw from N(0,error)
			mo    = np.random.normal(self.mexts[s],error)
			#Append dispersion term
			errors.append(error)
			#Append observation
			mobs.append(mo)

			#Get another error dispersion for m15B and repeat as above
			m15Bobserr  = truncnorm.rvs((Emin-Eloc)/Escal,(Emax-Eloc)/Escal,loc=Eloc,scale=Escal,size=1)
			dm15Boerr   = (m15Bobserr**2 + error[0]**2)**0.5
			dm15Bo      = np.random.normal(self.dm15Bs[s], dm15Boerr)
			dm15Bobs.append(np.array([dm15Bo,dm15Boerr]))

		self.errors   = errors
		self.mobs     = mobs
		self.dm15Bobs = dm15Bobs


def get_Lint_prior_samples(Nc=6,savefolder='products/Lint_sims/',generator_file='generator.stan',n_sampling=10000,n_chains=1,n_warmup=10):
	"""
	Get Lint_prior_samples

	Function to draw random correlation matrices from LKJ(1) prior using Stan

	Parameters
	----------
	Nc : int (optional; default=6)
		the number of filters (e.g. BVriJH, Nc=6)

	generator_file : str (optional; default='generator.stan')
		filename of stan to sample random correlation matrices from LKJ(1)

	n_sampling,n_chains,n_warmup : ints
		samples no.s for stan
	"""
	import stan

	Lint_sims_path = savefolder
	filename       = f'Lintsamples_Nc{Nc}.pkl'

	try:
		with open(Lint_sims_path+filename,'rb') as f:
			Lint_eta_samples = pickle.load(f)
		if len(Lint_eta_samples)<n_sampling:
			raise Exception('Generate more samples')
	except:
		with open(generator_file, 'r') as f:
			generator_file = f.read()
		posterior  = stan.build(generator_file, data={'Nc':Nc})
		fit        = posterior.sample(num_chains=n_chains, num_samples=n_sampling, num_warmup = n_warmup)
		df         = fit.to_frame()

		Lint_eta_samples = []
		for ind in df.index:
			if ind%100==0:
				print (ind,'/',n_sampling)
			L_int_eta = np.zeros((Nc,Nc))
			for i in range(Nc):
				for j in range(Nc):
					L_int_eta[i,j] = df[f'L_cint_eta.{i+1}.{j+1}'].loc[ind]
			Lint_eta_samples.append(L_int_eta)

		with open(Lint_sims_path+filename,'wb') as f:
			pickle.dump(Lint_eta_samples,f)

	return Lint_eta_samples



def get_DF_M_from_truths(truths,pblist,errstr,tref):
	"""
	Get DF_M from Truths

	Simple function that takes in simulated data,
	and outputs a DF_M that can be read by BIRDSNACK
	so as to do HBM fit

	Parameters
	----------
	truths : dict
		simulated data

	pblist : list
		list of passband strings

	errstr : str
		appended to each passband string to designate measurement error

	tref : float
		key in DF_M marking peak time

	Returns
	----------
	DF_M : dict
		{tref:df} where df is pandas df of peak magnitude measurements
	"""
	columns   = pblist+[pb+errstr for pb in pblist]
	data_rows = [np.hstack((mo,eo)) for mo,eo in zip(truths.mobs,truths.errors)]
	df        = pd.DataFrame(data=data_rows,columns=columns)

	DF_M = {tref:df}
	return DF_M
