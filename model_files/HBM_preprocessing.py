"""
HBM Preprocesser Module

Module containing HBM_preprocessor class
Used to take data products and prepare stan_data dictionry for HBM

Contains:
--------------------
HBM_preprocessor class:
	inputs: choices, DF_M

	Methods are:
		get_leff_rest()
		get_dustlaw()
		get_BV_ordered_DF_M()
		get_RVbin_data(flatRVsmin=1, flatRVsmax=6, a_sigma_beta0=1, a_sigma_beta1=1, a_sigma_muAVsigmoid=1, a_sigma_sigAVsigmoid=1)
		get_censored_data()
		get_intrinsic_model(FITSpath)
		get_dm15Bs()
		get_CM_transformation_matrix(DataTransformation=None, returner=False)
		get_CC_transformation_matrices()
		get_transformed_data()
		multiply_dataset(stan_data)
		get_init(stan_data)
		data_model_checks(stan_data)
		get_stan_file(stanpath)
		modify_stan_file()
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

import copy, pickle, re
import numpy as np
from snpy import fset
import pandas as pd
import spline_utils
from contextlib import suppress

class HBM_preprocessor:

	def __init__(self,choices,DF_M):
		"""
		Initialisation

		Parameters
		----------
		choices : dict
			.yaml dict of choices

		DF_M : dict of pandas df
			magnitudes in pblist at tilist, and also Extra Features like mB(t=15)
		"""
		self.choices  = choices
		self.DF_M     = DF_M
		self.DLAMpath = self.choices['rootpath']+'sbc/DLAMS/'

	def get_leff_rest(self):
		"""
		Get Effective Wavelengths in Rest Frame

		Method to return effective wavelengths using either passband central wavelengths, or lam_eff computed using simulations

		Returns
		----------
		l_eff_rest : array
			effective wavelengths in pblist
		"""
		pblist      = self.choices['preproc_parameters']['pblist']
		lam_choice  = self.choices['analysis_parameters']['lam_choice']

		l_eff_rest  = np.array([fset[pb].ave_wave for pb in pblist])

		if lam_choice=='central':
			l_eff_rest = np.array([fset[pb].ave_wave for pb in pblist])
			print (f'Using passband central wavelengths: {l_eff_rest}')
		else:
			if lam_choice=='effective_thetaON':
				l_eff_rest = pd.read_csv(f'{self.DLAMpath}DLAM_thetaON.txt',index_col=False)[pblist].median(axis=0).values
			elif lam_choice=='effective_thetaOFF':
				l_eff_rest = pd.read_csv(f'{self.DLAMpath}DLAM_thetaOFF.txt',index_col=False)[pblist].median(axis=0).values
			elif lam_choice=='effective_thetaOFFepsOFF':
				l_eff_rest = pd.read_csv(f'{self.DLAMpath}DLAM_thetaOFFepsOFF.txt',index_col=False)[pblist].median(axis=0).values
			print (f'Using effective wavelengths ({lam_choice}): {l_eff_rest}')

		self.l_eff_rest = l_eff_rest
		return self.l_eff_rest

	def get_dustlaw(self):
		"""
		Get Dustlaw

		Method to get dust law matrices used for constraining free RVs parameters

		End Product(s)
		----------
		self.matrice(s) : array(s)
			depends on dustlaw and data transformation
			for example,
		"""
		IntrinsicModel = self.choices['analysis_parameters']['IntrinsicModel']
		l_eff_rest = self.get_leff_rest()

		self.xk = np.array([0.0, 1e4/26500., 1e4/12200., 1e4/6000., 1e4/5470., 1e4/4670., 1e4/4110., 1e4/2700., 1e4/2600.])
		KD_x    = spline_utils.invKD_irr(self.xk)
		M_fitz_block = spline_utils.spline_coeffs_irr(1e4/l_eff_rest, self.xk, KD_x)
		if IntrinsicModel=='Deviations':
			self.M_fitz_block = M_fitz_block
		else:
			if IntrinsicModel=='Adjacent':	#Adjacent Colours
				dM_fitz_block   = np.array([ M_fitz_block[i,:]-M_fitz_block[i+1,:] for i in range(M_fitz_block.shape[0]-1) ])
			elif IntrinsicModel=='B-X':		#B-X colours
				dM_fitz_block   = np.array([ M_fitz_block[0,:]-M_fitz_block[i+1,:] for i in range(M_fitz_block.shape[0]-1) ])
			elif IntrinsicModel=='X-H':		#X-H colours
				dM_fitz_block   = np.array([ M_fitz_block[i,:]-M_fitz_block[-1,:]  for i in range(M_fitz_block.shape[0]-1) ])
			self.dM_fitz_block = dM_fitz_block

	def get_BV_ordered_DF_M(self):
		"""
		Get BV Ordered DF_M

		Method to take DF_M and order rows in increasing B-V colour

		End Product(s)/Returns
		----------
		Returns self.DF_M now ordered by B-V colour
		"""
		DF_M    = copy.deepcopy(self.DF_M)
		BVs     = DF_M[0]['B']-DF_M[0]['V']
		BVs.sort_values(ascending=True,inplace=True)
		self.BVs = BVs
		ordered_sns = list(BVs.index)#Order of SNs according to B-V

		#For each pandas df stored in DF_M dictionary
		for key in self.DF_M:
			if key not in ['extra','choices']:
				self.DF_M[key]= self.DF_M[key].loc[ordered_sns]#Re-order
			if key=='extra':
				for mini_key in self.DF_M[key]:
					self.DF_M[key][mini_key] = self.DF_M[key][mini_key].loc[ordered_sns]
		return self.DF_M

	def get_RVbin_data(self, flatRVsmin=1, flatRVsmax=6, a_sigma_beta0=1, a_sigma_beta1=1, a_sigma_muAVsigmoid=1, a_sigma_sigAVsigmoid=1):
		"""
		Get RV Bin Data

		Method takes in the B-V Colour Boundaries, and the ordered choices of RV distribution in each bin
		Then gets data products required to run Stan model fit

		Parameters
		----------
		input data: flatRVsmin=1, flatRVsmax=6, a_sigma_beta0=1, a_sigma_beta1=1, a_sigma_muAVsigmoid=1, a_sigma_sigAVsigmoid=1
			if these keys don't appear in .choices, then update to these default values

		End Product(s)
		----------
		self.N_RV_bins : int
			No. of different bins for different RV distributions

		self.N_GaussRV_dists : int
		 	No. of these Bins which have a Gaussian RV distribution (and hence No. of sets of muR,sigR hyperparameters)

		self.RV_bin_vec : array
			vector denoting Bin Number for each SN

		self.RVstyle_per_bin : list
			list maps string keys for different RV distributions to float values read by Stan
			'Gauss'->0, 'Flat'->1, 'Fix_?'->2

		self.map_RVBin_to_RVHyp : list
			for each Bin, either maps to null==-1, or the nth Gaussian RV bin

		self.fixed_RVs : array
			vector denoting value to fix RV to if required, null=-1
		"""
		#Assign default values if not already assigned
		kwargs = locals().copy()
		for key in ['flatRVsmin','flatRVsmax','a_sigma_beta0','a_sigma_beta1','a_sigma_muAVsigmoid','a_sigma_sigAVsigmoid']:
			if key not in self.choices['analysis_parameters']:
				self.choices['analysis_parameters'][key] = kwargs[key]
		#Where RVbin info will be collected
		BINS = pd.DataFrame()
		#Input Choices
		BVbinboundaries = self.choices['analysis_parameters']['BVbinboundaries']+[np.inf]
		RVstyles        = self.choices['analysis_parameters']['RVstyles']
		BVs             = self.BVs.values
		#No. of bins, and no. which are Gaussian RV dist.
		self.N_RV_bins           = int(len(BVbinboundaries))
		self.N_GaussRV_dists     = int(len([x for x in RVstyles if x=='Gauss']))
		self.N_AVRVBeta_dists    = int(len([x for x in RVstyles if x=='AVRVBeta']))
		self.N_AVRVSigmoid_dists = int(len([x for x in RVstyles if x=='AVRVSigmoid']))
		#Divide SNe into bins, ordered by BV colour
		binbool = [BVs < BVbinboundaries[0]]+[(BVs >= BVbinboundaries[i-1])&(BVs < BVbinboundaries[i]) for i in range(1, len(BVbinboundaries))]
		BINS['RV_bin_vec'] = np.select(condlist=binbool,choicelist=np.arange(self.N_RV_bins)+1)
		BINS['Input_RVstyles'] = np.select(condlist=binbool,choicelist=RVstyles)
		#Get info that will map float labels to specific RV distributions
		bin_to_style    = dict(zip(np.arange(self.N_RV_bins)+1, [x.split('_')[0] for x in RVstyles]))
		RVstylemapper   = dict(zip(['Gauss','Flat','Fix','AVRVBeta','AVRVSigmoid'],[0,1,2,3,4]))
		BINS['Stan_RVstyles']   = BINS['RV_bin_vec'].apply(lambda x : bin_to_style[x])
		BINS['RVstyle_per_bin'] = BINS['Stan_RVstyles'].apply(lambda x : RVstylemapper[x])
		self.RV_bin_vec      = BINS['RV_bin_vec'].astype(int).values
		self.RVstyle_per_bin = [RVstylemapper[bin_to_style[x]] for x in BINS['RV_bin_vec'].unique()]
		#Map bins to vector of RV hyperparameters
		map_RVBin_to_RVHyp = []; counter_gauss = 1 ; counter_AVRVBeta = 1; counter_AVRVSigmoid = 1
		for RVstyle in RVstyles:
			if RVstyle=='Gauss':
				map_RVBin_to_RVHyp.append(counter_gauss)
				counter_gauss += 1
			elif RVstyle=='AVRVBeta':
				map_RVBin_to_RVHyp.append(counter_AVRVBeta)
				counter_AVRVBeta += 1
			elif RVstyle=='AVRVSigmoid':
				map_RVBin_to_RVHyp.append(counter_AVRVSigmoid)
				counter_AVRVSigmoid += 1
			else:
				map_RVBin_to_RVHyp.append(-1)
		self.map_RVBin_to_RVHyp = map_RVBin_to_RVHyp

		#Get vector of values where RVs is fixed
		def get_fixed_RV(col):
			if 'Fix' in col: return float(col.split('_')[-1])
			else:	return -1
		BINS['fixed_RVs'] = BINS['Input_RVstyles'].apply(get_fixed_RV)
		self.fixed_RVs    = BINS['fixed_RVs'].values
		#No. of SNe in each Bin
		self.N_in_each_bin = BINS['RV_bin_vec'].value_counts().sort_index().values

	def get_censored_data(self):
		"""
		Get Censored Data

		Method to get Total No. of SNe, No. of Censored SNe, and BV errors of Censored SNe

		End Product(s)
		----------
		self.S, self.SC : integers
			Total No. of SNe (which include Censored SNe), and No. of Censored SNe

		self.BVerrs_Cens : array
			the B-V measurement errors of the Censored SNe

		self.RetainedSNe,CensoredSNe,ExcludedSNe : lists of str
			the SNe that are Retained in low-reddening sample, the censored cut, and those which are excluded from the analysis entirely
		"""
		DF_M    = copy.deepcopy(self.DF_M)
		self.S  = DF_M[0].shape[0]
		self.SC = 0
		self.RetainedSNe = list(DF_M[0].index)
		self.BVerrs_Cens = np.array([])
		self.BVcutval    = self.choices['additional_cut_parameters']['BVcutval']

		if self.choices['analysis_parameters']['CensoredData']:
			CensoredCut = self.choices['analysis_parameters']['CensoredCut']
			if CensoredCut=='inf':
				CensoredCut = np.inf
			else:
				CensoredCut = float(CensoredCut)
			errstr = self.choices['preproc_parameters']['errstr']

			DF_M[0]['BV']          =  DF_M[0]['B']-DF_M[0]['V']
			DF_M[0][f'BV{errstr}'] = (DF_M[0][f'B{errstr}']**2+DF_M[0][f'V{errstr}']**2)**0.5

			RetainedSNe = list(DF_M[0][DF_M[0]['BV'].abs()<self.BVcutval].index)
			CensoredSNe = list(DF_M[0][(DF_M[0]['BV'].abs()>=self.BVcutval) & (DF_M[0]['BV'].abs()<CensoredCut)].index)
			ExcludedSNe = list(DF_M[0][DF_M[0]['BV'].abs()>=CensoredCut].index)
			BVerrs_Cens = DF_M[0].loc[CensoredSNe][f"BV{errstr}"].values

			self.RetainedSNe = RetainedSNe
			self.BVerrs_Cens = BVerrs_Cens
			self.CensoredSNe = CensoredSNe
			self.ExcludedSNe = ExcludedSNe

			self.S  = int(len(RetainedSNe)+len(CensoredSNe))
			self.SC = int(len(CensoredSNe))

	def get_intrinsic_model(self, FITSpath):
		"""
		Get Intrinsic Model

		Method to get posterior median L_mint, FPC0 and FPC1 from a previous fit, to fix these values in a new fit

		Parameters
		----------
		FITSpath : str
			path/to/where .pkl FITs files are saved

		End Product(s)
		----------
		self.L_mint, self.FPC0, self.FPC1 : matrix, array, array
			The intrinsic hyperparameters
		"""
		#The fixint .yaml choice is the strname of the HBM_savekey
		savekey = self.choices['analysis_parameters']['fixint']
		if savekey is not False and type(savekey) is str:
			#Load up FIT
			filename = f'{FITSpath}FIT{savekey}.pkl'
			with open(filename,'rb') as f:
				FIT = pickle.load(f)
			#Ensure both old and new model is intrinsic deviations
			assert(self.choices['analysis_parameters']['IntrinsicModel']=='Deviations')
			assert(FIT['choices']['analysis_parameters']['IntrinsicModel']=='Deviations')
			#Perparations
			df     = FIT['df']
			Nm     = int(len(self.choices['preproc_parameters']['pblist']))
			L_mint = np.zeros((Nm,Nm))
			#Get posterior median intrinsic hyperparameters, and update hyperparameters
			for i in range(L_mint.shape[0]):
				for j in range(L_mint.shape[1]):
					L_mint[i,j] = df[f'L_mint.{i+1}.{j+1}'].median()
			self.L_mint = L_mint
			self.FPC0   = df[[col for col in df.columns if 'FPC0.' in col]].median(axis=0).values
			self.FPC1   = df[[col for col in df.columns if 'FPC1.' in col]].median(axis=0).values
		else:
			raise Exception('To fix the intrinsic model, set fixint in .yaml to the str name of a HBM_savekey')

	def get_dm15Bs(self):
		"""
		Get dm15Bs

		Get light curve shape parameters from DF_M, considering both:
		order of SNe if CensoredData cut applied,
		and also whether LC shape is included at all

		End Product(s)
		----------
		dm15Bs, dm15B_errs : arrays
			light curve shape parameter measurements, if not include_LCshape, set these to dummy values, and scale by gamma_shape=0 in model
		"""
		DF_M = copy.deepcopy(self.DF_M)
		include_LCshape = self.choices['analysis_parameters']['include_LCshape']
		errstr = self.choices['preproc_parameters']['errstr']
		def get_dm15B(DF_M,sn):
			dm15B    =  DF_M['extra'][15]['B'].loc[sn]-DF_M[0]['B'].loc[sn]
			dm15Berr = (DF_M['extra'][15][f'B{errstr}'].loc[sn]**2+DF_M[0][f'B{errstr}'].loc[sn]**2)**0.5
			return dm15B, dm15Berr

		if include_LCshape:
			dm15Bs = [] ; dm15B_errs = []
			if not self.choices['analysis_parameters']['CensoredData']:
				ordered_sns = list(DF_M[0].index)
			elif self.choices['analysis_parameters']['CensoredData']:
				ordered_sns = self.RetainedSNe+self.CensoredSNe
			#Append dm15Bs in correct order, Retained SNe then Censored SNe
			for sn in ordered_sns:
				dm15B,dm15Berr = get_dm15B(DF_M,sn)
				dm15Bs.append(dm15B)
				dm15B_errs.append(dm15Berr)
		else:
			dm15Bs     = np.ones(self.S)
			dm15B_errs = np.ones(self.S)

		self.dm15Bs     = dm15Bs
		self.dm15B_errs = dm15B_errs

	def get_CM_transformation_matrix(self, DataTransformation=None, returner=False):
		"""
		Get CM Transformation Matrix, where c = CM @ m

		Parameters
		----------
		DataTransformation : str (optional; default=None)
			defines CM transformation, if None, use the entry in self.choices

		returner : bool (optional; default=False)
			if True, return and don't set attribute
			if False, set attribute and don't return

		Returns/End Product(s)
		----------
		CM : matrix
			transforms vector of mags to vector of colours
		"""
		if DataTransformation is None:
			DataTransformation = self.choices['analysis_parameters']['DataTransformation']
		Nc = int(len(self.choices['preproc_parameters']['pblist'])-1)
		CM = np.zeros((Nc,Nc+1))
		if DataTransformation=='Adjacent':
			for _ in range(Nc):
				CM[_,_]   =  1
				CM[_,_+1] = -1
		elif DataTransformation=='B-X':
			for _ in range(Nc):
				CM[_,0]   =  1
				CM[_,_+1] = -1
		elif DataTransformation=='X-H':
			for _ in range(Nc):
				CM[_,Nc]  = -1
				CM[_,_]   =  1
		if returner:
			return CM
		else:
			self.CM = CM

	def get_CC_transformation_matrices(self):
		"""
		Get Colour_to_Colour Tranformation Matrices

		End Product(s)
		----------
		self.N : int
			No. of observations
		self.CC : matrix
			transformation from model colours to data colours
		self.CC_to_adj : matrix
			transformation from model colours to adjacent colours (where 1st adj colour is B-V, for use in Censoring)
		"""
		#No. of Passbands
		Nm    = int(len(self.choices['preproc_parameters']['pblist']))
		#Different Transformations from magnitudes to colours
		CMadj = self.get_CM_transformation_matrix(DataTransformation='Adjacent',returner=True)
		CMBX  = self.get_CM_transformation_matrix(DataTransformation='B-X', 	returner=True)
		CMXH  = self.get_CM_transformation_matrix(DataTransformation='X-H', 	returner=True)
		#Pad with ones to create inverse
		CMadj = np.vstack((CMadj,np.ones(Nm))) ; CMBX  = np.vstack((CMBX,np.ones(Nm))) ; CMXH  = np.vstack((CMXH,np.ones(Nm)))
		CMs     = dict(zip(['Adjacent','B-X','X-H'],[CMadj,CMBX,CMXH]))
		#CM_data is the Transformation from Mags to the Data being fitted ; CM_model is the Transformation from Mags to the Intrinsic Colours Model
		CM_data  = CMs[self.choices['analysis_parameters']['DataTransformation']]
		CM_model = CMs[self.choices['analysis_parameters']['IntrinsicModel']]
		#Build matrices to transform from one set of colours to another ; Create transformation to adjacent colours so CenosoredData incorporate in Model (i.e. B-V = AdjColour[0])
		CC        = (CM_data@np.linalg.inv(CM_model))[:-1,:-1]
		CC_to_adj = (CMs['Adjacent']@np.linalg.inv(CM_model))[:-1,:-1]
		CCS       = {'CC':CC,'CC_to_adj':CC_to_adj}
		#Sometimes transformations aren't exact, round small number deviations down to makes 0's or 1's
		for ccc in CCS:
			CCC = CCS[ccc]
			for i in range(CCC.shape[0]):
				for j in range(CCC.shape[1]):
					if abs(CCC[i,j]-1)<1e-10: 	CCC[i,j] = int(1)
					elif abs(CCC[i,j])<1e-10:	CCC[i,j] = int(0)
					else: pass
			CCS[ccc] = CCC

		self.CC         = CCS['CC']
		self.CC_to_adj  = CCS['CC_to_adj']

	def get_transformed_data(self):
		"""
		Get Transformed Data

		Takes in DF_M, and outputs vectors of magnitudes or colours, and measurement errors (for colours, this is a covariance matrix of errors)

		End Product(s)
		----------
		self.mags, mags_errs: arrays
		OR
		self.capps, capps_errs : vector, vector of matrices respectively
		"""
		DF_M   = copy.deepcopy(self.DF_M)
		pblist = self.choices['preproc_parameters']['pblist']
		tref   = self.choices['preproc_parameters']['tilist'][self.choices['preproc_parameters']['tref_index']]
		errstr = self.choices['preproc_parameters']['errstr']
		DataTransformation = self.choices['analysis_parameters']['DataTransformation']
		Nm     = int(len(pblist))

		mags       = DF_M[tref][[pb for pb in pblist]].loc[self.RetainedSNe].stack().transpose().values
		mags_errs  = DF_M[tref][[pb+errstr for pb in pblist]].loc[self.RetainedSNe].stack().transpose().values
		self.mags = mags
		self.mags_errs = mags_errs
		if DataTransformation!='mags':
			capps,capps_errs = [],[]
			for s in range(len(self.RetainedSNe)):
				capps.extend(self.CM @ mags[s*Nm:(s+1)*Nm])
				data_cov_matrix = self.CM @ np.diag( np.asarray( mags_errs[s*Nm:(s+1)*Nm] )**2 ) @ self.CM.T
				capps_errs.append(data_cov_matrix)
			self.capps      = capps
			self.capps_errs = capps_errs

	def multiply_dataset(self, stan_data):
		"""
		Mulitply Dataset

		Make N=copymode exact replicas of sample, and fit that sample

		Parameters
		----------
		stan_data : dict
			dictionary of input data

		Returns
		----------
		stan_data with additional copies
		"""
		copymode = self.choices['analysis_parameters']['copymode']
		if copymode is not False:
			if self.choices['analysis_parameters']['DataTransformation']!='mags':
				raise Exception('Not yet built code for replicating colour data')
			stan_data['S']  = int(stan_data['S'] *(1+copymode))
			stan_data['SC'] = int(stan_data['SC']*(1+copymode))
			OG = {key:copy.deepcopy(stan_data[key]) for key in ['mags','mags_errs','dm15Bs','dm15B_errs','BVerrs_Cens']}
			for _ in range(copymode):
				for key in OG:
					stan_data[key] = np.concatenate((stan_data[key],OG[key]))
		return stan_data

	def get_init(self, stan_data):
		"""
		Get Stan Initialisations

		Parameters
		----------
		stan_data : dict
			dictionary of data

		Returns
		----------
		stan_init : dict
			dictionary of initialisations
		"""
		n_chains = self.choices['analysis_parameters']['n_chains']
		if self.choices['analysis_parameters']['IntrinsicModel']=='Deviations':
			mu_guesses = np.array([np.average(mi) for mi in np.array_split(self.mags,stan_data['S']-stan_data['SC'])])
			stan_init  = [{'mus':np.random.normal(mu_guesses,0.5)} for _ in range(n_chains)]
		else:
			stan_init  = [{} for _ in range(n_chains)]
		return stan_init

	def data_model_checks(self,stan_data):
		"""
		Data Model Checks

		Simple method to assert lengths of data vectors match those asserted by S,SC and Nm integers
		Also removes any data not required by stan file

		Parameters
		----------
		stan_data : dict
			dictionary of input data

		Returns
		----------
		stan_data with any superfluous data removed
		"""
		if self.choices['analysis_parameters']['DataTransformation']=='mags':
			assert(len(stan_data['mags'])==(stan_data['S']-stan_data['SC'])*stan_data['Nm'])
			assert(len(stan_data['mags_errs'])==(stan_data['S']-stan_data['SC'])*stan_data['Nm'])
		else:
			assert(len(stan_data['capps'])==(stan_data['S']-stan_data['SC'])*stan_data['Nc'])
			assert(len(stan_data['capps_errs'])==(stan_data['S']-stan_data['SC']))
		assert(len(stan_data['dm15Bs'])==stan_data['S'])
		assert(len(stan_data['dm15B_errs'])==stan_data['S'])
		assert(len(stan_data['BVerrs_Cens'])==stan_data['SC'])

		if 'skew_RV' in self.choices['analysis_parameters'] and self.choices['analysis_parameters']['skew_RV']:
			assert(self.choices['analysis_parameters']['RVprior']=='Norm')

		DataTransformation = self.choices['analysis_parameters']['DataTransformation']
		IntrinsicModel     = self.choices['analysis_parameters']['IntrinsicModel']
		stan_data.pop('RVsmax')
		if IntrinsicModel=='Deviations':
			stan_data.pop('a_sigma_cint')
			if DataTransformation=='mags':
				stan_data.pop('Nc')
		else:
			stan_data.pop('a_sigma_mint')
			stan_data.pop('Nm')
			stan_data.pop('zero_index')
		return stan_data

	def get_stan_file(self,stanpath):
		"""
		Get Stan File

		Get path/to/filename.stan of Stan Model File

		Parameters
		----------
		stanpath : str
			'path/to/folder/with/stan_files/'

		End Product(s)
		----------
		self.stan_file
		"""
		#Check whether to use stan file with skewint, skewRV, studenttRV
		use_skew_stan_file = False
		with suppress(KeyError):
			if self.choices['analysis_parameters']['RVprior']!='Norm': 	use_skew_stan_file = True
		with suppress(KeyError):
			if self.choices['analysis_parameters']['skew_RV']: 			use_skew_stan_file = True
		with suppress(KeyError):
			if self.choices['analysis_parameters']['skew_int']: 	   	use_skew_stan_file = True
		#This more complicated file is only valid for Devations model fitted to magnitudes data
		if use_skew_stan_file:
			assert(self.choices['analysis_parameters']['IntrinsicModel']=='Deviations')
			assert(self.choices['analysis_parameters']['DataTransformation']=='mags')
		#Get stan file
		if self.choices['analysis_parameters']['IntrinsicModel']=='Deviations':
			if self.choices['analysis_parameters']['DataTransformation']=='mags':
				print ("Applying Intrinsic Deviations Model to fit Apparent Magnitudes Data") #stan_file = f"{stanpath}deviations_model_fit_mags_Dirichletmuint.stan"
				if use_skew_stan_file:	stan_file = f"{stanpath}Extra/deviations_model_fit_mags_Gaussianmuintref_Skew.stan"
				elif 'BinnedRVFit' in self.choices['analysis_parameters'] and self.choices['analysis_parameters']['BinnedRVFit']:
					if 'AVRVBeta' in self.choices['analysis_parameters']['RVstyles']:#NOTE, HC(0,1) Hyperparameters are reason for files not running, so can either score these parameters out when not in use, OR (as done here) just load in different stan files
						stan_file = f"{stanpath}Extra/deviations_model_fit_mags_Gaussianmuintref_RVBins_wAVRVBeta.stan"
					elif 'AVRVSigmoid' in self.choices['analysis_parameters']['RVstyles']:
						stan_file = f"{stanpath}Extra/deviations_model_fit_mags_Gaussianmuintref_RVBins_wAVRVSigmoid.stan"
					else:
						stan_file = f"{stanpath}Extra/deviations_model_fit_mags_Gaussianmuintref_RVBins.stan"
				elif self.choices['analysis_parameters']['fixint'] is not False:
					if self.choices['analysis_parameters']['AVprior']=='Gauss':
						stan_file = f"{stanpath}Extra/deviations_model_fit_mags_Gaussianmuintref_fixedint_GaussAV.stan"
					else:
						stan_file = f"{stanpath}Extra/deviations_model_fit_mags_Gaussianmuintref_fixedint.stan"
				else:	stan_file = f"{stanpath}deviations_model_fit_mags_Gaussianmuintref.stan"
			else:
				print (f"Applying Intrinsic Deviations Model to fit {self.choices['analysis_parameters']['DataTransformation']} Colours Data")
				stan_file = f"{stanpath}deviations_model_fit_colours_Gaussianmuintref.stan"
		else:
			print (f"Applying {self.choices['analysis_parameters']['IntrinsicModel']} Intrinsic Colours Model to fit {self.choices['analysis_parameters']['DataTransformation']} Colours Data")
			stan_file = f"{stanpath}colours_model.stan"

		self.stan_file = stan_file
		self.use_skew_stan_file = use_skew_stan_file

	def modify_stan_file(self):
		"""
		Modify Stan File

		Modify .stan file based on analysis choices

		End Product(s)
		----------
		stan_file: str
			long string of .stan code read by Stan functions
		"""

		stan_file = copy.deepcopy(self.stan_file)
		AVprior,RVprior,skew_RV,skew_int='Exp','Norm',False,False
		edit_file = False
		AVprior   = self.choices['analysis_parameters']['AVprior']
		if AVprior!='Exp':#If AV prior is not exponential on tau
			edit_file = True
		with suppress(KeyError):#If RV prior is not Normal
			RVprior   = self.choices['analysis_parameters']['RVprior']
			if RVprior!='Norm':	edit_file = True
		with suppress(KeyError):#If RVs distribution is skewed
			skew_RV   = self.choices['analysis_parameters']['skew_RV']
			if skew_RV:			edit_file = True
		with suppress(KeyError):#If Intrinsic Deviations are skewed
			skew_int  = self.choices['analysis_parameters']['skew_int']
			if skew_int:		edit_file = True

		if edit_file:
			with open(stan_file,'r') as f:
				stan_file = f.read().splitlines()
			for il,line in enumerate(stan_file):
				###Model Block
				#Get rid of exponential prior
				if AVprior!='Exp':
					if bool(re.match(r'\s*AVs\s*~\s*exponential\(inv\(tauA\)\)',line)):
						stan_file[il] = line.replace('AVs','//AVs')
				#Introduce gamma prior
				if AVprior == 'Gamma':
					if bool(re.match(r'\s*//\s*AVs\s*~\s*gamma\(nu,1/tauA\)',line)):
						stan_file[il] = line.replace('//','')

				if AVprior in ['Gamma']:#Introduce parameter nu in parameters block and prior on nu
					###Prior Block
					if bool(re.match(r'\s*//\s*nu_tform\s*~\s*uniform\(0,pi\(\)/2\)',line)):
						stan_file[il] = line.replace('//','')
					###Parameters Block
					if bool(re.match(r'\s*//\s*real<lower=0,upper=pi\(\)/2>\s*nu_tform',line)):
						stan_file[il] = line.replace('//','')
					###Transformed Parameters Block
					if bool(re.match(r'\s*//\s*real<lower=0>\s*nu;',line)):
						stan_file[il] = line.replace('//','')
					if bool(re.match(r'\s*//\s*nu\s*=\s*tan\(nu_tform\)',line)):
						stan_file[il] = line.replace('//','')

				if skew_RV or RVprior in ['StudentT']:
					###Data Block
					if skew_RV:
						if bool(re.match(r'\s*//real\s*skew_RV_disp_prior',line)):
							stan_file[il] = line.replace('//real','real')
					###Parameters Block
					#Remove this
					if bool(re.match(r'\s*vector<lower=0,upper=1>\[S\]\s*nuRVs',line)):
						stan_file[il] = line.replace('vector','//vector')
					#Add these
					if bool(re.match(r'\s*//vector<lower=RVsmin>\[S\]\s*RVs',line)) and 'Rescaled Parameters' not in stan_file[il-1]:
						stan_file[il] = line.replace('//vector','vector')
					if skew_RV:
						if bool(re.match(r'\s*//real\s*alpha_skew_RV',line)):
							stan_file[il] = line.replace('//real','real')
					if RVprior=='StudentT':
						if bool(re.match(r'\s*//real<lower=0>\s*nuR',line)):
							stan_file[il] = line.replace('//real','real')
					###Transformed Parameters Block
					#Remove these
					if bool(re.match(r'\s*vector<lower=RVsmin>\[S\]\s*RVs',line)) and 'Rescaled Parameters' in stan_file[il-1]:
						stan_file[il] = line.replace('vector','//vector')
					if bool(re.match(r'\s*real<upper=0>\s*alpha',line)):
						stan_file[il] = line.replace('real','//real')
					if bool(re.match(r'\s*alpha\s*=\s*\(RVsmin\s*-\s*mu_RV\s*\)/sig_RV',line)):
						stan_file[il] = line.replace('alpha','//alpha')
					if bool(re.match(r'\s*RVs\s*=\s*mu_RV\s*-\s*sig_RV\s*\*\s*\(\s*inv_Phi\s*\(\s*nuRVs\s*\*\s*Phi\s*\(-alpha\)\s*\)\s*\)',line)):
						stan_file[il] = line.replace('RVs','//RVs')
					###Model Block
					#Remove this
					if bool(re.match(r'\s*nuRVs\s*~\s*uniform\(0,1\)',line)):
						stan_file[il] = line.replace('nuRVs','//nuRVs')
					#Add these
					if skew_RV:
						if bool(re.match(r'\s*//RVs\s*~\s*skew_normal\(mu_RV,\s*sig_RV,\s*alpha_skew_RV\s*\)',line)):
							stan_file[il] = line.replace('//RVs','RVs')
						if bool(re.match(r'\s*//alpha_skew_RV\s*~\s*normal\(0,skew_RV_disp_prior\)',line)):
							stan_file[il] = line.replace('//alpha_skew','alpha_skew')
					if RVprior=='StudentT':
						if bool(re.match(r'\s*//RVs\s*~\s*student_t\(nuR,\s*mu_RV,\s*sig_RV\s*\)',line)):
							stan_file[il] = line.replace('//RVs','RVs')
						if bool(re.match(r'\s*//nuR\s*~\s*gamma\(2,0.1\)',line)):
								stan_file[il] = line.replace('//nuR','nuR')
				if skew_int:
					###Data Block
					if bool(re.match(r'\s*//real\s*skew_int_disp_prior',line)):
						stan_file[il] = line.replace('//real','real')
					###Parameters Block
					if bool(re.match(r'\s*//real\s*alpha_skew_int',line)):
						stan_file[il] = line.replace('//real','real')
					###Model Block
					#Remove These
					if bool(re.match(r'\s*eta_mint\s*~\s*std_normal\(\)',line)):
						stan_file[il] = line.replace('eta_mint','//eta_mint')
					if bool(re.match(r'\s*eta_mint_Cens\s*~\s*std_normal\(\)',line)):
						stan_file[il] = line.replace('eta_mint_Cens','//eta_mint_Cens')
					if bool(re.match(r'\s*//alpha_skew_int\s*~\s*normal\(0,skew_int_disp_prior\)',line)):#Introduce skewness
						stan_file[il] = line.replace('//alpha_skew','alpha_skew')
					#Add these
					if bool(re.match(r'\s*//eta_mint\s*~\s*skew_normal\(0,1,alpha_skew_int\)',line)):
						stan_file[il] = line.replace('//eta_mint','eta_mint')
					if bool(re.match(r'\s*//eta_mint_Cens\s*~\s*skew_normal\(0,1,alpha_skew_int\)',line)):
						stan_file[il] = line.replace('//eta_mint_Cens','eta_mint_Cens')




			#Stan Reads \n delimited text
			stan_file = '\n'.join(stan_file)
		else:
			#Read in as \n delimited text
			with open(stan_file, 'r') as f:
				stan_file = f.read()
		#Print and return (edited) model file
		print (stan_file)
		self.stan_text_file = stan_file
		return self.stan_text_file
