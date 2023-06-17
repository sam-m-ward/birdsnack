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
		get_censored_data()
		get_dm15Bs()
		get_CM_transformation_matrix(DataTransformation=None, returner=False)
		get_CC_transformation_matrices()
		get_transformed_data()
		multiply_dataset(stan_data)
		data_model_checks(stan_data)
		get_stan_file(stanpath)
		modify_stan_file()
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

import copy
import numpy as np
from snpy import fset
import pandas as pd
import spline_utils

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
		DataTransformation = self.choices['analysis_parameters']['DataTransformation']
		l_eff_rest = self.get_leff_rest()

		self.xk = np.array([0.0, 1e4/26500., 1e4/12200., 1e4/6000., 1e4/5470., 1e4/4670., 1e4/4110., 1e4/2700., 1e4/2600.])
		KD_x    = spline_utils.invKD_irr(self.xk)
		if DataTransformation=='mags':
			M_fitz_block = spline_utils.spline_coeffs_irr(1e4/l_eff_rest, self.xk, KD_x)
			self.M_fitz_block = M_fitz_block
		else:
			if DataTransformation=='Adjacent':	#Adjacent Colours
				dM_fitz_block   = np.array([ M_fitz_block[i,:]-M_fitz_block[i+1,:] for i in range(M_fitz_block.shape[0]-1) ])
			elif DataTransformation=='B-X':		#B-X colours
				dM_fitz_block   = np.array([ M_fitz_block[0,:]-M_fitz_block[i+1,:] for i in range(M_fitz_block.shape[0]-1) ])
			elif DataTransformation=='X-H':		#X-H colours
				dM_fitz_block   = np.array([ M_fitz_block[i,:]-M_fitz_block[-1,:]  for i in range(M_fitz_block.shape[0]-1) ])
			self.M_fitz_block = dM_fitz_block

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

		if self.choices['analysis_parameters']['CensoredData']:
			CensoredCut = self.choices['analysis_parameters']['CensoredCut']
			if CensoredCut=='inf':
				CensoredCut = np.inf
			else:
				CensoredCut = float(CensoredCut)
			errstr = self.choices['preproc_parameters']['errstr']

			DF_M[0]['BV']          =  DF_M[0]['B']-DF_M[0]['V']
			DF_M[0][f'BV{errstr}'] = (DF_M[0][f'B{errstr}']**2+DF_M[0][f'V{errstr}']**2)**0.5

			RetainedSNe = list(DF_M[0][DF_M[0]['BV'].abs()<0.3].index)
			CensoredSNe = list(DF_M[0][(DF_M[0]['BV'].abs()>=0.3) & (DF_M[0]['BV'].abs()<CensoredCut)].index)
			ExcludedSNe = list(DF_M[0][DF_M[0]['BV'].abs()>=CensoredCut].index)
			BVerrs_Cens = DF_M[0].loc[CensoredSNe][f"BV{errstr}"].values

			self.RetainedSNe = RetainedSNe
			self.BVerrs_Cens = BVerrs_Cens
			self.CensoredSNe = CensoredSNe
			self.ExcludedSNe = ExcludedSNe

			self.S  = int(len(RetainedSNe)+len(CensoredSNe))
			self.SC = int(len(CensoredSNe))

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

		self.N          = (stan_data['S']-stan_data['SC'])*(Nm-1)
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

	def data_model_checks(self,stan_data):
		"""
		Data Model Checks

		Simple method to assert lengths of data vectors match those asserted by S,SC and Nm integers

		Parameters
		----------
		stan_data : dict
			dictionary of input data
		"""
		if self.choices['analysis_parameters']['IntrinsicModel']=='Deviations':
			assert(len(stan_data['mags'])==(stan_data['S']-stan_data['SC'])*stan_data['Nm'])
			assert(len(stan_data['mags_errs'])==(stan_data['S']-stan_data['SC'])*stan_data['Nm'])
		assert(len(stan_data['dm15Bs'])==stan_data['S'])
		assert(len(stan_data['dm15B_errs'])==stan_data['S'])
		assert(len(stan_data['BVerrs_Cens'])==stan_data['SC'])

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
		if self.choices['analysis_parameters']['IntrinsicModel']=='Deviations':
			if self.choices['analysis_parameters']['DataTransformation']=='mags':
				print ("Applying Intrinsic Deviations Model to fit Apparent Magnitudes Data")
				stan_file = f"{stanpath}deviations_model_fit_mags_Gaussianmuintref.stan"
				#stan_file = f"{stanpath}deviations_model_fit_mags_Dirichletmuint.stan"
			else:
				print (f"Applying Intrinsic Deviations Model to fit {self.choices['analysis_parameters']['DataTransformation']} Colours Data")
				stan_file = f"{stanpath}mag_model_fit_colours_indRVs_popdist_FullModel_F99.stan"
		else:
			print (f"Applying {self.choices['analysis_parameters']['IntrinsicModel']} Intrinsic Colours Model to fit {self.choices['analysis_parameters']['DataTransformation']} Colours Data")
			stan_file = f"{stanpath}colour_model_indRVs_popdist_FullModel_F99.stan"

		self.stan_file = stan_file

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
		AVprior   = self.choices['analysis_parameters']['AVprior']

		if AVprior!='Exp':#If AV prior is not exponential on tau
			with open(stan_file,'r') as f:
				stan_file = f.read().splitlines()
			for il,line in enumerate(stan_file):
				###Model Block
				#Get rid of exponential prior
				if bool(re.match(r'\s*AVs\s*~\s*exponential\(inv\(tauA\)\)',line)):
					stan_file[il] = line.replace('AVs','//AVs')
				#Introduce gamma prior
				if AVprior == 'Gamma':
					if bool(re.match(r'\s*//\s*AVs\s*~\s*gamma\(nu,1/tauA\)',line)):
						stan_file[il] = line.replace('//','')

				###Parameters and Prior Blocks
				#Introduce parameter nu in parameters block and prior on nu
				if AVprior in ['Gamma']:
					#Prior
					if bool(re.match(r'\s*//\s*nu_tform\s*~\s*uniform\(0,pi\(\)/2\)',line)):
						stan_file[il] = line.replace('//','')
					#Parameters Block
					if bool(re.match(r'\s*//\s*real<lower=0,upper=pi\(\)/2>\s*nu_tform',line)):
						stan_file[il] = line.replace('//','')
					#Transformed Parameters Block
					if bool(re.match(r'\s*//\s*real<lower=0>\s*nu;',line)):
						stan_file[il] = line.replace('//','')
					if bool(re.match(r'\s*//\s*nu\s*=\s*tan\(nu_tform\)',line)):
						stan_file[il] = line.replace('//','')

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
