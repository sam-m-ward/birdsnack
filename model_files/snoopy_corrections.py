"""
SNooPy Corrections Module

Module containing SNOOPY_CORRECTIONS class and additional functions
Used to transform snpy object with corrections, saves Tmax_Kcorrdict, snpy_product, and snana lc files


Contains:
--------------------
get_Tmax_Kcorrdict_foldername(snpy_params):
	get folder name for Tmax and Kcorrections dictionary, pre-computed for all analysis variants in one go to save on time complexity

SNOOPY_CORRECTIONS class:
	inputs: choices, outputdir

	Methods are:
		get_fitted_snpy_Tmax_Kcorrdict(fitbands,obs_to_rest)
		correct_data_with_snpy_MWayExtinction_and_K_corrections()
		get_sn_paths()
		load_snpy_product(sn,survey)
		correct_sn(sn, SNsnpy)
		convert_snsnpy_to_snana()
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""
import copy
import numpy as np
from snpy import fset
import extinction
import os
import pickle
from miscellaneous import ensure_folders_to_file_exist
from snana_functions import get_snana_foldername
import snanaio as io

def get_Tmax_Kcorrdict_foldername(snpy_params):
	"""
	Get Tmax Kcorrdict Folder Name

	Function to get folder name for Tmax and Kcorrections dictionary

	Parameters
	----------
	snpy_params: dict
	Choices for SNooPy analysis
	  keys are [snrcut,snpy_model,snpy_shape,apply_EBVMW_corrections,RVMW,dustlaw,insert_pre_defined_Tmax,apply_K_corrections,mangle,interpflts]

	 Returns
	 ----------
	 str:
	 	foldername, used for multiple analysis variants, different sets of Kcorrections
	"""
	appender  = 'Tmax_Kcorrdicts_'

	#SNR Choices
	appender += f"apply_SNR_cut{snpy_params['apply_SNR_cut']}_"
	if snpy_params['apply_SNR_cut']:
		appender += f"snrcut{snpy_params['snrcut']}_"

	#Choice of SNooPy model
	appender += f"snpymodel{snpy_params['snpy_model']}_snpyshape{snpy_params['snpy_shape']}_"

	#Choice of Milky Way Extinction Corrections
	appender += f"apply_EBVMW_corrections{snpy_params['apply_EBVMW_corrections']}_"
	if snpy_params['apply_EBVMW_corrections']=='Presubtract':#Only if we are pre-subtracting out MWay Extinction does RVMW/dustlaw matter
		appender += f"RVMW{snpy_params['RVMW']}_dustlaw{snpy_params['dustlaw']}_"

	#Choice of Tmax
	if snpy_params['insert_pre_defined_Tmax']:#If Tmax is fixed at value during fit (e.g. from previous .fitMCMC())
		appender += f"predefinedTmax{snpy_params['insert_pre_defined_Tmax']}_"

	return appender[:-1]+'/'

class SNOOPY_CORRECTIONS:

	def get_fitted_snpy_Tmax_Kcorrdict(self,fitbands,obs_to_rest):
		"""
		Get Fitted SNooPy Tmax/Kcorrection Dictionary

		SNooPy Fitted Tmax and Kcorrections are uniquely defined by SNRcut, EBVMW correction, and the passbands used to perform SNooPy fit
		Thus compute this dictionary separately to avoid repeat computations

		Parameters
		----------
		fitbands: list
			list of SNooPy observer frame passbands, specific to each SN, that fall inside interpflts str

		obs_to_rest: dict
			keys are the fitbands observer frame passbands, values are the rest frame passbands each obs-pb maps to

		End Product(s)
		----------
		Tmax_Kcorrdict: dict
			Tmax_Kcorrdict = {'Tmax_snpy_fitted':None,'Mangle':{},'NoMangle':{},'snpy_cut_bools':None,'MWCorrections':{}}
			In each Kcorr dict is key,value pair, key is flt, value is snsnpy[filt].ks
			snpy_cut_bools is dictionary of whether sn was cut
			MWCorrections is integral over passband correcting for MWay extinction (gets applied when we dont pre-subtract MWay)

		"""
		obs_to_rest_str = '{'+','.join([f"{key}->{value}" for key,value in obs_to_rest.items()])+'}'
		filename        = f"{self.Tmax_Kcorrdict_folder}{self.snsnpy.name}_{self.survey}_{obs_to_rest_str}.pkl"
		snpy_params     = self.snpychoices
		try:
			with open(filename,'rb') as f:
				Tmax_Kcorrdict = pickle.load(f)

		except:
			snpy_cut_bools = {'Mangle':False,'NoMangle':False,'MWay':False}
			Tmax_Kcorrdict = {'Tmax_snpy_fitted':None,'Mangle':{},'NoMangle':{},'snpy_cut_bools':None,'MWCorrections':{}}

			try:
				fitted_snsnpy = copy.deepcopy(self.snsnpy)

				if snpy_params['insert_pre_defined_Tmax']:
					fitted_snsnpy.fit(fitbands,Tmax=fitted_snsnpy.Tmax)
				else:
					fitted_snsnpy.fit(fitbands)
				#Get kcorrected Tmax; Use SNooPy fit
				Tmax_Kcorrdict['Tmax_snpy_fitted'] = copy.deepcopy(fitted_snsnpy.Tmax)


			except Exception as e:
				print ('Fit fail:',e)
				snpy_cut_bools['Mangle']   = True
				snpy_cut_bools['NoMangle'] = True

			if snpy_params['apply_EBVMW_corrections']=='Model':
				try:
					#Get data without MWay Extinction
					no_mw_fitted_snsnpy        = copy.deepcopy(fitted_snsnpy)
					no_mw_fitted_snsnpy.EBVgal = 0
					for flt in fitbands:
						mags_no_MW_eval, success             = no_mw_fitted_snsnpy.data[flt].eval(no_mw_fitted_snsnpy.data[flt].MJD)
						mags_w_MW_eval, success              =       fitted_snsnpy.data[flt].eval(      fitted_snsnpy.data[flt].MJD)
						Tmax_Kcorrdict['MWCorrections'][flt] = mags_no_MW_eval-mags_w_MW_eval

				except Exception as e:
					print ('MWay Correction fail:',e)
					snpy_cut_bools['MWay']  = True

			mangle_snsnpy    = copy.deepcopy(fitted_snsnpy)
			no_mangle_snsnpy = copy.deepcopy(fitted_snsnpy)

			try:#Get K-corrections with Mangling
				mangle_snsnpy.kcorr(fitbands,mangle=True)
				for filt in mangle_snsnpy.ks:
					Tmax_Kcorrdict['Mangle'][filt] = mangle_snsnpy.ks[filt]
			except:
				snpy_cut_bools['Mangle'] = True
			try:#Get K-corrections without Mangling
				no_mangle_snsnpy.kcorr(mangle=False)
				for filt in no_mangle_snsnpy.ks:
					Tmax_Kcorrdict['NoMangle'][filt] = no_mangle_snsnpy.ks[filt]
			except:
				snpy_cut_bools['NoMangle'] = True

			Tmax_Kcorrdict['snpy_cut_bools'] = snpy_cut_bools

			with open(filename,'wb') as f:
				pickle.dump(Tmax_Kcorrdict,f)

		self.Tmax_Kcorrdict = Tmax_Kcorrdict

	def correct_data_with_snpy_MWayExtinction_and_K_corrections(self):
		"""
		Correct Data using SNooPy for Milky Way Extinction and Redshifting

		Returns corrected snpy object, and dictionary of important snpy values e.g. kcorrections, rest-frame filters, Tmaxs

		End product(s)
		----------
		snsnpy: snpy.sn.sn object
			snpy object of particular sn (corrected for MW-Extinction and Kcorrections)

		snpy_product: dict
			keys are ['kcorrections', 'snpy_cut', 'obsbands', 'bandmap', 'Tmax_snpy_fitted']
			kcorrections are dict={flt:np.array of kcorrections}
			snpy_cut is bool, if True, SNooPy procedure failed
			obsbands and bandmap are SNooPy filters, (key,value) pair are obs and rest bands
			Finally, Tmax float values from SNooPy fits
		"""
		snpy_params = self.snpychoices

		#Apply SNR Masking
		if snpy_params['apply_SNR_cut']:
			self.snsnpy.mask_SNR(snpy_params['snrcut'])

		#Change snpy model if required
		if snpy_params['snpy_model'] and snpy_params['snpy_shape']:#Default is simply False and False
			print ('Applying .choose_model')
			self.snsnpy.choose_model(snpy_model, stype=snpy_shape)

		def return_correct_restbands(snsnpy):
			"""
			Return Correct RestBands

			Function to map observer-frame bands to rest-frame bands
			Implements x2 Rules to boost sample size

			Parameters
			----------
			snsnpy: snpy object

			Returns
			----------
			snsnpy with updated restbands mapping dictionary
			"""
			MAPPER0 = dict(zip(['U','B','G','V','R','I','Y','J','H','K'],['u','B','g','V','r','i','Y','J','H','K']))

			#INCLUDED RULE 1
			#If J_K or H_K data available at all, transform to J,H
			for obs in snsnpy.restbands:
				rest = snsnpy.restbands[obs]
				if rest in ['J_K','H_K']:
					snsnpy.restbands[obs] = MAPPER0[rest[0].upper()]

			#INCLUDED RULE 2
			#If any of BVri not currently in restbands, then map Bs, Vs, Rs, Is to BVri
			not_currently_present = [f for f in 'BVri' if f not in list(snsnpy.restbands.values())]
			for obs,rest in dict(snsnpy.restbands).items():
				if rest in ['Bs','Vs','Rs','Is']:
					new_rest = MAPPER0[rest[0].upper()]
					if new_rest in not_currently_present:
						snsnpy.restbands[obs] = new_rest
			return snsnpy

		#Map all bands to rest-frame uBgVriYJH;
		self.snsnpy = return_correct_restbands(self.snsnpy)

		#If we are pre-subtracting out the MWay Extinction (and thus adopting the SED approximation)
		if snpy_params['apply_EBVMW_corrections']=='Presubtract':
			####################################
			EBVMW = copy.deepcopy(self.snsnpy.EBVgal)
			RVMW  = snpy_params['RVMW']
			#MWay Extinction Correction;
			#Apply this to observer-frame magnitude data
			lambdaC = {filt:fset[filt].ave_wave for filt in self.snsnpy.allbands()}

			if dustlaw=='fitzpatrick99':
				Alams   = dict( zip(lambdaC.keys(),extinction.fitzpatrick99( np.asarray(list(lambdaC.values())), EBVMW*RVMW, RVMW)))
			elif dustlaw=='ccm89':
				Alams   = dict( zip(lambdaC.keys(),extinction.ccm89(         np.asarray(list(lambdaC.values())), EBVMW*RVMW, RVMW)))
			elif dustlaw=='ccm89_o94':
				Alams   = dict( zip(lambdaC.keys(),extinction.odonnell94(    np.asarray(list(lambdaC.values())), EBVMW*RVMW, RVMW)))

			for filt in Alams:
				self.snsnpy.data[filt].mag += -Alams[filt]
			####################################
			self.snsnpy.EBVgal = 0 #Make sure this metadata is zero, because we've already corrected MWay extinction, so we don't want this modelled in the SED during sn.fit() (otherwise it affects K-corrections)
		elif not snpy_params['apply_EBVMW_corrections']:
			self.snsnpy.EBVgal = 0

		#Initialise SNooPy Product
		snpy_product  = {'kcorrections':{}}
		snpy_product['snpy_cut'] = False
		snpy_product['obsbands'] = copy.deepcopy(self.snsnpy.allbands())
		snpy_product['bandmap']  = copy.deepcopy(self.snsnpy.restbands)

		#Choose observer bands to fit when fitting SNooPy (to get e.g. fitted Tmax, Kcorr from SNooPy fit)
		if snpy_params['interpflts']=='all':	fitbands = [obsband for obsband in snpy_product['bandmap']]
		else:									fitbands = [obsband for obsband in snpy_product['bandmap'] if snpy_product['bandmap'][obsband] in snpy_params['interpflts']]

		#Mapper for obs->rest where rest is in interpflts
		obs_to_rest    = {obsband:self.snsnpy.restbands[obsband] for obsband in fitbands}

		#Get pre-computed dictionary of Kcorrections, MWay extinction corrections, and snpy Tmax estimate
		self.get_fitted_snpy_Tmax_Kcorrdict(fitbands,obs_to_rest)

		#Update snpy_product Tmax_snpy_fitted
		snpy_product['Tmax_snpy_fitted'] = self.Tmax_Kcorrdict['Tmax_snpy_fitted']

		#Apply Milky Way Extinction Corrections
		if snpy_params['apply_EBVMW_corrections']=='Model':#If we are not pre-subtracting MWay Extinction using SED approximation, then use fitted snpy to correct data to EBV_MW=0 by integrating SED over the passband
			for flt in self.Tmax_Kcorrdict['MWCorrections']:
				self.snsnpy.data[flt].mag  += self.Tmax_Kcorrdict['MWCorrections'][flt]

			if self.Tmax_Kcorrdict['snpy_cut_bools']['MWay']:
				print (f"MWay SNooPy Corrections Failed for {self.snsnpy.name}")#; {self.Tmax_Kcorrdict['MWCorrections']}")
				snpy_product['snpy_cut'] = bool(snpy_product['snpy_cut']+self.Tmax_Kcorrdict['snpy_cut_bools']['MWay'])

		#Apply K-corrections
		if snpy_params['apply_K_corrections']:
			####################################
			#Apply Kcorr to observer-frame magnitude data (that has been corrected for MWay Extinction)
			Key = {True:'Mangle',False:'NoMangle'}[snpy_params['mangle']]
			for filt in self.Tmax_Kcorrdict[Key]:
				self.snsnpy.data[filt].mag        += - self.Tmax_Kcorrdict[Key][filt]
				snpy_product['kcorrections'][filt] =   self.Tmax_Kcorrdict[Key][filt]
			####################################
			snpy_product['snpy_cut'] = bool(snpy_product['snpy_cut']+self.Tmax_Kcorrdict['snpy_cut_bools'][Key])

		if snpy_product['snpy_cut']:
			print (f'SNooPy Corrections Failed for {self.snsnpy.name}')

		with open(self.path_snpy_product,'wb') as f:
			pickle.dump(snpy_product,f)

		self.snpy_product = snpy_product


	def __init__(self, choices, outputdir):
		self.choices     = choices
		self.snpychoices = self.choices['snpy_parameters']
		self.outputdir   = outputdir

		self.Tmax_Kcorrdict_folder = self.outputdir+get_Tmax_Kcorrdict_foldername(self.snpychoices)
		self.snana_folder          = self.outputdir+get_snana_foldername(self.snpychoices)

	def get_sn_paths(self):
		"""
		Get SN paths

		Simple Method to set paths up for individual SN

		End Products(s)
		----------
		Paths to snana lc and snpy_product set up as attributes
		"""
		self.path_snana_product  = f"{self.snana_folder}{self.sn}_{self.survey}.snana.dat"
		self.path_snpy_product   = f"{self.snana_folder}{self.sn}_{self.survey}_snpy_product.pkl"

	def load_snpy_product(self, sn=None, survey=None):
		"""
		Load snpy_product

		Simple Method to load snpy_product

		Parameters
		----------
		sn : str (optional; default is None)
			name of sn

		survey : str (optional; default is None)
			name of survey

		Return
		----------
		snpy_product: dict
			Kcorrections, bandmapping, Tmax etc.
		"""
		if sn is not None:
			self.sn     = sn
		if survey is not None:
			self.survey = survey
		self.get_sn_paths()
		with open(self.path_snpy_product,'rb') as f:
			snpy_product = pickle.load(f)
		return snpy_product

	def correct_sn(self, sn, SNsnpy):
		"""
		Correct SN Method

		Sets up pathnames
		Then applies corrections to snsnpy object
		Either computes for the first time, or loads up appropriate MW/Kcorrections
		Final product is snana lc file

		Parameters
		----------
		sn : str
			name of SN
		SNSsnpy: dict
			keys are lc and survey, values are respectively the snsnpy object and the surveyname

		End Product(s)
		----------
		snpy_product: dict
			keys are ['kcorrections', 'snpy_cut', 'obsbands', 'bandmap', 'Tmax_snpy_fitted']

		lc: :py:class:`astropy.table.Table`
			light curve object
		"""
		self.sn       = sn
		self.snsnpy   = SNsnpy['lc']
		self.survey   = SNsnpy['survey']

		#Get paths to snana lc file and snpy_product file
		self.get_sn_paths()

		#Update raw snsnpy to apply Milky Way Extinction and K corrections
		self.correct_data_with_snpy_MWayExtinction_and_K_corrections()

		#Save snsnspy as an snana lc
		self.convert_snsnpy_to_snana()

	def convert_snsnpy_to_snana(self):
		"""
		Convert snsnpy to snana lc

		Method saves/loads snana lc from snsnpy object

		End Product(s)
		----------
		lc: :py:class:`astropy.table.Table`
			light curve object
		"""
		if not os.path.exists(self.path_snana_product) or self.choices['load_data_parameters']['rewrite_snana']:
			#Write lc from snsnpy; we see the filter is always stored as the rest frame filter...and there is question of whether we have include Kcorrections or not in snsnpy.data
			mjd,flt,mag,magerr = [],[],[],[]
			for filt in self.snsnpy.allbands():
				rest_filt = self.snsnpy.restbands[filt]
				#if rest_filt in interpflts and interpflts!='all':#Big Choice here, only allow bands where SNooPy has mapped to rest frame filters in interpflts
				mjd.extend(list(self.snsnpy.data[filt].MJD))
				#Notice here how we designate filter as rest frame, because the data itself has had Kcorrections added on
				flt.extend([rest_filt for _ in range(len(self.snsnpy.data[filt].MJD))])
				mag.extend(list(self.snsnpy.data[filt].mag))
				magerr.extend(list(self.snsnpy.data[filt].e_mag))
			try:
				tmax = snpy_product['Tmax_snpy_fitted']
			except:
				tmax = None
			if mjd==[]:#i.e. if no data falls inside interpflts
				mjd,flt,mag,magerr = [0],['None'],[0],[0]
			snname  = self.snsnpy.name; tmax=tmax; z_helio = self.snsnpy.z; z_cmb=self.snsnpy.get_zcmb(); z_cmb_err=0; ebv_mw = self.snsnpy.EBVgal
			io.write_snana_lcfile(self.snana_folder, snname, mjd, flt, mag, magerr, tmax, z_helio, z_cmb, z_cmb_err, ebv_mw,ra=self.snsnpy.ra,dec=self.snsnpy.decl, survey=self.survey)
