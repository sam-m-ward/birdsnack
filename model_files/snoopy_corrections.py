"""
Module containing functions to take in snpy object and return SNR-trimmed rest-frame MWay-Extinction-corrected snana lc object

Contains:
--------------------
correct_lc_data_SNR:
	take in snana lc and cut on SNR

correct_snsnpy_data_SNR:
	take in snpy object and mask on SNR

get_fitted_snpy_Tmax_Kcorrdict:
	computes and/or returns Tmax_Kcorrdict, uniquely defined by snr,EBVMW and SNooPy observer-frame passbands that fall inside interpflts

correct_data_with_snpy_MWayExtinction_and_K_corrections:
	uses SNooPy to correct for MWay extinction and apply K corrections

get_lc_from_snsnpy:
	take in snpy object and create snana lc file

correct_data:
	combines all the above functions, takes in snpy object,
	apply all of SNRcut, MWaydust subtraction, Kcorr (in that order)
	return snana lc object

--------------------
Functions take in various arguments

correct_lc_data_SNR(lc,snrcut,error_boost_dict)

correct_snsnpy_data_SNR(snsnpy,snrcut)

get_fitted_snpy_Tmax_Kcorrdict(snsnpy,snpy_params,fitbands,outputdir)

correct_data_with_snpy_MWayExtinction_and_K_corrections(snsnpy, snpy_params, path_snpy_product, method='spline')

get_lc_from_snsnpy(snsnpy,snrcut,snpy_product,outputdir,rewrite_snana,path,error_boost_dict)

correct_data(sn,snsnpy,snpy_params,outputdir,rewrite_snana,error_boost_dict,filts_pecking_order=None)

--------------------
Functions use simple operations

Written by Sam M. Ward: smw92@cam.ac.uk
"""
import copy
import numpy as np
from snpy import fset
import extinction
import os
import pickle
from miscellaneous import ensure_folders_to_file_exist

def get_Tmax_Kcorrdict_foldername(snpy_params):
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

def get_snana_foldername(snpy_params):
	appender = ''

	for key in snpy_params:
		#If apply_SNR_cut is False, don't appned snrcut
		if key=='snrcut' and not snpy_params['apply_SNR_cut']:
			pass

		#If apply_EBVMW_corrections is not Presubtract, don't append RVMW or dustlaw
		elif key in ['RVMW','dustlaw'] and snpy_params['apply_EBVMW_corrections']!='Presubtract':
			pass

		#If apply_K_corrections is False, don't append mangling choice
		elif key=='mangle' and not snpy_params['apply_K_corrections']:
			pass

		#Otherwise, append the snpy correction choices
		else:
			appender += f"{key}{str(snpy_params[key])}_"


	return appender[:-1]+'/'

class SNOOPY_CORRECTIONS:

	def correct_lc_data_SNR(lc,snrcut,error_boost_dict):
		"""
		Cut on SNR

		Returns lc object by applying a cut on Signal-to-Noise-Ratio

		Parameters
		----------
		lc: :py:class:`astropy.table.Table`
			light curve objects

		snrcut: float or int
			SNR>snrcut is kept

		error_boost_dict: dict
			keys are factors to multiply flux errors by, values are list of SNe which correction should be applied to

		Returns
		----------
		lc object with SNR>snrcut
		"""
		lc = lc[lc["fluxcal"]/lc["fluxcalerr"] > snrcut]
		for fac in error_boost_dict:
			if lc.meta['SNID'] in error_boost_dict[fac]:
				lc['fluxcalerr'] *= fac
		return lc

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
					#Can't really change RV_host, uses calibration number 6 from Folatelli2010, which is RVhost=1.01+/-0.31
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

	def correct_data_with_snpy_MWayExtinction_and_K_corrections(self):#, path_snpy_product, outputdir, filts_pecking_order):
		"""
		Correct Data using SNooPy for Milky Way Extinction and Redshifting

		Returns corrected snpy object, and dictionary of important snpy values e.g. kcorrections, rest-frame filters, Tmaxs

		Returns
		----------
		snsnpy: snpy.sn.sn object
			snpy object of particular sn (corrected for MW-Extinction and Kcorrections)

		snpy_product: dict
			keys are ['kcorrections', 'snpy_cut', 'obsbands', 'bandmap', 'Tmax_GP_obsframe', 'Tmax_snpy_obsframe', 'Tmax_snpy_restframe', 'Tmax_snpy_fitted']
			kcorrections are dict={flt:np.array of kcorrections}
			snpy_cut is bool, if True, SNooPy procedure failed
			obsbands and bandmap are SNooPy filters, (key,value) pair are obs and rest bands
			Finally, Tmax float values from SNooPy fits
		"""
		snpy_params = self.snpychoices
		if snpy_params['snpy_model'] and snpy_params['snpy_shape']:#Default is simply False and False
			print ('Applying .choose_model')
			self.snsnpy.choose_model(snpy_model, stype=snpy_shape)

		snpy_product  = {'kcorrections':{}}
		def return_correct_restbands(snsnpy):
			#Implements x2 Rules to boost sample size
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

		snpy_product['snpy_cut'] = False
		snpy_product['obsbands'] = copy.deepcopy(self.snsnpy.allbands())
		snpy_product['bandmap']  = copy.deepcopy(self.snsnpy.restbands)

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

		#Choose observer bands to fit when fitting SNooPy (to get e.g. fitted Tmax, Kcorr from SNooPy fit)
		if snpy_params['interpflts']=='all':
			fitbands = [obsband for obsband in snpy_product['bandmap']]
		else:
			fitbands = [obsband for obsband in snpy_product['bandmap'] if snpy_product['bandmap'][obsband] in snpy_params['interpflts']]
		obs_to_rest    = {obsband:self.snsnpy.restbands[obsband] for obsband in fitbands}
		self.get_fitted_snpy_Tmax_Kcorrdict(fitbands,obs_to_rest)

		snpy_product['Tmax_snpy_fitted'] = self.Tmax_Kcorrdict['Tmax_snpy_fitted']

		if snpy_params['apply_EBVMW_corrections']=='Model':#If we are not pre-subtracting MWay Extinction using SED approximation, then use fitted snpy to correct data to EBV_MW=0 by integrating SED over the passband
			for flt in self.Tmax_Kcorrdict['MWCorrections']:
				self.snsnpy.data[flt].mag  += self.Tmax_Kcorrdict['MWCorrections'][flt]

			if self.Tmax_Kcorrdict['snpy_cut_bools']['MWay']:
				print (f'MWay SNooPy Corrections Failed for {self.snsnpy.name}')
				snpy_product['snpy_cut'] = bool(snpy_product['snpy_cut']+self.Tmax_Kcorrdict['snpy_cut_bools']['MWay'])

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

	def correct_sn(self, sn, SNsnpy):

		self.sn       = sn
		self.snsnpy   = SNsnpy['lc']
		self.survey   = SNsnpy['survey']

		try:
			self.apply_corrections()
		except Exception as e:
			print (e)
			pass

		self.lc       = {sn:{'lc':None,'survey':self.survey}}

	def apply_corrections(self):
		"""
		Correct snpy object data

		Returns light curve object by correcting snpy object for: SNRcut, Milky Way Dust, and K-corrections (in that order)

		"""

		##############################################################
		self.path_snana_product  = self.snana_folder+self.sn+'.snana.dat';
		self.path_snpy_product   = self.snana_folder+self.sn+'_snpy_product.pkl'

		if self.choices['load_data_parameters']['rewrite_snana'] or not os.path.exists(self.path_snana_product) or not os.path.exists(self.path_snpy_product):
			if self.snpychoices['apply_SNR_cut']: self.snsnpy.mask_SNR(self.snpychoices['snrcut'])
			self.correct_data_with_snpy_MWayExtinction_and_K_corrections()

		snpy_product = self.snpy_product

		err=1/0
		self.get_lc_from_snsnpy(error_boost_dict)
		return lc, snpy_product

	def get_lc_from_snsnpy(snsnpy,snrcut,interpflts,snpy_product,outputdir,rewrite_snana,path,error_boost_dict):
		"""
		Get Light Curve from snpy object

		Function loads up an snana lc given a sn name, and creates this file from a snpy object if it doesn't already exist

		Parameters
		----------
		snsnpy: snpy.sn.sn object
			snpy object of particular sn

		snrcut: float or int
			SNR>snrcut is kept

		interpflts: str
			string of filters we allow SNooPy to use when fitting, e.g. 'BVrRiIYJH'
			in this function, trim on interpflts so only these bands appear in snana lc object
			e.g. if interpflts=='uBgVriYJH' then exclude Bs,Vs,J_K etc. (which are each rest-frame bands already naturally mapped to by SNooPy)

		snpy_product: dict
			keys are ['kcorrections', 'snpy_cut', 'obsbands', 'bandmap', 'Tmax_GP_obsframe', 'Tmax_snpy_obsframe', 'Tmax_snpy_restframe']
			kcorrections are dict={flt:np.array of kcorrections}
			snpy_cut is bool, if True, SNooPy procedure failed
			obsbands and bandmap are SNooPy filters, (key,value) pair are obs and rest bands
			Finally, Tmax float values from SNooPy fits

		outputdir: str
			path/to/dir/where_snana_lcs_are_saved

		rewrite_snana: bool
			if True, create new snan lc object

		path: str
			path/to/snanalc.dat

		error_boost_dict: dict
			keys are factors to multiply flux errors by, values are list of SNe which correction should be applied to

		Returns
		----------
		lc: :py:class:`astropy.table.Table`
			light curve object
		"""
		#path = outputdir+snsnpy.name+'.snana.dat'
		if not os.path.exists(path) or rewrite_snana:
			#################################
			#Write lc from snsnpy; we see the filter is always stored as the rest frame filter...and there is question of whether we have include Kcorrections or not in snsnpy.data
			mjd,flt,mag,magerr = [],[],[],[]
			for filt in snsnpy.allbands():
				rest_filt = snsnpy.restbands[filt]
				#if rest_filt in interpflts and interpflts!='all':#Big Choice here, only allow bands where SNooPy has mapped to rest frame filters in interpflts
				mjd.extend(list(snsnpy.data[filt].MJD))
				#Notice here how we designate filter as rest frame, because the data itself has had Kcorrections added on
				flt.extend([rest_filt for _ in range(len(snsnpy.data[filt].MJD))])
				mag.extend(list(snsnpy.data[filt].mag))
				magerr.extend(list(snsnpy.data[filt].e_mag))
			try:
				tmax = snpy_product['Tmax_snpy_restframe']
			except:
				tmax = None
			if mjd==[]:#i.e. if no data falls inside interpflts
				mjd,flt,mag,magerr = [0],['None'],[0],[0]
			snname  = snsnpy.name; tmax=tmax; z_helio = snsnpy.z; z_cmb=snsnpy.get_zcmb(); z_cmb_err=0; ebv_mw = snsnpy.EBVgal
			io.write_snana_lcfile(outputdir, snname, mjd, flt, mag, magerr, tmax, z_helio, z_cmb, z_cmb_err, ebv_mw,ra=snsnpy.ra,dec=snsnpy.decl)
			#################################
		sn , lc = io.read_snana_lcfile(path)
		lc      = correct_lc_data_SNR(lc,snrcut,error_boost_dict)
		lc      = get_mag(lc)
		return lc
