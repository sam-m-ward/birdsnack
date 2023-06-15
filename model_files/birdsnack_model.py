"""
BIRDSNACK Module

Module contains BIRDSNACK class for default analysis: apply corrections, interpolate rest-frame data, analyse with HBM

Contains
--------------------
BIRDSNACK class
	inputs: configname='analysis_config.yaml',loader={},edit_dict = {},dfmeta=None,update_vars=False

	Methods are:
		load_and_preprocess_snana_lcs()
		trim_sample(apply_trims=True, return_trimmed=False)
		load_empty_DF_M()
		get_peak_mags(savekey='Default', overwrite=False)

--------------------
"""

import copy,os,pickle,yaml
from snoopy_corrections import SNOOPY_CORRECTIONS
from miscellaneous import *
from LC_object import *

class BIRDSNACK:
	"""
	BIRDSNACK Class Object

	Master class object that in snana light curve files
	Performs preprocessing cuts
	Uses GP to get peak apparent colours
	Uses stan to fit for intrinsic/extrinsic components

	Parameters
	----------
	configname: str (optional; default='analysis_config.yaml')
		name of config.yaml file used for analysis

	loader : dict (optional; default={})
		dictionary of information important to analysis
		e.g. loader['path_file_survey'] = [path,file,survey]

	edit_dict: dict (optional; default={})
		used to edit individual key,value pairs in the config.yaml choices file

	dfmeta : pandas df (optional; default = None)
		DataFrame of SN host galaxy stellar masses and spectroscopic subclassifications

	update_vars : bool (optional; default=False)
		if True, reload BIRDSNACK with updated .yaml choices

	Methods
	----------

	"""

	def __init__(self,configname='analysis_config.yaml',loader={},edit_dict = {},dfmeta=None,update_vars=False):
		#Set Configname
		self.configname   = configname

		#Load up config.yaml choices
		with open(self.configname) as f:
			self.choices = yaml.load(f, Loader=yaml.FullLoader)

		#Make edits to config.yaml choices
		for key,value in edit_dict.items():
			self.choices[key] = value

		#Set Pathnames
		self.rootpath     = self.choices['rootpath']
		self.analysispath = self.rootpath+'analysis/'
		self.datapath     = self.rootpath+'data/'
		self.plotpath     = self.rootpath+'plots/'
		self.productpath  = self.rootpath+'products/'
		self.DFpath	      = self.productpath+'DFs/'
		self.SNSpath      = self.productpath+'snpy_SNS/'
		self.snanapath    = self.SNSpath+'snana_copies/'
		for path in [self.analysispath,self.datapath,self.plotpath,self.productpath,self.DFpath,self.SNSpath,self.snanapath]: ensure_folders_to_file_exist(path)

		#Load up metadata
		self.dfmeta = dfmeta

		#If in preprocessing mode, load up the SNSsnpy dictionary, and load/create the snana lcs
		if self.choices['load_data_parameters']['load_mode']=='preproc':
			if 'path_file_survey' in loader:
				path,file,survey = loader['path_file_survey']
				with open(self.SNSpath+file,'rb') as f:
					self.SNSsnpy = pickle.load(f)
				self.survey  = survey
			elif 'SNSsnpy' in loader:
				self.SNSsnpy = loader['SNSsnpy']

			#Load up snana lcs and apply SNR cut, error_boost_dict, and trim to interpflts
			self.load_and_preprocess_snana_lcs()


	def load_and_preprocess_snana_lcs(self):
		"""
		Load and Pre-process snana lcs

		Method takes in list of sns from SNSsnpy, then outputs dictionary of snana lcs={sn:lc}
		Each lc is preprocessed to apply SNRcut, error boosting, create magnitudes column, trim to interpflts, and update metadata with mass and spectype
		If snana lc object with SNooPy corrections applied doesn't exist to start with, it will create this

		End Product(s)
		----------
		lcs : dict
			{sn:lc} where lc is snana lc objects
		"""
		#Initialise SNooPy Corrections
		snpycorr = SNOOPY_CORRECTIONS(self.choices, self.snanapath)

		#Import snana_functions
		from snana_functions import get_snana_foldername, correct_lc_data_SNR_boostdict, get_mag, trim_filters, set_lcmeta_ordered_flts, update_lcmetadata
		import snanaio as io

		#Initialise products
		self.lcs = {} ; self.snpy_products = {}
		self.sns = list(self.SNSsnpy.keys())

		#Folder where snana lcs are loaded from and/or saved into
		self.snana_folder = self.snanapath+get_snana_foldername(self.choices['snpy_parameters'])

		for sn in self.SNSsnpy:
			survey = self.SNSsnpy[sn]['survey']
			path_snana_product  = f"{self.snana_folder}{sn}_{survey}.snana.dat"
			if not os.path.exists(path_snana_product) or self.choices['load_data_parameters']['rewrite_snana']:
				print ('###'*3)
				snpycorr.correct_sn(sn, self.SNSsnpy[sn])

			#Load File
			sn , lc      = io.read_snana_lcfile(path_snana_product)
			#Apply SNR cut, boost flux errors, get magnitudes column
			lc           = get_mag(correct_lc_data_SNR_boostdict(lc, self.choices['snpy_parameters']['snrcut'], self.choices['preproc_parameters']['error_boost_dict']))
			#Trim to filters in interpflts
			lc           = set_lcmeta_ordered_flts(trim_filters(lc, self.choices['snpy_parameters']['interpflts']))
			#Update lcmeta with mass, spectral type, and SNooPy Tmax estimate
			snpy_product = snpycorr.load_snpy_product(sn=sn,survey=survey)
			lc           = update_lcmetadata(lc,self.dfmeta[self.dfmeta['SN']==sn],snpy_product)
			#Append lcs and snpy_products to dictionary
			self.lcs[sn]           = lc
			self.snpy_products[sn] = snpy_product

	def trim_sample(self, apply_trims=True, return_trimmed=False):
		"""
		Trim Sample

		Work through self.lcs dictionary and trim based on: SNooPy cut, GPTmax error, availability of data.

		Parameters
		----------
		apply_trims : bool (optional; default=True)
			if True, trim self.lcs

		return_trimmed : bool (optional; default=False)
			if True, return trimmed_lcs and the reasons dictionary

		Returns
		----------
		if return_trimmed, returns:
			trimmed_lcs : dict
				{sn:lc} of lcs trimmed
			REASONS : dict
				{reason:[list of sne]} reasons for being cut from sample

		End Product(s)
		----------
		if apply_trims:
			update self.lcs and self.sns; excludes trimmed SNe
		"""

		def trims_on_SNooPy_fit_and_Tmax(reasons):
			"""
			Trims on SNooPy Fit and Tmax

			Uses lc.meta data to trim on whether SNooPy fit succeeded/failed, and whether Tmax estimate succeeded/failed

			Parameters
			----------
			reasons : dict
				reasons[reason] = [list of SNe]

			Returns
			----------
			trim : bool
				if True, SN was trimmed

			updated reasons dict
			"""
			#Initialise trim bool and get snpy_products
			trim = False ; snpy_product = self.snpy_products[sn]
			use_SNooPy_fitted_Tmax = trim_choices['Tmaxchoicestr']=='Tmax_snpy_fitted'
			#If SNooPy corrections failed, record SN trim
			if snpy_product['snpy_cut']:
				reasons['snpy_cut'].append(sn) ; trim = True
			#If SNooPy Tmax estimate failed and using SNooPy Tmax, record SN trim
			if use_SNooPy_fitted_Tmax:
				if lc.meta['Tmax_snpy_fitted'] is None:
					reasons['None_Tmax_snpy_fitted'].append(sn) ; trim = True
			#If no GP Tmax estimate, or its error is too large, record SN trim
			elif not use_SNooPy_fitted_Tmax:
				if lc.meta['Tmax_GP_restframe'] is None:
					reasons['None_Tmax_GP_restframe'].append(sn) ; trim = True

				if lc.meta['Tmax_GP_restframe'] is not None:
					if lc.meta['Tmax_GP_restframe_std']>Tmax_stdmax:
						reasons['Large_Tmax_GP_restframe_std'].append(sn) ; trim = True

			return trim, reasons

		trim_choices = self.choices['preproc_parameters']
		pblist = trim_choices['pblist'] ; tilist = trim_choices['tilist'] ; Extra_Features = trim_choices['Extra_Features'] ; Tmax_stdmax = trim_choices['Tmax_stdmax']

		#Initialise reasons dictionary for trimming each SN
		#SNooPy fit failed, bad Tmax estimate, or not enough points around reference band
		reasons  = {'snpy_cut':[],'None_Tmax_GP_restframe':[],'None_Tmax_snpy_fitted':[],'Large_Tmax_GP_restframe_std':[],'Points_Before_After_RefBand_Max':[]}
		#Not enough points around non-reference passbands
		reasons2 = {f'{pb}_{ti}':[] for pb in pblist for ti in tilist}
		#Not enough points around special bands/phases, e.g. B-band 15days
		reasons3 = {f'Extra_{pb}_{ti}':[] for ti,pbmini in Extra_Features.items() for pb in pbmini}


		new_lcs = {} ; trimmed_lcs = {}

		if trim_choices['trim_on_refband_max']:
			Nlist = [trim_choices[x] for x in ['Nbeforemax','Nbeforegap','Naftermax','Naftergap']]
			Nbeforemax,Nbeforegap,Naftermax,Naftergap = Nlist[:]
			print ('############################')
			print (f"Trimming lcs so in {trim_choices['ref_band']} there {'are' if Nbeforemax!=1 else 'is'}: \n{Nbeforemax} point{'s' if Nbeforemax!=1 else ''} within {Nbeforegap} days before max, and \n{Naftermax} point{'s' if Naftermax!=1 else ''} within {Naftergap} days after max \n...")
			for sn,lc in self.lcs.items():
				lcobj = LCObj(lc,trim_choices)
				lcobj.get_Tmax()
				#With GP Tmax estimated, can perform generic trims
				trim_bool, reasons = trims_on_SNooPy_fit_and_Tmax(reasons)
				#Check to see Points_Before_After_RefBand_Max
				if lc.meta[trim_choices['Tmaxchoicestr']] is not None:#If Tmax is not estimated, can't compute phase, therefore can't assess availability of points
					trim = lcobj.test_data_availability(ref_band = trim_choices['ref_band'], tref = trim_choices['tilist'][trim_choices['tref_index']], Nlist = Nlist)
					if trim: reasons['Points_Before_After_RefBand_Max'].append(sn)
					trim_bool += trim
				if trim_bool: trimmed_lcs[sn] = lc

		if trim_choices['trim_on_pblist']:
			Nlist = [trim_choices[x] for x in ['Nbefore_local','Nbeforegap_local','Nafter_local','Naftergap_local']]
			Nbefore_local,Nbeforegap_local,Nafter_local,Naftergap_local = Nlist[:]
			print ('############################')
			print (f"Trimming lcs so in each of {pblist} there {'are' if Nbefore_local!=1 else 'is'}:\n{Nbefore_local} point{'s' if Nbefore_local!=1 else ''} within {Nbeforegap_local} days before {tilist}, and\n{Nafter_local} point{'s' if Nafter_local!=1 else ''} within {Naftergap_local} days after {tilist} \n...")

			for sn,lc in self.lcs.items():
				lcobj = LCObj(lc,trim_choices)
				lcobj.get_Tmax()
				#With GP Tmax estimated, can perform generic trims
				trim_bool, reasons = trims_on_SNooPy_fit_and_Tmax(reasons)
				#Check to see pbtilist points
				if lc.meta[trim_choices['Tmaxchoicestr']] is not None:#If Tmax is not estimated, can't compute phase, therefore can't assess availability of points
					for pb in pblist:
						for ti in tilist:
							trim = lcobj.test_data_availability(ref_band = pb, tref = ti, Nlist = Nlist)
							if trim: reasons2[f'{pb}_{ti}'].append(sn)
							trim_bool += trim
				if trim_bool: trimmed_lcs[sn] = lc

		if trim_choices['trim_on_extras']:
			Nlist = [trim_choices[x] for x in ['Nbefore_local_extra','Nbeforegap_local_extra','Nafter_local_extra','Naftergap_local_extra']]
			Nbefore_local_extra,Nbeforegap_local_extra,Nafter_local_extra,Naftergap_local_extra = Nlist[:]
			print ('############################')
			print (f"Trimming lcs using {Extra_Features} there {'are' if Nbefore_local_extra!=1 else 'is'}: \n{Nbefore_local_extra} point{'s' if Nbefore_local_extra!=1 else ''} within {Nbeforegap_local_extra} days before phases, and \n{Nafter_local_extra} point{'s' if Nafter_local_extra!=1 else ''} within {Naftergap_local_extra} days after phases \n...")
			for sn,lc in self.lcs.items():
				lcobj = LCObj(lc,trim_choices)
				lcobj.get_Tmax()
				#With GP Tmax estimated, can perform generic trims
				trim_bool, reasons = trims_on_SNooPy_fit_and_Tmax(reasons)
				#Check to see Extras points
				if lc.meta[trim_choices['Tmaxchoicestr']] is not None:
					for ti,pbmini in Extra_Features.items():
						for pb in pbmini:
							trim = lcobj.test_data_availability(ref_band = pb, tref = ti, Nlist = Nlist, local=True)
							if trim: reasons3[f'Extra_{pb}_{ti}'].append(sn)
							trim_bool += trim
				if trim_bool: trimmed_lcs[sn] = lc

		#Bring all reasons together
		REASONS  = {**reasons, **reasons2, **reasons3}

		#Create trimmed_lcs and complement==new_lcs
		trimmed_lcs = {sn:trimmed_lcs[sn] for sn in self.sns if sn     in list(trimmed_lcs.keys())} #Correct the order
		new_lcs     = {sn:self.lcs[sn]    for sn in self.sns if sn not in list(trimmed_lcs.keys())}


		print (f'(Original LCs : {len(self.lcs)})')
		print (f' Retained LCs : {len(new_lcs)}')
		print (f'  Trimmed LCs : {len(trimmed_lcs)}')
		for sn in trimmed_lcs:
			sreasons = [reason for reason in REASONS if sn in REASONS[reason]]
			print (f'{sn} trimmed for: {sreasons}')
		print ('############################')
		self.reasons     = REASONS
		self.trimmed_lcs = trimmed_lcs
		if apply_trims:
			self.total_sns = copy.deepcopy(list(self.lcs.keys()))
			self.lcs = new_lcs
			self.sns = [sn for sn in list(self.lcs.keys())]
		if return_trimmed:
			return trimmed_lcs, REASONS

	def load_empty_DF_M(self):
		"""
		Load Empty DF

		Method to create an empty DF that can then be filled up

		End Product(s)
		----------
		DF: dict of pandas.df
			DF = {ti:df for ti in tilist}; df is GP interpolation data at a given phase, e.g. apparent mag, colours etc.
			Also includes {'extra':DF_extra} for e.g. dm15B
		"""
		pblist = self.choices['preproc_parameters']['pblist']
		tilist = self.choices['preproc_parameters']['tilist']
		errstr = self.choices['preproc_parameters']['errstr']

		#df for specific phase across all bands in pblist
		DF_m_t = pd.DataFrame(columns = pblist+[pb+errstr for pb in pblist])
		#dictionary for all phases in tilist
		DF_m   = {ti:copy.deepcopy(DF_m_t) for ti in tilist}

		#Additional df for special bands/phases
		DF_extra = {}
		for ti,pbmini in self.choices['preproc_parameters']['Extra_Features'].items():
			for pb in pbmini:
				DF_extra[ti] = pd.DataFrame(columns = pbmini+[pb+errstr for pb in pbmini])

		#Combine
		DF_M = {**DF_m,**{'extra':DF_extra}}
		return DF_M

	def get_peak_mags(self, savekey='Default', overwrite=False):
		"""
		Get Peak Magnitudes

		Interpolate rest-frame data to reference time(s) and extract pandas df measurements and errors

		Parameters
		----------
		savekey : str (optional; default='Default')
			name for DF_M file

		overwrite : bool (optional; default=False)
			if True, recompute DF_M, regardless of whether or not file exists

		Returns
		----------
		DF_M: dict of pandas.df
			DF = {ti:df_ti for ti in tilist,**{'extra':df_ti_extra}}; df_ti is GP interpolation data at a given phase for pblist, df_ti_extra is same for Extra_Features
		"""
		DF_M_name = f'{self.DFpath}DF_M_{savekey}.pkl'

		if overwrite or not os.path.exists(DF_M_name):
			print (f'Creating DF_M for lcs(N={len(self.lcs)})')
			#Initialise empty sample-level dataframe of interpolations
			DF_M = self.load_empty_DF_M()
			for sn,lc in self.lcs.items():
				#Initialise LCobj of lc
				lcobj = LCObj(lc,self.choices['preproc_parameters'])
				#Get Tmax
				lcobj.get_Tmax()
				#Trim on phase
				lcobj.get_phase()
				#Interpolate data
				lcobj.get_GP_interpolation()
				#Extract Interpolations at times and passbands of interest
				lcobj.extract_interpolations()
				#Update sample-level DF with sn-specific interpolations
				DF_M = lcobj.update_DF_M(DF_M, sn)
				#Plot light curve interpolation
				#lcobj.plot_lc(PLOTTER(self.choices['plotting_parameters'],self.plotpath))
		self.DF_M = DF_M

	def additional_cuts(self):
		"""
		Additional Cuts

		Methods to apply additional cuts sample, based on metadata, e.g. B-V colour, Host galaxys stellar mass, Spectroscopic type etc., and also measurement errors

		End Product(s)
		----------
		DF_M with appropriate cuts applied
		"""
		DF_M    = self.DF_M ; tref = self.choices['preproc_parameters']['tilist'][self.choices['preproc_parameters']['tref_index']]
		choices = self.choices['additional_cut_parameters']

		tilist = self.choices['preproc_parameters']['tilist']
		pblist = self.choices['preproc_parameters']['pblist']
		Extra_Features = self.choices['preproc_parameters']['Extra_Features']

		print (f"Original sample size:{DF_M[tref].shape[0]}")

		#Remove NaNs
		print ("###"*5+"\nRemoving NaNs")
		for ti in tilist:
			DF_M[ti].dropna(inplace=True)
		for ti in Extra_Features:
			DF_M['extra'][ti].dropna(inplace=True)
		DF_M = trim_to_common_SNS(DF_M)


		#Remove spectrosopically peculiar SNe
		if choices['spectype']=='normal':
			print ("###"*5+f"\nCutting to spectype={choices['spectype']}")
			drop_sns = []
			for sn in list(DF_M[tref].index):
				if self.lcs[sn].meta['Spectral_Type']!='normal':
					if self.lcs[sn].meta['Spectral_Type'] is None:
						print (f'Dropping SN{sn} because no spec data')
					else:
						print (f"Dropping SN{sn} because spec data is {self.lcs[sn].meta['Spectral_Type']}")
					drop_sns.append(sn)
			DF_M = trim_to_common_SNS(DF_M, drop_sns=drop_sns)

		#Trim on magnitude measurement errors
		if choices['magerrcut']:
			print ("###"*5+f"\nCutting so magnitude measurement errors are <{choices['magerrcut']}")
			for ti in tilist:
				DF_M[ti] = DF_M[ti][DF_M[ti][[col for col in DF_M[ti].columns if self.choices['preproc_parameters']['errstr'] in col]]<choices['magerrcut']]
			for ti in Extra_Features:
				DF_M['extra'][ti] =  DF_M['extra'][ti][DF_M['extra'][ti][[col for col in DF_M['extra'][ti].columns if self.choices['preproc_parameters']['errstr'] in col]]<choices['magerrcut']]
			DF_M = trim_to_common_SNS(DF_M)

		#Trim pre-selected SNe and quote reason why
		if choices['extra_drop_SNe']!={}:
			print ("###"*5+"\nRemoving pre-selected SNe")
			for sn,reason in choices['extra_drop_SNe'].items():
				with suppress(KeyError):
					DF_M[tref].drop([sn], axis=0, inplace=True)
					print (reason)
			DF_M = trim_to_common_SNS(DF_M)

		#Apply B-V<0.3 mag cosmology cut
		if choices['BVcut']:
			print ("###"*5+"\nApplying B-V<0.3 mag cut")
			DF_M[0] = DF_M[0][(DF_M[0]['B']-DF_M[0]['V']).abs()<0.3]
			DF_M = trim_to_common_SNS(DF_M)

		#Apply host galaxy stellar mass cut
		if choices['mass_mode']!='all_masses':
			print ("###"*5+f"\nCutting on host galaxy stellar mass to retain: {choices['mass_mode']}; defined at logMcut={choices['logMcut']}")
			drop_sns = []
			for sn in list(DF_M[tref].index):
				Mbest = self.lcs[sn].meta['Mbest']
				if Mbest is None:
					drop_sns.append(sn)
				elif Mbest>=logMcut and choices['mass_mode']=='low_masses':
					drop_sns.append(sn)
				elif Mbest<logMcut and choices['mass_mode']=='high_masses':
					drop_sns.append(sn)
			DF_M = trim_to_common_SNS(DF_M, drop_sns=drop_sns)

		#Reset attribute with cuts to sample
		self.DF_M = DF_M
