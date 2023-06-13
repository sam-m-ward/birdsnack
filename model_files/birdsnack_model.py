"""
BIRDSNACK Module

Module contains BIRDSNACK class for default analysis: apply corrections, interpolate rest-frame data, analyse with HBM

Contains
--------------------
BIRDSNACK class
    inputs: configname='analysis_config.yaml',loader={},edit_dict = {},dfmeta=None,update_vars=False

    Methods are:
        load_and_preprocess_snana_lcs()

--------------------
"""

import os,pickle,yaml
from snoopy_corrections import SNOOPY_CORRECTIONS
from miscellaneous import ensure_folders_to_file_exist

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
		self.productpath  = self.rootpath+'products/'
		self.SNSpath      = self.productpath+'snpy_SNS/'
		self.snanapath    = self.SNSpath+'snana_copies/'
		for path in [self.analysispath,self.datapath,self.productpath,self.SNSpath,self.snanapath]: ensure_folders_to_file_exist(path)

		#Load up metadata
		self.dfmeta = dfmeta

		#If in preprocessing mode, load up the SNSsnpy dictionary, and load/create the snana lcs
		if self.choices['load_data_parameters']['load_mode']=='preproc':
			path,file,survey = loader['path_file_survey']
			with open(self.SNSpath+file,'rb') as f:
				self.SNSsnpy = pickle.load(f)

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

		#Folder where snana lcs are loaded from and/or saved into
		self.snana_folder = self.snanapath+get_snana_foldername(self.choices['snpy_parameters'])

		for sn in self.SNSsnpy:
			path_snana_product  = f"{self.snana_folder}{sn}_{self.SNSsnpy[sn]['survey']}.snana.dat"
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
			snpy_product = snpycorr.load_snpy_product(sn=sn,survey=self.SNSsnpy[sn]['survey'])
			lc           = update_lcmetadata(lc,self.dfmeta[self.dfmeta['SN']==sn],snpy_product)
			#Append lcs and snpy_products to dictionary
			self.lcs[sn]           = lc
			self.snpy_products[sn] = snpy_product
