import pickle,yaml
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
	configname: str
		name of config.yaml file used for analysis

	edit_dict: dict
		used to edit individual key,value pairs in the config.yaml choices file

	Methods
	----------

	"""

	def __init__(self,loader={},configname='analysis_config.yaml',edit_dict = {},load_dictionary={},dfmeta=None,update_vars=True):
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

		#If loading up snpytxtfiles, load up the SNSsnpy dictionary
		if self.choices['load_data_parameters']['load_mode']=='snpytxtfiles':
			path,file,survey = loader['path_file_survey']
			with open(self.SNSpath+file,'rb') as f:
				self.SNSsnpy = pickle.load(f)

		#Apply snpy corrections to SNSsnpy snpy files if load_mode=='snpytxtfiles' to create corrected snana lcs dictionary
		if self.choices['load_data_parameters']['apply_snpy_corrections']:
			if self.choices['load_data_parameters']['load_mode']!='snpytxtfiles':
				raise Exception("If applying SNooPy corrections via apply_snpy_corrections=True, make sure load_mode is set to 'snpytxtfiles' ")
			else:
				snpycorr = SNOOPY_CORRECTIONS(self.choices, self.snanapath)
				for sn in self.SNSsnpy:
					snpycorr.correct_sn(sn, self.SNSsnpy[sn])
					print (snpycorr.lc)
