"""
BIRDSNACK Module

Module contains BIRDSNACK class for default analysis: apply corrections, interpolate rest-frame data, analyse with HBM

Contains
--------------------
BIRDSNACK class
	inputs: configname='analysis_config.yaml',loader={},edit_dict = {},dfmeta=None,update_vars=False

	Methods are:
		load_and_preprocess_snana_lcs()
		inspect_data_availability(plot=False, get_outside_range=True, return_others=False)
		trim_sample(apply_trims=True, return_trimmed=False)
		load_empty_DF_M()
		get_peak_mags(savekey='Default', overwrite=False)
		plot_lcs(sns=None,**kwargs)
		additional_cuts(cutter=True)
		fit_stan_model(save=True,Rhat_threshold=1.02,Ntrials=3,additive=False,adder=1000)
		plot_posterior_samples()
		plot_mag_deviations()
		plot_colour_corner(use_intrinsic_mean=False)
		plot_int_model_posteriors(PAR='muRV',FITS=None)

Functions include:
	get_edit_dict(choices,CYCLE_DICT,HBM_savekey)
	map_fit_summary(fitsummary)
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

import arviz as az
import copy,os,pickle,yaml
from HBM_preprocessing import *
from LC_object import *
from miscellaneous import *
from snoopy_corrections import *
from plotting_functions import get_parlabels, get_Lines
from posterior_plotter import *
import stan

class BIRDSNACK:
	"""
	BIRDSNACK Class Object

	Master class object that takes in snana light curve files
	Performs preprocessing cuts
	Uses GP to get peak apparent colours
	Uses stan to fit for intrinsic/extrinsic components
	"""

	def __init__(self,configname='analysis_config.yaml',loader={},edit_dict = {},dfmeta=None,update_vars=False):
		"""
		Initialisation

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
		"""
		#Set Configname
		self.configname   = configname

		#Load up config.yaml choices
		with open(self.configname) as f:
			self.choices = yaml.load(f, Loader=yaml.FullLoader)

		#Make edits to config.yaml choices
		for glob_key in edit_dict:
			if glob_key=='rootpath':
				self.choices[glob_key] = edit_dict[glob_key]
			else:
				for key,value in edit_dict[glob_key].items():
					if key in self.choices[glob_key]:#If this choice is already in .yaml
						if type(self.choices[glob_key][key])==dict:
							self.choices[glob_key][key] = {**self.choices[glob_key][key],**value}
						else:
							self.choices[glob_key][key] = value
					else:#If this is a new entry
						self.choices[glob_key][key] = value

		#Set Pathnames
		self.rootpath     = self.choices['rootpath']
		self.datapath     = self.rootpath+'data/'
		self.plotpath     = self.rootpath+'plots/'
		self.productpath  = self.rootpath+'products/'
		self.modelpath    = self.rootpath+'model_files/'
		self.DFpath	      = self.productpath+'DFs/'
		self.SNSpath      = self.productpath+'snpy_SNS/'
		self.FITSpath     = self.productpath+'stan_fits/FITS/'
		self.snanapath    = self.SNSpath+'snana_copies/'
		for path in [self.datapath,self.plotpath,self.productpath,self.modelpath,self.DFpath,self.SNSpath,self.FITSpath,self.snanapath]:
			ensure_folders_to_file_exist(path)

		#Load up metadata
		self.dfmeta = dfmeta

		#Load up the SNSsnpy dictionary, and load/create the snana lcs
		if 'path_file_survey' in loader:
			path,file,survey = loader['path_file_survey']
			with open(self.SNSpath+file,'rb') as f:
				self.SNSsnpy = pickle.load(f)
			self.survey  = survey
		elif 'SNSsnpy' in loader:
			self.SNSsnpy = loader['SNSsnpy']

		if 'SNSsnpy' in list(self.__dict__.keys()):
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

	def inspect_data_availability(self, plot=False, get_outside_range=True, return_others=False):
		"""
		Inspect Data Availability

		Method to identify phase of data point closest to peak time in each band

		Parameters
		----------
		plot : bool (optional; default=False)
			if True, plot histogram of data availability

		get_outside_range : bool (optional; default=True)
			if True, return dictionary of {sn:t_pb} which should be cut as a result of not having data point within phase window

		return_others : bool (optional; default=False)
			if True, return other dicionaries with more details

		Returns
		----------
		POINTS : dict
			{sn:{pb:t_pb for pb in pblist}}

		Npbs : dict
			{pb:{t_pb:sn for sne in sns}}
		"""
		POINTS = {} ; trim_choices = self.choices['preproc_parameters']
		pblist = trim_choices['pblist'] ; tilist = trim_choices['tilist']
		#For each sn
		for sn,lc in self.lcs.items():
			lcobj = LCObj(lc,trim_choices)
			try:	lcobj.lc.meta[trim_choices['Tmaxchoicestr']] = self.DF_M['Tmax'][trim_choices['Tmaxchoicestr']].loc[sn]
			except:	lcobj.get_Tmax()
			points = {}
			#Find phase of data point closest to peak time
			if lc.meta[trim_choices['Tmaxchoicestr']] is not None:#If Tmax is not estimated, can't compute phase, therefore can't assess availability of points
				for pb in pblist:
					for ti in tilist:
						t_pb = lcobj.get_closest_data_point(ref_band = pb, tref = ti)
						points[f"{pb}_{ti}"] = t_pb
			#Store dictionary
			POINTS[sn] = points

		#Transform Dictionary for Plotting
		Npbs = {}
		for sn in POINTS:
			for pb in POINTS[sn]:
				if pb in Npbs.keys():
					Npbs[pb] = {**Npbs[pb],**{POINTS[sn][pb]:sn}}
				else:
					Npbs[pb] = {POINTS[sn][pb]:sn}
		if plot:#Plot histograms of closest phases of data to tref
			print ('Plotting Data Availability; Set up for N=6 passbands')
			phase_max = self.choices['additional_cut_parameters']['phase_max'] ; FS=14
			fig,ax = pl.subplots(3,2,figsize=(8,10))
			for ipb,pb in enumerate(Npbs):
			    iy,ix = int(ipb//3),int(ipb%3)
			    ts    = list(Npbs[pb].keys())
			    ax[ix,iy].hist(ts,label=pb[0],color=f'C{ipb}')
			    ax[ix,iy].annotate(f"Median:{np.median(ts):.2}\nStd:{np.std(ts):.2}\nMin:{min(ts):.2}\nMax:{max(ts):.2}",xy=(0.05,0.35),fontsize=FS,xycoords='axes fraction',ha='left')
			    ax[ix,iy].legend(fontsize=FS)
			    ax[ix,iy].tick_params(labelsize=FS)
			fig.subplots_adjust(wspace=0.15,hspace=0.15)
			fig.text(0.5, 0.05, 'Phase of Closest Data Point (days)', rotation='horizontal',fontsize=FS,ha='center')
			fig.text(0.5, 0.9, f'Data Near Peak in Fiducial Sample of {int(len(POINTS))} SNe', rotation='horizontal',fontsize=FS,ha='center')
			pl.savefig(self.plotpath+'Ncuts.pdf',bbox_inches='tight')
			#pl.show()

		if get_outside_range:#Return SNe outside phase window
			outside_range = {} ; phase_max = self.choices['additional_cut_parameters']['phase_max']
			for sn in POINTS:
			    for pb in POINTS[sn]:
			        if np.abs(POINTS[sn][pb])>phase_max:
			            if sn not in list(outside_range.keys()):
			                outside_range[sn] = {pb:POINTS[sn][pb]}
			            else:
			                outside_range[sn] = {**outside_range[sn],**{pb:POINTS[sn][pb]}}
			print (f'Remove SNe given phase window=={phase_max}:',list(outside_range.keys()))
			if return_others:
				return POINTS, Npbs, outside_range
			else:
				return outside_range
		else:
			if return_others:
				return POINTS, Npbs


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

		def trims_on_SNooPy_fit_Tmax_and_data_availability(reasons):
			"""
			Trims on SNooPy Fit, Tmax and Data Availability

			Uses lc.meta data to trim on whether SNooPy fit succeeded/failed,
			whether Tmax estimate succeeded/failed,
			and if there is at least 1 datum per lc in pblist

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
			#Ensure at least one data point somewhere in light curve for each of passbands in pblist
			for pb in pblist:
				if pb not in lc.meta['flts']:
					reasons['Has_pblist'].append(sn) ; trim = True
					break

			return trim, reasons

		trim_choices = self.choices['preproc_parameters']
		pblist = trim_choices['pblist'] ; tilist = trim_choices['tilist'] ; Extra_Features = trim_choices['Extra_Features'] ; Tmax_stdmax = trim_choices['Tmax_stdmax']

		#Initialise reasons dictionary for trimming each SN
		#SNooPy fit failed, bad Tmax estimate, has at least one datum per passband in pblist, or not enough points around reference band
		reasons  = {'snpy_cut':[],'None_Tmax_GP_restframe':[],'None_Tmax_snpy_fitted':[],'Large_Tmax_GP_restframe_std':[],'Points_Before_After_RefBand_Max':[],'Has_pblist':[]}
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
				trim_bool, reasons = trims_on_SNooPy_fit_Tmax_and_data_availability(reasons)
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
				trim_bool, reasons = trims_on_SNooPy_fit_Tmax_and_data_availability(reasons)
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
				trim_bool, reasons = trims_on_SNooPy_fit_Tmax_and_data_availability(reasons)
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

		#Save Tmax estimates
		DF_Tmax = pd.DataFrame(columns=['Tmax_GP_restframe','Tmax_GP_restframe_std','Tmax_snpy_fitted'])

		#Save pre-processing choices made to create the peak mags list
		choices_dict = {glob_key:self.choices[glob_key] for glob_key in self.choices if glob_key in ['snpy_parameters','preproc_parameters']}

		#Combine
		DF_M = {**DF_m,**{'extra':DF_extra},**{'Tmax':DF_Tmax},**{'choices':choices_dict}}
		return DF_M

	def get_peak_mags(self, savekey=None, overwrite=False):
		"""
		Get Peak Magnitudes

		Interpolate rest-frame data to reference time(s) and extract pandas df measurements and errors

		Parameters
		----------
		savekey : str (optional; default is None)
			name for DF_M file

		overwrite : bool (optional; default=False)
			if True, recompute DF_M, regardless of whether or not file exists

		Returns
		----------
		DF_M: dict of pandas.df
			DF = {ti:df_ti for ti in tilist,**{'extra':df_ti_extra}}; df_ti is GP interpolation data at a given phase for pblist, df_ti_extra is same for Extra_Features
		"""
		if savekey is None: savekey = self.choices['preproc_parameters']['DF_savekey']
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
			#Save DF_M
			with open(DF_M_name,'wb') as f:
				pickle.dump(DF_M, f)
		else:#Load DF_M
			print (f'Loading pre-created DF_M=={DF_M_name}')
			with open(DF_M_name,'rb') as f:
				DF_M = pickle.load(f)

			######################
			#Check to ensure DF_M that is being pre-loaded actually agrees with the current pre-processing choices
			#Designed to catch instances where .yaml choices have been changed but DF_savekey is set to load up DF_M from different .yaml choices
			disagreements   = {}
			for glob_key in DF_M['choices']:
				for key,previous_choice in DF_M['choices'][glob_key].items():
					if previous_choice!=self.choices[glob_key][key]:
						#If differeing choices is pblist, and pblist is a subset of original, then all is OK
						if glob_key=='preproc_parameters' and key=='pblist' and len([pb for pb in self.choices[glob_key][key] if pb in previous_choice])==len(self.choices[glob_key][key]):
							pass
						else:#Otherwise, there is a disgreement
							disagreements[key] = [self.choices[glob_key][key],previous_choice]
			if disagreements!={}:
				print (('###!!!'*3+'###')*5)
				for choice in disagreements:
					print (f"Current '{choice}' is '{disagreements[choice][0]}' but DF_savekey='{savekey}' had '{choice}' set to '{disagreements[choice][1]}'")
					print (('###!!!'*3+'###')*5)
				raise Exception('As precaution, please ensure DF_savekey is chosen to align with current .yaml choices')
			######################

		self.DF_M = DF_M

	def plot_lcs(self,sns=None,**kwargs):
		"""
		Plot Light Curves

		Parameters
		----------
		sns : list of str (optional; default=None)
			if list, plot only sns in the list

		**kwargs : any
			optional to change e.g. plotting_parameters

		End Product(s)
		----------
		Saved and/or Displayed plot of light curve
		"""
		#Update kwargs e.g. choice to show plot, save plot etc.
		choices = copy.deepcopy(self.choices)
		for kkey,value in kwargs.items():
			for glob_key in choices:
				for key in choices[glob_key]:
					if key==kkey:
						choices[glob_key][key] = value
		#Choice to include a line at t=0
		include_peak_time_line = bool('include_peak_time_line' in kwargs.keys() and kwargs['include_peak_time_line'])

		#List of sns to plot
		if sns is None:
			sns = list(self.lcs.keys())

		for sn in sns:
			#Initialise LCobj of lc
			lcobj = LCObj(self.lcs[sn],choices['preproc_parameters'])
			#Get Tmax
			lcobj.get_Tmax()
			#Trim on phase
			lcobj.get_phase()
			#Interpolate data
			lcobj.get_GP_interpolation()
			#Plot light curve interpolation
			lcobj.plot_lc(PLOTTER(choices['plotting_parameters'],self.plotpath),include_peak_time_line=include_peak_time_line)


	def additional_cuts(self, cutter=True):
		"""
		Additional Cuts

		Methods to apply additional cuts sample, based on metadata, e.g. B-V colour, Host galaxys stellar mass, Spectroscopic type etc., and also measurement errors

		Parameters
		----------
		cutter : bool (optional; default=True)
			if False, don't cut on these items
			(For example, Kcorrections without mangling leads to large measurement errors, but we want to fit same sample, so apply cutter bool to this cut)

		End Product(s)
		----------
		DF_M with appropriate cuts applied
		"""
		DF_M    = self.DF_M ; tref = self.choices['preproc_parameters']['tilist'][self.choices['preproc_parameters']['tref_index']]
		choices = self.choices['additional_cut_parameters']

		tilist = self.choices['preproc_parameters']['tilist']
		pblist = self.choices['preproc_parameters']['pblist']
		Extra_Features = self.choices['preproc_parameters']['Extra_Features']
		print ('###'*5)
		print (f"Original sample size:{DF_M[tref].shape[0]}")

		#Remove NaNs
		print ("###"*5+"\nRemoving NaNs")
		if self.choices['preproc_parameters']['trim_on_pblist']:
			for ti in tilist:
				DF_M[ti].dropna(inplace=True)
		if self.choices['preproc_parameters']['trim_on_extras']:
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
		if choices['magerrcut'] and cutter:
			print ("###"*5+f"\nCutting so magnitude measurement errors are <{choices['magerrcut']}")
			if self.choices['preproc_parameters']['trim_on_pblist']:
				for ti in tilist:
					errcols = [col for col in DF_M[ti].columns if self.choices['preproc_parameters']['errstr'] in col]
					for col in errcols:
						DF_M[ti] = DF_M[ti][DF_M[ti][col]<choices['magerrcut']]
			if self.choices['preproc_parameters']['trim_on_extras']:
				for ti in Extra_Features:
					errcols = [col for col in DF_M['extra'][ti].columns if self.choices['preproc_parameters']['errstr'] in col]
					for col in errcols:
						DF_M['extra'][ti] = DF_M['extra'][ti][DF_M['extra'][ti][col]<choices['magerrcut']]
			DF_M = trim_to_common_SNS(DF_M)

		#Trim based on phase window
		if choices['phase_max'] is not False:
			outside_range = self.inspect_data_availability()
			print (f"Removing {list(outside_range.keys())} for not having data within phase window {choices['phase_max']}days")
			DF_M = trim_to_common_SNS(DF_M, drop_sns=list(outside_range.keys()))

		#Trim pre-selected SNe and quote reason why
		if choices['extra_drop_SNe']!={}:
			print ("###"*5+"\nRemoving pre-selected SNe")
			for sn,reason in choices['extra_drop_SNe'].items():
				with suppress(KeyError):
					DF_M[tref].drop([sn], axis=0, inplace=True)
					print (f"{sn}:{reason}")
			DF_M = trim_to_common_SNS(DF_M)

		#Apply B-V<0.3 mag cosmology cut
		if choices['BVcut']:
			print ("###"*5+"\n"+r"Applying B-V<%s mag cut"%choices['BVcutval'])
			DF_M[0] = DF_M[0][(DF_M[0]['B']-DF_M[0]['V']).abs()<choices['BVcutval']]
			DF_M = trim_to_common_SNS(DF_M)

		#Apply host galaxy stellar mass cut
		if choices['mass_mode']!='all_masses':
			print ("###"*5+f"\nCutting on host galaxy stellar mass to retain: {choices['mass_mode']}; defined at logMcut={choices['logMcut']}")
			drop_sns = [] ; logMcut = self.choices['additional_cut_parameters']['logMcut']
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
		self.lcs  = {sn:self.lcs[sn] for sn in self.lcs if sn in list(DF_M[list(DF_M.keys())[0]].index)}
		self.sns  = list(self.lcs.keys())
		print ('###'*5)

	def fit_stan_model(self,save=True,Rhat_threshold=1.02,Ntrials=3,additive=False,adder=1000):
		"""
		Fit Stan Model

		Method to perform hierarchical Bayesian inference on magnitude measurements to infer pop. dists. in intrinsic and extrinsic components

		Parameters
		----------
		save : bool (optional; default=True)
			if save, save the posterior samples

		Rhat_threshold : float (optional; default=1.02)
			default value of Rhat for Rhat_check_params which cannot be exceeded else fit will be done again with more samples

		Ntrials : float (optional; default=3)
			keep doing fit, multiplying up n_sampling, until Rhats<threshold

		additive : bool (optional; default=False)
			if True, add on to n_sampling, rather than multiplying up, for repeat trials when Rhat too large

		adder : float (optional; default=1000)
			number of samples to add on to each chain if additive=True

		Returns
		----------
		"""

		DF_M    = self.DF_M
		savekey = self.choices['analysis_parameters']['HBM_savekey']
		pblist  = self.choices['preproc_parameters']['pblist']
		tref    = self.choices['preproc_parameters']['tilist'][self.choices['preproc_parameters']['tref_index']]
		DataTransformation = self.choices['analysis_parameters']['DataTransformation']
		IntrinsicModel     = self.choices['analysis_parameters']['IntrinsicModel']
		n_warmup,n_sampling,n_chains,n_thin,random_seed = self.choices['analysis_parameters']['n_warmup'],self.choices['analysis_parameters']['n_sampling'],self.choices['analysis_parameters']['n_chains'],self.choices['analysis_parameters']['n_thin'],self.choices['analysis_parameters']['random_seed']
		self.Rhat_check_params = get_parlabels(self.choices,return_parnames_only=True)

		#Initialisation of stan_data
		stan_data = {}
		stan_data['Nm'] = len(pblist)
		stan_data['Nc'] = stan_data['Nm']-1
		stan_data['zero_index']  = self.choices['analysis_parameters']['zero_index']
		stan_data['gamma_shape'] = 1 if self.choices['analysis_parameters']['include_LCshape']   else 0
		stan_data['gamma_res']   = 1 if self.choices['analysis_parameters']['include_residuals'] else 0
		analysis_parameters_list = ['muRVmin','muRVmax','RVsmin','RVsmax','disp_sigmaRV','a_sigma_mint','a_sigma_cint']
		for par in analysis_parameters_list:
			stan_data[par] = self.choices['analysis_parameters'][par]

		modelloader = HBM_preprocessor(self.choices, DF_M)

		#Add data for binned RV fit
		if 'BinnedRVFit' in modelloader.choices['analysis_parameters'] and modelloader.choices['analysis_parameters']['BinnedRVFit']:
			print ('Building Stan data from BinnedRVFit Parameters'+'\n'+'###'*10)
			DF_M = modelloader.get_BV_ordered_DF_M() #print (DF_M[0]['B']-DF_M[0]['V'])
			modelloader.get_RVbin_data()
			for key in ['N_RV_bins','N_GaussRV_dists','N_AVRVBeta_dists','N_AVRVSigmoid_dists','RV_bin_vec','RVstyle_per_bin','map_RVBin_to_RVHyp','fixed_RVs']:
				stan_data[key] = modelloader.__dict__[key]
				print (f"{key}:{stan_data[key]}")
			for key in ['flatRVsmin','flatRVsmax','a_sigma_beta0','a_sigma_beta1','a_sigma_muAVsigmoid','a_sigma_sigAVsigmoid']:
				stan_data[key]                           = modelloader.choices['analysis_parameters'][key]
				self.choices['analysis_parameters'][key] = modelloader.choices['analysis_parameters'][key]#If these values are only called upon as default, save in .choices for plotting
			print (f"B-V Boundaries are:{modelloader.choices['analysis_parameters']['BVbinboundaries']}\nNumber of SNe in each bin is: {modelloader.N_in_each_bin}")
			self.choices['analysis_parameters']['N_in_each_bin'] = modelloader.N_in_each_bin#Save this info for plotting

		#Incorporate Censored Data
		modelloader.get_censored_data()
		stan_data['S']  = modelloader.S
		stan_data['SC'] = modelloader.SC
		stan_data['BVerrs_Cens'] = modelloader.BVerrs_Cens
		stan_data['BVcutval']    = modelloader.BVcutval
		if self.choices['analysis_parameters']['CensoredData']:
			if 'BinnedRVFit' in modelloader.choices['analysis_parameters'] and modelloader.choices['analysis_parameters']['BinnedRVFit']:
				if modelloader.__dict__['N_RV_bins']>1:
					raise Exception('Havent yet full considered combination of RVbinning models and Censored Data. \nWhen applying censored data, best to use single set of population distributions')
			print (f"Incorporating censored data: Fitting {stan_data['S']-stan_data['SC']} SNe and {stan_data['SC']} Censored SNe")
			print (f"Censored SNe are: {modelloader.CensoredSNe}")
			print (f"These have {self.choices['additional_cut_parameters']['BVcutval']}<=B-V<{self.choices['analysis_parameters']['CensoredCut']}")
			print (f"Completely Excluded SNe are: {modelloader.ExcludedSNe}")

		#If the intrinsic model is fixed, load this up from another fit with HBM_savekey=fixint
		if self.choices['analysis_parameters']['fixint'] is not False:
			modelloader.get_intrinsic_model(self.FITSpath)
			stan_data['L_mint'] = modelloader.L_mint
			stan_data['FPC0']	= modelloader.FPC0
			stan_data['FPC1']	= modelloader.FPC1
			if self.choices['analysis_parameters']['AVprior']=='Gauss':
				stan_data['disp_sigmaAV'] = self.choices['analysis_parameters']['disp_sigmaAV']

		#Incorporate LC shape data measurements
		modelloader.get_dm15Bs()
		stan_data['dm15Bs']     = modelloader.dm15Bs
		stan_data['dm15B_errs'] = modelloader.dm15B_errs

		#Get dust law items
		modelloader.get_dustlaw()
		stan_data['xk'] = modelloader.xk
		if IntrinsicModel=='Deviations': stan_data['Mmatrix'] = modelloader.M_fitz_block
		else:							 stan_data['DelM']    = modelloader.dM_fitz_block

		#Get transformation matrix from mags to colours
		modelloader.get_CM_transformation_matrix()
		if IntrinsicModel=='Deviations' and DataTransformation!='mags':
			stan_data['CM'] = modelloader.CM

		#Get transformation matrix from 1 colours model to another set of colours
		if IntrinsicModel!='Deviations':
			modelloader.get_CC_transformation_matrices()
			stan_data['CC'] = modelloader.CC
			stan_data['CC_to_adj'] = modelloader.CC_to_adj

		#Get list of peak mags or colours, and measurement errors
		modelloader.get_transformed_data()
		if DataTransformation=='mags':
			stan_data['mags']      = modelloader.mags
			stan_data['mags_errs'] = modelloader.mags_errs
		else:
			stan_data['capps']      = modelloader.capps
			stan_data['capps_errs'] = modelloader.capps_errs

		#If skew_RV=True, update stan_data with skew_RV_disp_prior
		if 'skew_RV' in modelloader.choices['analysis_parameters'] and modelloader.choices['analysis_parameters']['skew_RV']:
				stan_data['skew_RV_disp_prior'] = modelloader.choices['analysis_parameters']['skew_RV_disp_prior']

		#If skew_int=True, update stan_data with skew_int_disp_prior
		if 'skew_int' in modelloader.choices['analysis_parameters'] and modelloader.choices['analysis_parameters']['skew_int']:
				stan_data['skew_int_disp_prior'] = modelloader.choices['analysis_parameters']['skew_int_disp_prior']

		#Make replica copies of data and fit multiplied sample
		stan_data = modelloader.multiply_dataset(stan_data)

		#Initialise Stan Model Samples with Sample Average Mags
		stan_init = modelloader.get_init(stan_data)

		#Assert size of data vectors matches asserted integer lengths, remove data not neeeded
		stan_data = modelloader.data_model_checks(stan_data)

		#Get Stan File
		modelloader.get_stan_file(stanpath=self.modelpath+'stan_files/')
		#Modify Stan File, e.g. changes to AVprior
		stan_file = modelloader.modify_stan_file()

		print (f'Beginning HBM_savekey fit: {savekey}')
		#Build Stan Model
		posterior = stan.build(stan_file, data=stan_data, random_seed=random_seed)
		for itrial in range(Ntrials):
			BREAK = True
			#Fit Model
			fit = posterior.sample(num_chains=n_chains, num_samples=n_sampling*(1+itrial)*(additive==False)+(n_sampling+adder)*(additive==True), num_warmup = n_warmup,init=stan_init)
			#Get samples as pandas df
			df  = fit.to_frame()

			if df.shape[0]>n_thin*n_chains:#Thin samples
				print (f'Thinning samples down to {n_thin} per chain to save on space complexity')
				Nthinsize = int(n_thin*n_chains)
				df = df.iloc[0:df.shape[0]:int(df.shape[0]/Nthinsize)]						#thin to e.g. 1000 samples per chain
				dfdict = df.to_dict(orient='list')											#change to dictionary so key,value is parameter name and samples
				fit = {key:np.array_split(value,n_chains) for key,value in dfdict.items()}	#change samples so split into chains

			#Fitsummary including Rhat valuess
			fitsummary = map_fit_summary(az.summary(fit))#feed dictionary into arviz to get summary stats of thinned samples
			#Check for high Rhat values
			for par in self.Rhat_check_params:
				if fitsummary.loc[par]['r_hat']>Rhat_threshold:
					BREAK = False
			if BREAK:	break
			else:		print (f'Completed {itrial+1}/{Ntrials}; Repeating fit with more samples because of Rhats:',{par:fitsummary.loc[par]['r_hat'] for par in self.Rhat_check_params})
		#Products of HMC fit
		FIT        = {'df':df,'fitsummary':fitsummary,'stan_data':stan_data,'choices':self.choices}
		#Save FIT
		if save:
			with open(f'{self.FITSpath}FIT{savekey}.pkl','wb') as f:
				pickle.dump(FIT,f)
		#Set attribute
		self.FIT   = FIT

		#Print Posterior summaries
		print ('~~~'*5+f'\n{savekey} FIT summary')
		for par in self.Rhat_check_params:
			with suppress(KeyError):
				print ('###'*5)
				print ('Par, median, std, 5%, 68%, 95% quantiles')
				print (par,df[par].median().round(2),df[par].std().round(2), df[par].quantile(0.05).round(2), df[par].quantile(0.68).round(2),df[par].quantile(0.95).round(2))
				print (f"Rhat = {fitsummary.loc[par]['r_hat']}")

	def plot_posterior_samples(self, returner=False):
		"""
		Plot Posterior Samples

		Method to take FIT.pkl object and plot posterior samples

		Parameters
		----------
		returner : bool (optional; default=False)
			if True, return Summary_Strs

		Returns
		----------
		Summary_Strs: dict
			key is parameter
			value is string of posterior summary statistic for use in LaTeX Tables
		"""
		savekey  = self.choices['analysis_parameters']['HBM_savekey']
		filename = f'{self.FITSpath}FIT{savekey}.pkl'
		with open(filename,'rb') as f:
			FIT = pickle.load(f)

		df         = FIT['df']
		fitsummary = FIT['fitsummary']
		NSNe       = FIT['stan_data']['S']
		NCens      = FIT['stan_data']['SC']

		pars,parnames,parlabels,bounds = get_parlabels(FIT['choices'])
		Rhats      = {par:fitsummary.loc[par]['r_hat'] for par in parnames}
		samples    = {par:np.asarray(df[par].values)   for par in parnames}
		print ('Rhats:',Rhats)

		#Corner Plot
		postplot = POSTERIOR_PLOTTER(samples, parnames, parlabels, bounds, Rhats, self.choices['plotting_parameters'])
		Summary_Strs = postplot.corner_plot()#Table Summary
		Summary_Strs['NSNe(NCens)'] = f"{NSNe-NCens}({NCens})"

		plotpath        =  self.plotpath+'Corner_Plot/'
		save,quick,show = [self.choices['plotting_parameters'][x] for x in ['save','quick','show']][:]
		finish_corner_plot(postplot.fig,postplot.ax,get_Lines(FIT['choices'],NSNe-NCens, NCens),save,quick,show,plotpath,savekey)

		if returner: return Summary_Strs

	def plot_mag_deviations(self):
		"""
		Plot Mag Deviations

		Plot to show the combination of dust and intrinsic colour in a population

		End Product(s)
		----------
		MagDeviations.pdf in self.plotpath
		"""
		#Initialisations
		DF_M = self.DF_M
		modelloader = HBM_preprocessor(self.choices, DF_M)
		pblist      = self.choices['preproc_parameters']['pblist']
		try:		l_eff_rest  = modelloader.get_leff_rest()
		except:		l_eff_rest = np.array([fset[pb].ave_wave for pb in pblist])
		errstr      = self.choices['preproc_parameters']['errstr']
		tref        = self.choices['preproc_parameters']['tilist'][self.choices['preproc_parameters']['tref_index']]
		logMcut     = self.choices['additional_cut_parameters']['logMcut']
		BVcutval    = self.choices['additional_cut_parameters']['BVcutval']
		FS          = self.choices['plotting_parameters']['FS']
		Nm          = len(pblist)
		linestyles  = ['-',':','-.','--']

		#Set up figure
		pl.figure(figsize=(9.6,7.2))
		pl.title("Magnitude Deviations from each SN's Mean Apparent Magnitude:\n"+r'$m^s_i-\langle m_i^s \rangle_i = \delta N^s_i - \langle \delta N_i^s \rangle_i + A^s_V\delta \xi^s_i$', fontsize=FS)
		mags = DF_M[tref][pblist].stack().transpose().values; magerrs = DF_M[tref][[pb+errstr for pb in pblist]].stack().transpose().values
		Nl = 100 ; RVs = np.array([1.5,2.5,3.5,4.5]) ; lams = np.linspace(l_eff_rest[0],l_eff_rest[-1],Nl) ; NCOL = 4
		#Plot deviations for each SN
		Bind = pblist.index('B') ; Vind = pblist.index('V')
		for s in range(DF_M[tref].shape[0]):
			sn   = DF_M[tref].index[s] ; mass = self.lcs[sn].meta['Mbest'] ; colour, highmass = get_mass_label(mass,self.choices)
			DL   = 100-200*(highmass is False)-100*(highmass is None)
			mean_mag = np.average(mags[s*Nm:(s+1)*Nm])
			vec      = (np.array([mags[s*Nm+_] for _ in range(Nm)])-mean_mag)
			mean_mag_err = (sum(magerrs[s*Nm:(s+1)*Nm]**2)**0.5)/Nm
			vec_err      = (np.array([magerrs[s*Nm+_]**2 for _ in range(Nm)])+mean_mag_err**2)**0.5
			for _ in range(Nm):
				pl.errorbar(l_eff_rest[_]+DL,vec[_],c=colour,linestyle='None',marker={True:'o',False:'x'}[mags[s*Nm+Bind]-mags[s*Nm+Vind]<BVcutval],alpha=0.45,capsize=2)
				if s==0:	pl.annotate(pblist[_],xy=(l_eff_rest[_]-200,0.79),weight='bold',fontsize=FS+4)
			pl.plot(l_eff_rest+DL,vec,c='black',linestyle='-',alpha=0.05)
		#Plot RV lines
		for iRV,RV in enumerate(RVs):
			xibar = np.average(extinction.fitzpatrick99(l_eff_rest,1,RV))
			xi    = extinction.fitzpatrick99(lams,1,RV)
			pl.plot(lams,(xi-xibar),label=f'$R_V$={RV}',linestyle=linestyles[iRV])
		#Plot masses legend
		pl.xlim(pl.gca().get_xlim())
		pl.scatter(1000,0,c=get_mass_label(logMcut*2,self.choices)[0],alpha=0.45,label='$\log M/M_{\odot}\geq 10$')
		pl.scatter(1000,0,c=get_mass_label(logMcut*0.5,self.choices)[0],alpha=0.45,label='$\log M/M_{\odot}<10$')
		try: pl.scatter(1000,0,c=get_mass_label(None,self.choices)[0],alpha=0.45,label='N/A')
		except: NCOL += -1
		pl.annotate(
			r"      $|B-V|<%s\,$mag        High Reddening  "%BVcutval,xy=(0.025, 0.05),xycoords="axes fraction",fontsize=FS-2,
			bbox=dict(boxstyle="round",alpha=0.25,facecolor = "white",edgecolor = "grey"))
		#Plot high-reddening legend
		pl.errorbar(4300,-0.7775,c='black',marker='o',alpha=0.45)
		pl.errorbar(8275,-0.7775,c='black',marker='x',alpha=0.45)
		pl.annotate('Fitzpatrick99 Dust Law',xy=(0.130,0.925),xycoords='axes fraction',fontsize=FS-1)
		pl.annotate('Host Galaxy Stellar Mass',xy=(0.575,0.925),xycoords='axes fraction',fontsize=FS-1)
		pl.legend(fontsize=FS-2,loc='upper center',ncol=NCOL,bbox_to_anchor=(0.54,0.925),columnspacing=0.985)
		pl.plot(list(pl.gca().get_xlim()),[0,0],c='black')
		#Finish Plot
		plotter = PLOTTER(self.choices['plotting_parameters'], self.plotpath)
		plotter.finish_plot(r'$\lambda (\AA)$',r'$\delta N_i - \langle \delta N_i \rangle_i + A_V \delta \xi_i$ (mag)',savename='MagDeviations.pdf')

	def plot_colour_corner(self,use_intrinsic_mean=False):
		"""
		Plot Colour Colour

		Function to plot colour colour corner plot

		Parameters
		----------
		use_intrinsic_mean : bool (optional; default=False)
			if False, centre dust projections and sample's 5% quantiles.
			Otherwise, specify the string name of the file to load intrinsic hyperparameters from

		Returns
		-----------
			Either saves plot, or shows plot, or both
		"""

		DF_M = self.DF_M
		modelloader = HBM_preprocessor(self.choices, DF_M)
		pblist      = self.choices['preproc_parameters']['pblist']
		try:		l_eff_rest  = modelloader.get_leff_rest()
		except:		l_eff_rest = np.array([fset[pb].ave_wave for pb in pblist])
		tilist      = self.choices['preproc_parameters']['tilist']
		errstr      = self.choices['preproc_parameters']['errstr']
		tref        = self.choices['preproc_parameters']['tilist'][self.choices['preproc_parameters']['tref_index']]
		Style       = self.choices['plotting_parameters']['ColourCornerStyle']
		Input       = self.choices['plotting_parameters']['ColourCornerInput']
		parstr      = 'parnames'

		def lam(f): return dict(zip(pblist,l_eff_rest))[f]
		def fitz(lam,RV): return extinction.fitzpatrick99(np.array([lam]),1,RV)

		grid_sets, grid_names = get_list_of_colours_for_corner(pblist,Style,Input=Input,parstr=parstr)
		#Plot Hyperparams
		colour    = 'C0'
		alph      = 0.2
		FS        = 15
		Nsim      = 100
		DC        = 0.1

		##########################
		#Get Intrinsic Mean Hyperparamter from Fiducial Analysis if use_intrinsic_mean!=False
		if use_intrinsic_mean is not False:
			file = copy.deepcopy(use_intrinsic_mean)
			use_intrinsic_mean = True
			assert(type(file)==str)
			with open(f"{self.FITSpath}FIT{file}.pkl",'rb') as f:
				FIT = pickle.load(f)
			df     = FIT['df']
			mu_int = np.array([df[f'FPC0.{i}'].median() for i in range(1,7)])
			CM = np.zeros((5,6))
			for _ in range(5):
				CM[_,_]   =  1
				CM[_,_+1] = -1
			c_int  = CM@mu_int
		##########################

		for ti in tilist:
			for gname in grid_names:
				#Set up grid
				col_list = grid_sets[gname]
				parnames = grid_sets[f"{gname}_{parstr}"]
				sfigx,sfigy = 2.5*len(col_list),2.5*len(col_list)
				fig,ax = pl.subplots(len(col_list),len(col_list),figsize=(sfigx,sfigy),sharex='col',sharey=False)
				#Initialisations
				Cmeans = {}
				dfBV        = (DF_M[ti]['B']-DF_M[ti]['V']).dropna().copy().values
				highBVbools = [not BV<self.choices['additional_cut_parameters']['BVcutval'] for BV in dfBV]
				for iax,col in enumerate(col_list):
					parname     = parnames[iax]
					df          = DF_M[ti].dropna().copy()
					c1s, c1errs = (df[col[0]]-df[col[1]]).values, ((df[col[0]+errstr]**2+df[col[1]+errstr]**2)**0.5).values
					chain       = np.sort(np.array([np.random.normal(c1,c1err) for c1,c1err in zip(c1s,c1errs) for _ in range(Nsim)]))
					data = PARAMETER(chain,col,r'$%s-%s$'%(col[0],col[1]),[min(chain)-2*np.std(chain),max(chain)+2*np.std(chain)],[None,None],1,iax,self.choices,3)
					data.get_xgrid_KDE()
					#Get Quantiles
					i50,x50,K50 = data.get_KDE_values(value=data.samp_median)
					i16         = data.get_KDE_values(value=chain[int((len(chain)-1)*0.16)],return_only_index=True)
					i84         = data.get_KDE_values(value=chain[int((len(chain)-1)*0.84)],return_only_index=True)+1
					i05,x05,K05 = data.get_KDE_values(value=chain[int((len(chain)-1)*0.05)])
					i95         = data.get_KDE_values(value=chain[int((len(chain)-1)*0.95)],return_only_index=True)+1
					#Plot KDEs
					ax[iax,iax].set_title(f"{str(parname).replace(' (mag)','')} = "+r"$%s^{+%s}_{-%s}$"%(round(data.samp_median,2),round(data.dfchain.par.quantile(0.84)-data.samp_median,2),round(data.samp_median-data.dfchain.par.quantile(0.16),2)),fontsize=FS)
					ax[iax,iax].plot(data.xgrid,data.KDE, color=colour)
					ax[iax,iax].plot(np.ones(2)*data.samp_median,[0,data.KDE[i50]],c=colour)
					ax[iax,iax].fill_between(data.xgrid[i16:i84],np.zeros(i84-i16),data.KDE[i16:i84],color=colour,alpha=alph)
					ax[iax,iax].fill_between(data.xgrid[i05:i95],np.zeros(i95-i05),data.KDE[i05:i95],color=colour,alpha=alph)
					ax[iax,iax].set_ylim([0,None])
					ax[iax,iax].set_xlim([min(c1s)-0.1, max(c1s)+0.1])
					ax[iax,iax].set_yticks([])
					if use_intrinsic_mean:	Cmeans[col]  = [c_int[iax], min(c1s)-DC, max(c1s)+DC ]
					else:					Cmeans[col]  = [x05,        min(c1s)-DC, max(c1s)+DC ]
					ax[len(col_list)-1,iax].set_xlabel(parname,fontsize=FS)
					if iax>0:
						ax[iax,0].set_ylabel(parname,fontsize=FS)
					for iaxbelow in range(iax):
						ax[iaxbelow,iax].axis('off')
					if iax==len(col_list)-1:
						ax[iax,iax].tick_params(labelsize=FS)
				#Plot 2D colour scatter panels
				lobool=True ; hibool=True ; highmassbool = True; lowmassbool = True
				for ix,c1 in enumerate(col_list[:-1]):
					for ic2,c2 in enumerate(col_list[ix+1:]):
						if ix==0 and ic2==0 and use_intrinsic_mean:
							ax[0,0].scatter(Cmeans[c1][0]-10, Cmeans[c2][0],marker='*', color='black',s=100,zorder=100,label='Intrinsic \nPop. Mean')
							ax[0,0].legend(fontsize=FS,bbox_to_anchor=(5.225, -2.17),loc='upper right',handletextpad=0.1)
						iy  = ic2+ix+1
						df  = DF_M[ti].dropna().copy()
						sns = list(df.index)
						c1s, c1errs = (df[c1[0]]-df[c1[1]]).values, ((df[c1[0]+errstr]**2+df[c1[1]+errstr]**2)**0.5).values
						c2s, c2errs = (df[c2[0]]-df[c2[1]]).values, ((df[c2[0]+errstr]**2+df[c2[1]+errstr]**2)**0.5).values
						if ix==0:
							ax[iy,ix].tick_params(labelsize=FS)
						else:
							ax[iy,ix].set_yticks(ax[iy,0].get_yticks())
							ax[iy,ix].set_yticklabels([])
							ax[iy,ix].set_ylim(ax[iy,0].get_ylim())
						if iy==len(col_list)-1:
							ax[iy,ix].tick_params(labelsize=FS)
						for cc1,cc2,c1e,c2e,highBVbool,sn in zip(c1s,c2s,c1errs,c2errs,highBVbools,sns):
							ccc = 'C0'
							ccc,highmass = get_mass_label(self.lcs[sn].meta['Mbest'],self.choices)
							if lobool or hibool:
								if highBVbool and hibool and iy==1 and ix==0:
									ax[iy,ix].errorbar(cc1,cc2,xerr=c1e,yerr=c2e,color=ccc,marker={False:'o',True:'s'}[highBVbool],capsize=0,alpha={False:0.09,True:0.45}[highBVbool],markersize=6,linestyle='none',label='High Reddening')
									hibool = False
								if not highBVbool and lobool and iy==1 and ix==0 and not hibool:
									ax[iy,ix].errorbar(cc1,cc2,xerr=c1e,yerr=c2e,color=ccc,marker={False:'o',True:'s'}[highBVbool],capsize=0,alpha={False:0.09,True:0.45}[highBVbool],markersize=6,linestyle='none',label=r'$|B-V|<%s\,$mag'%self.choices['additional_cut_parameters']['BVcutval'])
									lobool = False
							else:
								if highmass and highmassbool and iy==2 and ix==0 and highBVbool:
									print (iy,ix,'high')
									ax[iy,ix].errorbar(cc1,cc2,xerr=c1e,yerr=c2e,color=ccc,marker={False:'o',True:'s'}[highBVbool],capsize=0,alpha={False:0.09,True:0.45}[highBVbool],markersize=6,linestyle='none',label='High Mass')
									highmassbool = False
								elif not highmass and lowmassbool and iy==2 and ix==0 and highBVbool:
									print (iy,ix,'low')
									ax[iy,ix].errorbar(cc1,cc2,xerr=c1e,yerr=c2e,color=ccc,marker={False:'o',True:'s'}[highBVbool],capsize=0,alpha={False:0.09,True:0.45}[highBVbool],markersize=6,linestyle='none',label='Low Mass')
									lowmassbool  = False
								else:
									ax[iy,ix].errorbar(cc1,cc2,xerr=c1e,yerr=c2e,color=ccc,marker={False:'o',True:'s'}[highBVbool],capsize=0,alpha={False:0.09,True:0.45}[highBVbool],markersize=6,linestyle='none')
						ax[1,0].legend(fontsize=FS,bbox_to_anchor=(5.225,-0.01),loc='upper right')
						ax[2,0].legend(fontsize=FS,bbox_to_anchor=(5.225, 1.512),loc='upper right')
						X = np.linspace(Cmeans[c1][1],Cmeans[c1][2],10)
						RVS = [1.5,2.5,3.5][::-1]
						for iRV,RV in enumerate(RVS):
							gradient = (fitz(lam(c2[0]),RV)-fitz(lam(c2[1]),RV))/(fitz(lam(c1[0]),RV)-fitz(lam(c1[1]),RV))
							Y = Cmeans[c2][0] + (X-Cmeans[c1][0])*gradient
							if use_intrinsic_mean: ax[iy,ix].scatter(Cmeans[c1][0], Cmeans[c2][0], marker='*', color='black',s=100,zorder=100)
							ax[iy,ix].plot(X,Y,color=f"C{3-iRV}",linestyle=[':','-','-.'][::-1][iRV],label=f'$R_V={RV}$'*(iy==len(col_list)-1 and ix==len(col_list)-2))#*(ix==0 and iy==1))
							ax[iy,ix].plot([Cmeans[c1][0],Cmeans[c1][0]+fitz(lam(c1[0]),RV)-fitz(lam(c1[1]),RV)],[Cmeans[c2][0],Cmeans[c2][0]+fitz(lam(c2[0]),RV)-fitz(lam(c2[1]),RV)],color=f"C{3-iRV}",linestyle=[':','-','-.'][iRV-iRV+1],linewidth=5,alpha=0.4)
							if iy==len(col_list)-1 and ix==len(col_list)-2:
								ax[iy,ix].legend(fontsize=FS,bbox_to_anchor=(2,2.75),loc='upper right')
				#Fix up axes
				for col in range(len(col_list)):
					increments = ax[len(col_list)-1,col].get_xticks()
					step_order = -int(np.floor(np.log10(np.average(increments[1:]-increments[:-1]))))
					ax[len(col_list)-1,col].set_xticklabels(np.round(ax[len(col_list)-1,col].get_xticks(),step_order))
				fig.subplots_adjust(top=0.9)
				fig.subplots_adjust(wspace=0.075, hspace=0.075)
				#Lines annotations
				delta_y = 0.15
				Lines = get_Lines(self.choices, DF_M[tref].shape[0], 0, posterior=False)
				pl.annotate(r'Peak Colour Measurements', xy=(1,len(ax)+delta_y),xycoords='axes fraction',fontsize=FS+2,color='black',weight='bold',ha='right')
				pl.annotate(Lines[0], xy=(1,len(ax)),xycoords='axes fraction',fontsize=FS+2,color='black',weight='bold',ha='right')
				for counter,line in enumerate(Lines[1:]):
					pl.annotate(line, xy=(1,len(ax)-delta_y*(counter+1)),xycoords='axes fraction',fontsize=FS+1,color=colour,weight='bold',ha='right')
				#Save and/or show plot
				if self.choices['plotting_parameters']['save']:
					savefile = f"{self.plotpath}ColourCorner/{gname}.pdf"
					ensure_folders_to_file_exist(savefile)
					pl.savefig(savefile,bbox_inches='tight')
				if self.choices['plotting_parameters']['show']:
					pl.show()

	def plot_int_model_posteriors(self,PAR='muRV',FITS=None):
		"""
		Plot Intrinsic Model Posteriors

		Plotting function to show different Hyperparameter Posteriors under different Intrinsic Models

		Parameters
		----------
		par : str
			Name of parameter to plot
		FITS : dict (optional; default is None)
			posteriors to loop over
			{key:value} pairs are 'filename':{'label':label}
			if None, defaults to x4 models Deviations,Adj,BX,XH

		End Product(s)
		----------
		Plot of muRV posteriors
		"""
		if FITS is None:
			#FITS = {
			#'AVexp_Cens1.0':{'label':'Deviations'},
			#'AdjCols_Cens1.0':{'label':'Adjacent Colours'},
			#'BXCols_Cens1.0':{'label':'$B-X$ Colours'},
			#'XHCols_Cens1.0':{'label':'$X-H$ Colours'},
			#}
			FITS = {
			'AVExp_BVcut1.0':{'label':'Deviations'},
			'AdjCols_BVcut1.0':{'label':'Adjacent Colours'},
			'BXCols_BVcut1.0':{'label':'$B-X$ Colours'},
			'XHCols_BVcut1.0':{'label':'$X-H$ Colours'},
			}
		FS = self.choices['plotting_parameters']['FS']
		pl.figure(figsize=(8,6))
		pl.title('Choice of Intrinsic SN Model',fontsize=FS)
		pl.ylabel('Posterior Densities',fontsize=FS)
		for counter,file in enumerate(FITS):
			colour  = ['C0','C1','C2','C3'][counter]
			ls      = ['-','--',':','-.'][counter]
			with open(f"{self.FITSpath}FIT{file}.pkl",'rb') as f:
				FIT = pickle.load(f)
			df  = FIT['df']
			pars,parnames,parlabels,bounds = get_parlabels(FIT['choices'])
			ipar = pars.index(PAR)
			par  = parnames[ipar]
			lims = [df[par].min(),df[par].max()]
			if par=='mu_RV': lims = [1.5,5]
			samps = PARAMETER(df[par],par,parlabels[ipar],lims,bounds[ipar],0,0,{},2)
			samps.get_xgrid_KDE()
			pl.plot(samps.dfchain.quantile(0.5).values[0]*np.ones(2),[0,np.amax(samps.KDE)],c=colour	,linewidth=2)
			pl.plot(samps.xgrid,samps.KDE,alpha=1,color=colour,linestyle=ls,linewidth=3,
					#label=FITS[file]['label']+'\n'+r'$%s\pm%s$'%(round(samps.dfchain.quantile(0.5).values[0],2),round(samps.dfchain.std().values[0],2))
					label=FITS[file]['label']+'\n'+r'$%s^{+%s}_{-%s}$'%(round(samps.dfchain.quantile(0.5).values[0],2),round(samps.dfchain.quantile(0.84).values[0]-samps.dfchain.quantile(0.5).values[0],2),round(samps.dfchain.quantile(0.50).values[0]-samps.dfchain.quantile(0.16).values[0],2)  )
					)
			pl.xlim(lims)
		pl.tick_params(labelsize=FS)
		pl.yticks([])
		pl.legend(fontsize=FS-2,framealpha=0.8,loc='upper right',title = r'%s Inference'%(parlabels[ipar])+False*'Reference frame for modelling intrinsic chromatic variations',title_fontsize=FS-2,ncol=1,columnspacing=0.7)
		pl.gca().set_ylim([0,1*pl.gca().get_ylim()[1]])
		pl.xlabel(parlabels[ipar],fontsize=FS)
		pl.tight_layout()
		pl.savefig(f"{self.plotpath}SensitivityToIntrinsicModel.pdf",bbox_inches='tight')
		pl.show()

def get_edit_dict(choices,CYCLE_DICT,HBM_savekey):
	"""
	Get Edit Dictionary

	Function for running multiple analysis choices
	Edits the Bird-Snack choices
	Takes in changes common to all CYCLE_DICT['RUNS']
	and Individual changes specific to HBM_savekey : {'newdict':{changes}}

	Parameters
	----------
	choices : dict
		default dictionary of choices from .yaml
	CYCLE_DICT :dict
		{'COMMON_CHANGES':{'newdict':{key,value pairs}, 'HBMappender':str },
		 'RUNS' : {
			'HBM_savekey' : {'newdict' : {key,value pairs},
									'label': str}
		 }}
	HBM_savekey : str
		Name of an analysis variant in CYCLE_DICT['RUNS']

	Returns
	----------
	edit_dict : dict
		contains updated analysis choices and HBM_savekey
	"""
	#Initialise edit_dict
	edit_dict = {glob_key:{} for glob_key in choices if glob_key!='rootpath'}
	#Get changes
	COMMON_CHANGES       = CYCLE_DICT['COMMON_CHANGES']['newdict']
	HBM_savekey_appender = CYCLE_DICT['COMMON_CHANGES']['HBMappender']
	RUNS                 = CYCLE_DICT['RUNS']
	newdict              = RUNS[HBM_savekey]['newdict']
	#Firstly Implement Common Changes
	for key in COMMON_CHANGES:
		for glob_key in edit_dict:
			for kkey,vvalue in choices[glob_key].items():
				if kkey==key:
					edit_dict[glob_key][key] = COMMON_CHANGES[key]
	#Then Implement Specific Changes; where adding to dictionary make sure both common and specific are retained
	for key in newdict:
		for glob_key in edit_dict:
			for kkey,vvalue in choices[glob_key].items():
				if kkey==key:
					if type(newdict[key])==dict and key in edit_dict[glob_key]:
						edit_dict[glob_key][key] = {**edit_dict[glob_key][key],**newdict[key]}
					else:
						edit_dict[glob_key][key] = newdict[key]
	if HBM_savekey_appender!='':	edit_dict['analysis_parameters']['HBM_savekey'] = f"{HBM_savekey}_{HBM_savekey_appender}"
	else:							edit_dict['analysis_parameters']['HBM_savekey'] =    HBM_savekey
	#Return changes
	return edit_dict

def map_fit_summary(fitsummary):
	"""
	Map Fit Summary

	Stan chains save e.g. mu_RV.1, mu_RV.2,
	Whereas fitsummary save as e.g. mu_RV[0], mu_RV[1]
	This function maps fitsummary to be similar to Stan chains

	Parameters
	----------
	fitsummary : pandas df
		contains e.g. Rhat values of MCMC chain

	End Product(s)
	----------
	fitsummary with renamed index
	"""
	def map_str(col):
		pattern = r'(.*)\[([0-9]+)\]'
		try:
			x = list(re.search(pattern,col).groups())
			return f"{x[0]}.{int(x[1])+1}"
		except:
			return col

	fitsummary = fitsummary.rename(index={column:map_str(column) for column in list(fitsummary.index)})
	return fitsummary
