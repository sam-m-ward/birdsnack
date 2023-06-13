"""
SNANA Functions Module

Module containing functions to preprocess snana lcs


Contains:
--------------------
get_snana_foldername(snpy_params):
	get foldername for specific analysis variant where snana lc is saved

correct_lc_data_SNR_boostdict(lc,snrcut,error_boost_dict):
	correct snana lc by 1) applying SNR cut 2) Boosting errors where defined

get_mag(lc):
	create magnitudes column from flux column

trim_filters(lc, interpflts):
	trim filters in lc to interpolation filters

set_lcmeta_ordered_flts(lc):
	set lc.meta['flts'] and lc.meta['lams'] in order of wavelengths

update_lcmetadata(lc,dfsn,snpy_product)
	update lc.meta with mass, spectroscopic subtype, Tmax measurements etc.


--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""
import numpy as np

def get_snana_foldername(snpy_params):
	"""
	Get snana LC Folder Name

	Function to get folder name for snana LC

	Parameters
	----------
	snpy_params: dict
	Choices for SNooPy analysis
	  keys are [snrcut,snpy_model,snpy_shape,apply_EBVMW_corrections,RVMW,dustlaw,insert_pre_defined_Tmax,apply_K_corrections,mangle,interpflts]

	 Returns
	 ----------
	 str:
		foldername, used for specific analysis variant
	"""
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

def correct_lc_data_SNR_boostdict(lc,snrcut,error_boost_dict):
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

def get_mag(lc):
	"""
	Get Mag Columns

	Returns light curve object including new mag and magerr columns

	Parameters
	----------
	lc: :py:class:`astropy.table.Table`
		light curve object

	Returns
	----------
	lc with mag and magerr columns
	"""
	lc['mag']    = -2.5*np.log10(lc['fluxcal']) + lc['zptmag']
	lc['magerr'] = np.fabs((2.5/np.log(10))*lc["fluxcalerr"]/lc["fluxcal"])
	return lc

def trim_filters(lc, interpflts):
	"""
	Trim Filters

	Retain only filters in interpolation filters

	Parameters
	----------
	lc: :py:class:`astropy.table.Table`
		light curve object

	interpflts : str
		string of interpolation filters, either 'all' or e.g. 'uBgVriYJH'

	Returns
	----------
	lc with only interpflts subset
	"""
	if interpflts!='all':
		for flt in set(lc["flt"]):
			if flt not in interpflts:
				lc = lc[lc['flt']!=flt]
	return lc

def set_lcmeta_ordered_flts(lc):
	"""
	Set LCmeta Ordered Filters

	Sets lc.meta['flts'] and lc.meta['lams'] so that they are ordered by wavelength from blue to red

	Parameters
	----------
	lc: :py:class:`astropy.table.Table`
		light curve object

	Returns
	----------
	lc with updated .meta
	"""
	from snpy import fset
	lamC = {fset[flt].ave_wave:flt for flt in list(set(lc['flt']))}
	lams = list(lamC.keys())
	lams.sort()
	ordered_lamC = {lamC[key]:key for key in lams}
	lc.meta['flts'] = list(ordered_lamC.values())
	lc.meta['lams'] = list(ordered_lamC.keys())
	return lc

def update_lcmetadata(lc,dfsn,snpy_product):
	"""
	Update LC metadata

	Updates lc.meta to include e.g. mass, spectroscopic sub-classification, Tmax estimates etc.

	Parameters
	----------
	lc: :py:class:`astropy.table.Table`
		light curve object

	dfsn: pandas df
		row(s) where dfmeta['SN']==sn

	snpy_product: dict
		dictionary with various SNooPy correction items, including Tmax estimates

	Returns
	----------
	lc with updated .meta
	"""
	try:
		lc.meta['Mbest'] = dfsn['Mbest'].mean()
	except:
		Ms = [float(m) for m in dfsn['Mbest'] if m!='--']
		if Ms==[]: 	lc.meta['Mbest'] = None
		else:		lc.meta['Mbest'] = np.average(Ms)


	spectypes = list(dfsn['SubtypeIa'].values)
	if spectypes==[]:
		lc.meta['Spectral_Type']=None
	elif spectypes.count(spectypes[0])==len(spectypes):
		lc.meta['Spectral_Type']=spectypes[0]
	else:
		lc.meta['Spectral_Type']='_'.join(spectypes)

	#Set GP restframe Tmax estimate to None, not yet estimated
	lc.meta['Tmax_GP_restframe']     = None
	lc.meta['Tmax_GP_restframe_std'] = None
	#Load up the SNooPy fitted Tmax estimate
	try:
		lc.meta['Tmax_snpy_fitted']  = snpy_product['Tmax_snpy_fitted']
	except:
		lc.meta['Tmax_snpy_fitted']  = None

	return lc
