"""
SNooPy Corrections Module

Module containing SNOOPY_CORRECTIONS class and additional functions
Used to transform snpy object with corrections, saves Tmax_Kcorrdict, snpy_product, and snana lc files


Contains:
--------------------
get_Tmax_Kcorrdict_foldername:
	get folder name for Tmax and Kcorrections dictionary, pre-computed for all analysis variants in one go to save on time complexity

get_snana_foldername:
	get foldername for specific analysis variant where snana lc is saved

SNOOPY_CORRECTIONS class:
	inputs: choices, outputdir

	Methods are:
		get_fitted_snpy_Tmax_Kcorrdict(fitbands,obs_to_rest)
		correct_data_with_snpy_MWayExtinction_and_K_corrections()
		correct_sn(sn, SNsnpy)
		convert_snsnpy_to_snana()
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
	if interpflts!='all':
		for flt in set(lc["flt"]):
			if flt not in interpflts:
				lc = lc[lc['flt']!=flt]
	return lc

def set_lcmeta_ordered_flts(lc):
	from snpy import fset
	lamC = {fset[flt].ave_wave:flt for flt in list(set(lc['flt']))}
	lams = list(lamC.keys())
	lams.sort()
	ordered_lamC = {lamC[key]:key for key in lams}
	lc.meta['flts'] = list(ordered_lamC.values())
	lc.meta['lams'] = list(ordered_lamC.keys())
	return lc

def update_lcmetadata(lc,dfsn):
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
	'''
	lc.meta['Tmax_GP_restframe']     = None
			lc.meta['Tmax_GP_restframe_std'] = None
			try:
				lc.meta['Tmax_snpy_obsframe']   = snpy_products[sn]['Tmax_snpy_obsframe']
			except:
				lc.meta['Tmax_snpy_obsframe']   = None
			try:
				lc.meta['Tmax_GP_obsframe']   = snpy_products[sn]['Tmax_GP_obsframe']
			except:
				lc.meta['Tmax_GP_obsframe']   = None
			try:
				lc.meta['Tmax_snpy_restframe']  = snpy_products[sn]['Tmax_snpy_restframe']
			except:
				lc.meta['Tmax_snpy_restframe']  = None
			try:
				lc.meta['Tmax_snpy_fitted']  = snpy_products[sn]['Tmax_snpy_fitted']
			except:
				lc.meta['Tmax_snpy_fitted']  = None
	'''
