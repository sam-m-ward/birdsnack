"""
SNANA Functions Module

Module containing functions to preprocess/manipulate snana lcs

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

update_lcmetadata(lc,dfsn,snpy_product):
	update lc.meta with mass, spectroscopic subtype, Tmax measurements etc.

get_lcf(lc,flt):
	trims lc to a single filter, flt

create_phase_column(lc,Tmax,phasemin=None,phasemax=None):
	creates phase column using Tmax and phase boundaries

get_data_availability_bool(lcf,Nbeforemax,Nbeforegap,Naftermax,Naftergap,tref=0,tcol='phase',local=False):
	check data availability around time of reference

get_time_lc_arrays(lcf,mjd_or_phase='mjd',flux_or_mag='flux'):
	get time, bright, brighterr arrays
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
	#Insert Mass metadata
	try:
		lc.meta['Mbest'] = dfsn['Mbest'].mean()
	except:
		Ms = [float(m) for m in dfsn['Mbest'] if m!='--']
		if Ms==[]: 	lc.meta['Mbest'] = None
		else:		lc.meta['Mbest'] = np.average(Ms)

	#Insert Spectype metadata
	spectypes = list(dfsn['SubtypeIa'].values)
	if spectypes==[]:
		lc.meta['Spectral_Type']=None
	elif spectypes.count(spectypes[0])==len(spectypes):
		lc.meta['Spectral_Type']=spectypes[0]
	else:
		lc.meta['Spectral_Type']='_'.join(spectypes)

	#Initialise GP restframe Tmax estimate to None, not yet estimated
	lc.meta['Tmax_GP_restframe']     = None
	lc.meta['Tmax_GP_restframe_std'] = None

	#Insert SNooPy fitted Tmax estimate
	try:
		lc.meta['Tmax_snpy_fitted']  = snpy_product['Tmax_snpy_fitted']
	except:
		lc.meta['Tmax_snpy_fitted']  = None

	return lc


def get_lcf(lc,flt):
	"""
	Get Light Curve in filter f

	Returns the light curve by applying a filter mask

	Parameters
	----------
	lc: :py:class:`astropy.table.Table`
		light curve object

	flt: str
		passband to keep

	Returns
	----------
	lcf: :py:class:`astropy.table.Table`
		light curve object in passband
	"""
	return lc[lc["flt"]==flt]

def create_phase_column(lc,Tmax,phasemin=None,phasemax=None):
	"""
	Create Phase Column

	Returns light curve object including new phase column, cut on phasemin<=phase<=phasemax

	Parameters
	----------
	lc: :py:class:`astropy.table.Table`
		light curve object

	Tmax: float
		Time of maximum brightness in measured in some reference frame band

	phasemin: float or int (optional; default=None)
		Apply cut on minimum phase

	phasemax: float or int (optional; default=None)
		Apply cut on maximum phase

	Returns
	----------
	lc with phase column, cut so that phasemin<=phase<=phasemax
	"""
	if Tmax is not None:
		lc['phase'] = (lc['mjd']-Tmax)/(1+lc.meta['REDSHIFT_HELIO'])
	if phasemin is not None:
		lc = lc[lc['phase']>=phasemin]
	if phasemax is not None:
		lc = lc[lc['phase']<=phasemax]
	return lc


def get_data_availability_bool(lcf,Nbeforemax,Nbeforegap,Naftermax,Naftergap,tref=0,tcol='phase',local=False):
    """
    Trim Light Curve in Single Filter f

    Returns True to trim or False to Retain

    Parameters
    ----------
    lcf: :py:class:`astropy.table.Table`
        light curve objects

    Nbeforemax,Naftermax: floats or ints or lists
        if float or int, the Number of data points before or after tref required to retain Ndata<Nmax is trim==True
        if list, iterate over according to Ngaps (below)

    Nbeforegap,Naftergap: floats or ints or lists
        if float or int, is the time window (e.g. 5days) where data must reside in (before or after tref)
        if list, iterate over (e.g. 1 point within 5 days, 3 points within 10 days)

    tref: float or int (optional, default=0)
        time or phase e.g. t=0 is Tmax

    tcol: str (optional; default='phase')
        either phase or mjd, choose column for delta t

    local: bool (optional; default=False)
        if False, then e.g. Nbeforegap=-10 is -10 days in absolute terms
        if True, then e.g. Nbeforegap=-10 is 10 days before the reference phase, e.g. tref=15, so require data point in window 5-15 days

    """
    if local:
        lcf[tcol] = copy.deepcopy(lcf[tcol]-tref)
        tref      = 0
    def Nwindow_before(lcf,L,U=tref,colname=tcol): return len(lcf[ (-L<lcf[colname]) & (lcf[colname]<U) ])
    def Nwindow_after( lcf,U,L=tref,colname=tcol): return len(lcf[ ( L<lcf[colname]) & (lcf[colname]<U) ])

    if type(Nbeforegap)==float or type(Nbeforegap)==int:
        Ndata_beforemax = Nwindow_before(lcf,Nbeforegap)
        if Ndata_beforemax<Nbeforemax:
            #print ('Hey1')
            return True
    elif type(Nbeforegap)==list:
        Ndata_beforemax_lens = [Nwindow_before(lcf,Nbg) for Nbg in Nbeforegap]
        for i,Ndb in enumerate(Ndata_beforemax_lens):
            if Ndb<Nbeforemax[i]:
                #print ('Hey2', Ndb,Nbeforemax[i])
                return True

    if type(Naftergap)==float or type(Naftergap)==int:
        Ndata_aftermax  = Nwindow_after(lcf,Naftergap)
        if Ndata_aftermax<Naftermax: #e.g. if less than 1 data point before max within 5 days
            #print ('Hey3')
            return True

    elif type(Nbeforegap)==list:
        Ndata_aftermax_lens = [Nwindow_after(lcf,Nag) for Nag in Naftergap]
        for i,Ndb in enumerate(Ndata_aftermax_lens):
            if Ndb<Naftermax[i]:
                #print ('Hey4',Ndb,Naftermax[i], Ndata_aftermax_lens)
                return True
    return False

def get_time_lc_arrays(lcf,mjd_or_phase='mjd',flux_or_mag='flux'):
	"""
	Get time lc arrays

	Simple function to get light curve object columns

	Parameters
	----------
	lcf: :py:class:`astropy.table.Table`
		light curve object

	mjd_or_phase: 'mjd' or 'phase' or 'both' or 'neither' (optional; default='mjd')
		string to indicate what time arrays to output

	flux_or_mag: 'flux' or 'mag' or 'both' or 'neither' (optional; default='flux')
		string to indicate what brightness (and brightness error) arrays to output

	Returns
	----------
		some combination of mjd,phase,flux,fluxerr,mag,magerr (in that order)
	"""
	if mjd_or_phase=='mjd':
		yield lcf["mjd"]
	elif mjd_or_phase=='phase':
		yield lcf["phase"]
	elif mjd_or_phase=='both':
		yield lcf["mjd"],lcf["phase"]
	elif mjd_or_phase=='neither':
		yield
	else:
		raise Exception ("Enter either 'mjd','phase','both','neither' for keyword mjd_or_phase")

	if flux_or_mag=='flux':
		yield from (lcf["fluxcal"],lcf["fluxcalerr"])
	elif flux_or_mag=='mag':
		yield from (lcf["mag"],lcf["magerr"])
	elif flux_or_mag=='both':
		yield from (lcf["fluxcal"],lcf["fluxcalerr"],lcf["mag"],lcf["magerr"])
	elif flux_or_mag=='neither':
		yield
	else:
		raise Exception ("Enter either 'flux','mag','both','neither' for keyword flux_or_mag")
