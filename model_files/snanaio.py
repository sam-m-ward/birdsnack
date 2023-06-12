# -*- coding: UTF-8 -*-
"""
This module has been copied directly from the BayeSN public release github on 25/08/22
see https://github.com/bayesn/bayesn-public; All credits to BayeSN team

I/O methods. All the submodules of the BayeSNmodel package use this module for
almost all I/O operations.
"""

from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
import warnings
from collections import OrderedDict
import os
import json
import numpy as np
import astropy.table as at
import pkg_resources
import h5py
from six.moves import range
import sncosmo
import glob
import time

def _read_ascii(filename, format='ascii', **kwargs):
    """
    Read ASCII files

    Read space separated ASCII file, with column names provided on first line
    (leading ``#`` optional). ``kwargs`` are passed along to
    :py:func:`astropy.table.Table.read`.

    Parameters
    ----------
    filename : str
        Filename of the ASCII file. Column names must be provided on the first
        line.
    format : str, optional
        Format string for the file - ascii normally sufficies
    kwargs : dict
        Extra options, passed directly to :py:func:`astropy.table.Table.read`

    Returns
    -------
    out : :py:class:`numpy.recarray`
        Record array with the data. Field names correspond to column names in
        the file.

    See Also
    --------
    :py:func:`astropy.table.Table.read`
    """

    indata = at.Table.read(filename, format=format, **kwargs)
    return indata


# create some aliases
# these exist so we can flesh out full functions later
# with different formats if necessary for different sorts of data
read_pbmap     = _read_ascii
"""Read passband obsmode mapping table - wraps :py:func:`_read_ascii`"""


def get_pkgfile(infile):
    """
    Returns the full path to a file inside the :py:mod:`BayeSNmodel` package

    Parameters
    ----------
    infile : str
        The name of the file to set the full package filename for

    Returns
    -------
    pkgfile : str
        The path to the file within the package.

    Raises
    ------
    IOError
        If the ``pkgfile`` could not be found inside the :py:mod:`BayeSNmodel` package.

    Notes
    -----
        This allows the package to be installed anywhere, and the code to still
        determine the location to a file included with the package, such as the
        model grid file.
    """

    pkgfile = pkg_resources.resource_filename('BayeSNmodel',infile)

    if not os.path.exists(pkgfile):
        message = 'Could not find package file {}'.format(pkgfile)
        raise IOError(message)
    return pkgfile


def read_model_grid(grid_file=None, grid_name=None):
    """
    Read the Hsiao grid file


    Parameters
    ----------
    grid_file : None or str
        Filename of the Hsiao model grid HDF5 file. If ``None`` reads the
        ``hsiao.hdf5`` file included with the :py:mod:`BayeSNmodel`
        package.
    grid_name : None or str
        Name of the group name in the HDF5 file to read the grid from. If
        ``None`` uses ``default``

    Returns
    -------
    phase : array-like
        The phase array of the grid with shape ``(nphase,)``
    wave : array-like
        The wavelength array of the grid with shape ``(nwave,)``
    flux : array-like
        The DA white dwarf model atmosphere flux array of the grid.
        Has shape ``(nphase, nwave)``

    Notes
    -----
        There are no easy command line options to change this deliberately
        because changing the grid file essentially changes the entire model,
        and should not be done lightly, without careful comparison of the grids
        to quantify differences.
    """

    if grid_file is None:
        grid_file = os.path.join('templates','hsiao.h5')

    # if the user specfies a file, check that it exists, and if not look inside the package directory
    if not os.path.exists(grid_file):
        grid_file = get_pkgfile(grid_file)

    if grid_name is None:
        grid_name = "default"

    with h5py.File(grid_file, 'r') as grids:
        try:
            grid = grids[grid_name]
        except KeyError as e:
            message = '{}\nGrid {} not found in grid_file {}. Accepted values are ({})'.format(e, grid_name,\
                    grid_file, ','.join(list(grids.keys())))
            raise ValueError(message)

        phase = grid['phase'][()].astype('float64')
        wave  = grid['wave'][()].astype('float64')
        flux  = grid['flux'][()].astype('float64')

    return phase, wave, flux


def read_snana_lcfile(lcfile, sampname=None, metafile=None):
    '''
    Read light curves for a specific SN sample

    Parameters
    ----------
    lcfile : str
        Full path to the SNANA-format data file with the light curve
    sampname : str, optional
        Specify sample name of SNe to read (CfA/CSP/LOSS/LCO/misc)
        This allows specific pre-processing based on sample name
    metafile : str, optional
        Full path to a file which may contain additional metadata
        about this supernova. This file must have an `SNID` column.
        If our supernova's `SNID` matches a row in `metafile`, all
        other columns in that row of `metafile` will be added to
        the metadata of the returned astropy Table. Any info found
        in `metafile` will take precedence over metadata in the
        SNANA file. This mechanism can be used to provide
        custom/improved metadata.

    Returns
    -------
    snname  : a string name for the SN, parsed from the file, or filename
    outdata : astropy Table with the light curve data, corresponding to ``snname``

    '''
    try:
        meta, lcdata = sncosmo.read_snana_ascii(lcfile, default_tablename='OBS')
    except Exception as e:
        # if we cannot load something, give up and move on
        message = 'Could not load lcfile {}'.format(lcfile)
        warnings.warn(message)
        return None, None

    # set a default name for the object in case SNID isn't present in the SNANA file
    default_snname = os.path.basename(lcfile).split('_')[0]
    if not default_snname.startswith('snf'):
        default_snname = default_snname.replace('SN','').replace('sn','')
    snname = meta.get('SNID', default_snname)

    # add a zeropoint column
    outdata = lcdata['OBS']
    nobs = len(outdata)
    outdata['zptmag'] = np.repeat(27.5, nobs)

    # some of the columns are Null/zeroes - don't keep them around needlessly
    bad_cols = ['MAG','MAGERR', 'FIELD']
    for bad_col in bad_cols:
        if bad_col in outdata.dtype.names:
            outdata.remove_column(bad_col)

    # make column names lowercase
    colnames = outdata.dtype.names
    for col in colnames:
        outdata.rename_column(col, col.lower())

    if sampname is not None:
        meta["SOURCE"] = sampname
        if sampname == 'Foundation_DR1':
            outdata['flt'] = np.array([x +'_PS1' for x in outdata['flt']])

    # get additional meta if possible - this will overwrite any existing key if it has the same name
    if metafile is not None:
        if os.path.exists(metafile):
            allmeta = at.Table.read(metafile, format="ascii")
            if "SNID" in allmeta.colnames and snname in allmeta["SNID"]:
                thismeta = allmeta[allmeta["SNID"] == snname][0]
                for col in thismeta.colnames:
                    if col != "SNID":
                        meta[col] = thismeta[col]
            else:
                print("WARNING: SN {} could not be succesfully cross-matched to metadata file. If this is unexpected, please check the provided file contains a SNID column with an entry for {}.".format(snname, snname))
        else:
            print("WARNING: Metadata file cannot be found - continuing without.")

    # save the meta and data together
    outdata.meta = meta
    return snname, outdata

def write_snana_lcfile(output_dir, snname, mjd, flt, mag, magerr, tmax, z_helio, z_cmb, z_cmb_err, ebv_mw, ra=None, dec=None, author="anonymous", survey=None, paper=None, filename=None):
    '''
    Write user data to an SNANA-like light curve file

    Parameters
    ----------
    output_dir : str
        Path to a directory where the file will be written. A default filename
        will be used, but you can specify your own with the `filename` argument.
        Default name format is `snname[_survey][_paper].snana.dat`, with the
        survey and/or paper being appended to the name if provided.
    snname : str
        Name of the supernova
    mjd : list or :py:class:`numpy.array`
        Modified Julian Dates of observations
    flt : list or :py:class:`numpy.array` of str
        Filter idenitifiers of observations
    mag : list or :py:class:`numpy.array`
        Magnitudes of observations
    magerr : list or :py:class:`numpy.array`
        Magnitude errors of observations
    tmax : float
        Estimated time of maximum
    z_helio : float
        Heliocentric redshift
    z_cmb : float
        CMB-frame redshift
    z_cmb_err : float
        Error on CMB-frame redshift (excluding peculiar velocity uncertainty contribution)
    ebv_mv : float
        E(B-V) reddening due to the Milky Way
    ra : float, optional
        Right Ascension, to be writen to the header if desired
    dec :  float, optional
        Declination, to be written into the header if desired
    author : str, optional
        Who is creating this file? Will be printed into the header's
        preamble, if desired
    survey : str, optional
        Optional argumanet specifying the survey the data came from. Will be
        written into the header and filename if provided.
    paper : str, optional
        Optional argument specifying the paper the data came from. Will be
        written into the filename if provided.
    filename : str, optional
        Custom filename to save as within `output_dir`. If not provided,
        a default format will be used. Do not provide an extension, as
        this will be added automatically.

    Returns
    -------
    path : str
        Full path to the generated light curve file.

    Notes
    -----
    This will write a user's data to the SNANA-like file format readable by
    out I/O routines. It will write the provided metadata into the file
    header, so this will be read in and used correctly by BayeSN. All vital
    metadata are required as inputs to this function.
    '''
    if not (len(mjd) == len(flt) == len(mag) == len(magerr)):
        raise ValueError("Provided columns are not the same length!")

    if not os.path.exists(output_dir):
        raise ValueError("Requested output directory does not exist!")

    tab = at.Table([mjd, flt, mag, magerr], names=["MJD", "FLT", "MAG", "MAGERR"])
    #Compute fluxcal and fluxcalerr
    tab["FLUXCAL"] = 10**((27.5 - tab["MAG"])/2.5)
    tab["FLUXCALERR"] = tab["FLUXCAL"]*tab["MAGERR"]*np.log(10)/2.5
    #Column which designates observations
    tab["VARLIST:"] = ["OBS:"]*len(tab)
    #Round fluxes and flux errors
    tab["FLUXCAL"] = np.round(tab["FLUXCAL"],4)
    tab["FLUXCALERR"] = np.round(tab["FLUXCALERR"],4)
    #Reorder columns
    tab = tab["VARLIST:", "MJD", "FLT", "FLUXCAL", "FLUXCALERR", "MAG", "MAGERR"]

    #Divider for the header
    divider = "-"*59

    #Write a preamble to the metadata dictionary
    datestamp = time.strftime("%Y.%m.%d",time.localtime())
    timestamp = time.strftime("%H.%M hrs (%Z)",time.localtime())
    preamble = ("\n# SNANA-like file generated from user-provided data\n" +
        "# Zeropoint of the converted SNANA file: 27.5 mag\n" +
        "# {}\n".format(divider) +
        "# Data table created by: {}\n".format(author) +
        "# On date: {} (yyyy.mm.dd); {}.\n".format(datestamp, timestamp) +
        "# Script used: BayeSNmodel.io.write_snana_lcfile.py\n" +
        "# {}".format(divider))
    tab.meta = {"# {}".format(snname): preamble}

    #Add metadata
    tab.meta["SNID:"] = snname
    if survey is not None:
        tab.meta["SOURCE:"] = survey
    if ra is not None:
        tab.meta["RA:"] = ra
    if dec is not None:
        tab.meta["DEC:"] = dec
    filters = ",".join(at.unique(tab, keys="FLT")["FLT"])
    tab.meta.update({"MWEBV:": ebv_mw, "REDSHIFT_HELIO:": z_helio, "REDSHIFT_CMB:": z_cmb, "REDSHIFT_CMB_ERR:": z_cmb_err, "PEAKMJD:": tmax, "FILTERS:": filters, "#": divider, "NOBS:": len(tab), "NVAR:": 6})

    #Write to file
    if filename is None:
        filename = snname + (survey is not None)*"_{}".format(survey) + (paper is not None)*"_{}".format(paper) + ".snana.dat"
    sncosmo.write_lc(tab, output_dir + filename, fmt="salt2", metachar="")

    #Write terminating line
    with open(output_dir + filename, "a") as f:
        f.write("END:")

    #Completion message
    print("File written to {}".format(output_dir + filename))

    #Return filename
    return output_dir + filename


def sort_lcs(lcs):
    # sort the objects
    outlc_names =  sorted(lcs.keys())
    outlcs = OrderedDict()
    for key in outlc_names:
        outlcs[key] = lcs[key]
    return outlcs


def read_sn_sample_file(sampfile='M20_training_set', metafile='M20_training_set_meta', input_dir="bayesn-data/lcs"):
    '''
    Read light curves pre-defined in a sample file

    Parameters
    ----------
    sampfile : str
        Specify sample of SNe to read. Sample files based on the Mandel+20
        and Thorp+21 training sets are provided in bayesn-data/lcs/tables.
        Unlike `read_sn_sample` the sample file defines the light curve for
        each object whereas the other relies of sample defined by directory
        therefore each SN can be from a different survey and each SN
        can potentially have multiple files defined, and these will be
        concatenated
    metafile : str, optional
        A file which may contain additional metadata about the sample.
        Can be either a full path, or the name of a file in
        input_dir/meta.
        This file must have an `SNID` column. SNe in the sample will be
        cross-matched against the `SNID` column in `metafile`.
        Data in the `metafile` will be added to
        the metadata of the returned astropy Tables. Any info found
        in `metafile` will take precedence over metadata in the
        SNANA files. This mechanism can be used to provide
        custom/improved metadata.
    input_dir : str, optional
        Path to a directory containing the sampfile, metafile and
        light curves. Defaults to bayesn-data/lcs

    Returns
    -------
    lcs : dict
        ``lcs`` is a dictionary of light curves
        for the sample specified in ``sampfile``
        ``lcs`` is indexed by SNID of the objects
        each element of ``lcs`` is a astropy Table
        ``lcs[snid].meta`` has the meta information for each object

    '''
    # we'll try to load ``sampfile`` from input_dir/tables/
    # if that doesn't work, we'll just try treating ``sampfile`` as a path
    if os.path.exists("{}/tables/{}.txt".format(input_dir, sampfile)):
        lcsamp = at.Table.read("{}/tables/{}.txt".format(input_dir, sampfile), comment='#', names=('sn','source','lcfile'), format='ascii.no_header')
    elif os.path.exists("{}/tables/{}".format(input_dir, sampfile)):
        lcsamp = at.Table.read("{}/tables/{}".format(input_dir, sampfile), comment='#', names=('sn','source','lcfile'), format='ascii.no_header')
    elif os.path.exists(sampfile):
        lcsamp = at.Table.read(sampfile, comment='#', names=('sn','source','lcfile'), format='ascii.no_header')
    else:
        raise OSError("Could not find sample file {}".format(sampfile))

    # check the metafile can be found, and warn once if it can't be
    # we'll look in input_dir/meta/ first, and try as a path
    # if that doesn't work
    if metafile is not None:
        if os.path.exists("{}/meta/{}.txt".format(input_dir, metafile)):
            metafile = "{}/meta/{}.txt".format(input_dir, metafile)
        elif os.path.exists("{}/meta/{}".format(input_dir, metafile)):
            metafile = "{}/meta/{}".format(input_dir, metafile)
        elif os.path.exists(metafile):
            pass
        else:
            print("WARNING: Metadata file cannot be found - continuing without.")
            metafile = None

    # load the data
    lcs = {}
    for entry in lcsamp:

        candfiles = entry['lcfile']
        lcfiles = candfiles.split(',')
        source = entry['source']
        for thisfile in lcfiles:
            if os.path.exists(thisfile):
                lcfile = thisfile
            elif os.path.exists(os.path.join(input_dir, source, thisfile)):
                lcfile = os.path.join(input_dir,source, thisfile)
            else:
                raise IOError("Couldn't find file at {} or {}!".format(thisfile, os.path.join(input_dir,source, thisfile)))

            snname, outdata = read_snana_lcfile(lcfile, sampname=source, metafile=metafile)

            if snname is None:
                continue

            # check if we've already loaded photometry from this sample for this object
            # in this case, we concatenate
            if snname in lcs:
                sndata = lcs[snname]
                outdata = at.vstack((sndata, outdata), metadata_conflicts='silent')
                outdata.meta['NOBS'] = len(outdata)
            else:
                pass

            lcs[snname] =  outdata
    return sort_lcs(lcs)


def read_sn_sample(sampname, input_dir='bayesn-data/lcs', extension='dat', sourcename=None, metafile=None):
    '''
    Read light curves for a specific SN sample

    Parameters
    ----------
    input_dir : None or str
        input directory containing the light curves. Defaults to ``bayesn-data/lcs``
    sampname : str
        Specify sample name of SNe to read (CfA/CSP/LCO/misc).
        This should be the name of a directory under ``input_dir``
    sourcename : None or str
        Specify a source name of SNe to read within the sample
        By default this is just based on ``sampname`` and interpreted as ``*sampname*``
        If specified, it is used instead of the sampname-based value as ``*sourcename*``
        This provides fine-grained control to read only objects from one paper
    metafile : str, optional
        A file which may contain additional metadata about the sample.
        Can be either a full path, or the name of a file in
        input_dir/meta.
        This file must have an `SNID` column. SNe in the sample will be
        cross-matched against the `SNID` column in `metafile`.
        Data in the `metafile` will be added to
        the metadata of the returned astropy Tables. Any info found
        in `metafile` will take precedence over metadata in the
        SNANA files. This mechanism can be used to provide
        custom/improved metadata.

    Returns
    -------
    lcs : dict
        ``lcs`` is a dictionary of light curves found in ``input_dir``
        for sample ``sampname``
        ``lcs`` is indexed by SNID of the objects
        each element of ``lcs`` is a astropy Table
        ``lcs[snid].meta`` has the meta information for each object

    Notes
    -----
        ``lcs`` cannot have duplicates - these are skipped with a warning
        Any failure to read a file will also result in a warning

    See also
    --------
        ``read_snana_lcfile`` which reads individual files

    '''

    # check that the input dir exists
    if not os.path.isdir(input_dir):
        message = '{} is not a directory'.format(input_dir)
        raise ValueError(message)

    # check the metafile can be found, and warn once if it can't be
    # we'll look in input_dir/meta/ first, and try as a path
    # if that doesn't work
    if metafile is not None:
        if os.path.exists("{}/meta/{}.txt".format(input_dir, metafile)):
            metafile = "{}/meta/{}.txt".format(input_dir, metafile)
        elif os.path.exists("{}/meta/{}".format(input_dir, metafile)):
            metafile = "{}/meta/{}".format(input_dir, metafile)
        elif os.path.exists(metafile):
            pass
        else:
            print("WARNING: Metadata file cannot be found - continuing without.")
            metafile = None

    # get a list of the lcfiles in this sample from input_dir
    sampname = str(sampname)

    if sourcename is None:
        if sampname == 'misc':
            sourcename = '*'
        else:
            sourcename = '*{}*'.format(sampname)
    else:
        sourcename = '*{}*'.format(str(sourcename))

    basepattern = '{}{}'.format(sourcename, extension)
    filepattern = os.path.join(input_dir, sampname, basepattern)
    lcfiles = glob.glob(filepattern)

    # load the data
    lcs = {}
    for lcfile in lcfiles:

        snname, outdata = read_snana_lcfile(lcfile, sampname=sampname, metafile=metafile)

        if snname is None:
            continue

        # check if we've already loaded photometry from this sample for this object
        if snname in lcs:
            message = 'Object {} from sample {} already exists. Skipping {}'.format(snname, sampname, lcfile)
            warnings.warn(message)
            continue

        lcs[snname] =  outdata

    return sort_lcs(lcs)
