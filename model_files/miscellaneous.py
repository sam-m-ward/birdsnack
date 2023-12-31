"""
Miscellaneous Module

Module contains simple functions that are called upon by various scripts

Contains
--------------------
ensure_folders_to_file_exist(savename)
trim_to_common_SNS(DF_M, verbose=True, drop_sns=None)
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""
import os
import numpy as np
from contextlib import suppress

def ensure_folders_to_file_exist(savename):
    """
    Ensure folders to file exist

    Makes sure all the directories exist up to where savefile will be saved as savename='path/to/dir/savefile.ext'

    Parameters
    ----------
    savename: str
        savename example is path/to/dir/savefile.ext

    Returns
    ----------
        Checks all subdirectories exist, and if they don't, create those directories
    """
    dirs = savename.split('/')[:-1]
    successive_dirs = ['/'.join(dirs[:id+1])+'/' for id in range(len(dirs))]
    for dir in successive_dirs:
        if not os.path.exists(dir):
            os.mkdir(dir)

def trim_to_common_SNS(DF_M, verbose=True, drop_sns=None):
    """
    Trim to Common SNs

    Takes in DF_M, which has multiple key,value entries, some of which may have SNe dropped
    Pools info together and ensures all values have the same SNe in them

    Parameters
    ----------
    DF_M : dict
        keys are tilist, extra, Tmax, choices

    verbose : bool (optional; default=True)
        if True, print new sample size

    drop_sns : list of str (optional; default = None)
        list of specific SNs to drop

    Returns
    ----------
    DF_M where all entries have same common SNS
    """
    columns = [c for c in DF_M if c not in ['extra','Tmax','choices']]

    if drop_sns is not None:
        DF_M[columns[0]].drop(drop_sns, axis=0, inplace=True)

    common_SNS = list(DF_M[columns[0]].index)
    for ti in columns:
        common_SNS = [sn for sn in common_SNS if sn in list(DF_M[ti].index)]
    for ti in DF_M['extra']:
        common_SNS = [sn for sn in common_SNS if sn in list(DF_M['extra'][ti].index)]

    ordered_SNS = []
    for sn in common_SNS:
        if ordered_SNS.count(sn)==0: ordered_SNS.append(sn)

    common_SNS  = [sn for sn in ordered_SNS if sn in list(set(common_SNS))]

    for ti in columns:
        DF_M[ti] = DF_M[ti].loc[common_SNS]
    for ti in DF_M['extra']:
        DF_M['extra'][ti] = DF_M['extra'][ti].loc[common_SNS]
    DF_M['Tmax'] = DF_M['Tmax'].loc[common_SNS]

    if verbose: print (f'New sample size in DF_M is {len(common_SNS)}')

    return DF_M
