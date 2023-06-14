"""
Miscellaneous Module

Module contains simple functions that are called upon by various scripts

Contains
--------------------
ensure_folders_to_file_exist(savename):
    simple function to create directories to where file will be saved if dirs do not already exist

--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""
import os
import numpy as np

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
