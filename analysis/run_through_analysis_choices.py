#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *
from CYCLE_DICTS import *
from birdsnack_model import BIRDSNACK
import argparse
parser = argparse.ArgumentParser(description="BirdSnack Analysis Cycle Script")
parser.add_argument("--mode",default='Science',help='Options are Science,CensoredData')
MODE   = parser.parse_args().mode
CYCLE_DICT = dict(zip(['Science','CensoredData'],[CYCLE_DICT_Science,CYCLE_DICT_CensoredData]))[MODE]

#Load up light curves of fiducial sample of SNe, and metadata
dataloader = LOAD_DATA()
SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')
dfmeta           = dataloader.dfmeta

COMMON_CHANGES = CYCLE_DICT['COMMON_CHANGES']['newdict']
RUNS           = CYCLE_DICT['RUNS']

for HBM_savekey in RUNS:
    print ('###'*10)
    print (f'MODE{MODE}; Performing analysis on {HBM_savekey}')
    def get_edit_dict(choices,newdict):
        edit_dict = {glob_key:{} for glob_key in choices if glob_key!='rootpath'}
        for key in newdict:
            for glob_key in edit_dict:
                for kkey,vvalue in choices[glob_key].items():
                    if kkey==key:   edit_dict[glob_key][key] = newdict[key]
        return edit_dict

    edit_dict = get_edit_dict(dataloader.choices,{**RUNS[HBM_savekey]['newdict'],**COMMON_CHANGES})
    #Load into Bird-Snack
    bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dfmeta, edit_dict=edit_dict)
    #Get peak magnitudes
    bs.get_peak_mags()
    #Perform additional sample cuts
    bs.additional_cuts(cutter=False)
    #Fit HBM to data
    #bs.fit_stan_model()
