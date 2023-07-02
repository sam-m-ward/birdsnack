#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *
from CYCLE_DICTS import *
from birdsnack_model import BIRDSNACK, get_edit_dict
import argparse, shutil
parser = argparse.ArgumentParser(description="BirdSnack Analysis Cycle Script")
parser.add_argument("--mode",default='Science',help='Options are Science,CensoredData')
MODE   = parser.parse_args().mode
CYCLE_DICT = dict(zip(['Science','CensoredData'],[CYCLE_DICT_Science,CYCLE_DICT_CensoredData]))[MODE]

###############
#FOR LOW RVs Analysis
CYCLE_DICT['COMMON_CHANGES']['newdict'] = {**CYCLE_DICT['COMMON_CHANGES']['newdict'],**{'extra_drop_SNe':{sn:'Low RVs' for sn in ['2009ds','16abc']}}}
CYCLE_DICT['COMMON_CHANGES']['HBMappender'] = 'LowRVs_Cens1.0'
###############

#Load up light curves of fiducial sample of SNe, and metadata
dataloader = LOAD_DATA()
SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')
dfmeta           = dataloader.dfmeta

#Settings for deleting Stan build periodically to conserve laptop memory
#periodically_delete = False
periodically_delete = 4
stan_build_dir = '/Users/samward/Library/Caches/httpstan/'

Summary_Strs = {} ; ISIM = 0
for HBM_savekey in CYCLE_DICT['RUNS']:
    edit_dict = get_edit_dict(dataloader.choices,CYCLE_DICT,HBM_savekey)
    print ('###'*10)
    print (f"MODE{MODE}; Performing analysis on {edit_dict['analysis_parameters']['HBM_savekey']}")
    #Load into Bird-Snack
    bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dfmeta, edit_dict=edit_dict)
    #Get peak magnitudes
    bs.get_peak_mags()
    #Perform additional sample cuts; set cutter=False (so SNe with magerrs>0.3mag are retained, ensures common sample fitted for consistency across pre-proc analysis variants)
    bs.additional_cuts(cutter=False)
    #Fit HBM to data
    if not os.path.exists(f"{bs.FITSpath}FIT{bs.choices['analysis_parameters']['HBM_savekey']}.pkl"):
        ISIM += 1
        bs.fit_stan_model()
        if periodically_delete:
            if ISIM%periodically_delete==0:
                delete_tmp_files = glob(f'{stan_build_dir}*')#Need to do this manually otherwise files stack up in tmp/ and for loop dies
                for file in delete_tmp_files:
                    shutil.rmtree(file)
    #Plot posterior samples, and get posterior summary string
    Summary_Str = bs.plot_posterior_samples(returner=True)
    Summary_Strs[CYCLE_DICT['RUNS'][HBM_savekey]['label']] = Summary_Str

#Print Results
for label, summary in Summary_Strs.items():
    line = [label]+[summary['NSNe(NCens)']]+[value for key,value in summary.items() if key!='NSNe(NCens)']
    line = ' & '.join(line) + ' \\\\ '
    print (line)
