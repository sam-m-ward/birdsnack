#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *
#Load up and save raw data as snpyfile dictionaries, load up metadata
dataloader = LOAD_DATA(load_all_SNSsnpy=True)

from birdsnack_model import BIRDSNACK

#Try loading fiducial sample of SNSsnpy; fiducial dictionary of snpy files from multiple literature surveys
if not os.path.exists(f'{dataloader.SNSpath}SNSsnpy_fiducial.pkl'):
    #Try loading up SNS retained over multiple surveys, highest in pecking order
    if os.path.exists(f'{dataloader.SNSpath}SNSsnpy_combined.pkl'):
        with open(f'{dataloader.SNSpath}SNSsnpy_combined.pkl','rb') as f:
            SNSsnpy_combined = pickle.load(f)
    #If not available, compute this dictionary
    else:
        #Loop through surveys, load up snsnpy, apply corrections, convert to snana lc, estimate Tmax, trim on phase, check data availability, record whether trimmed or not
        SURVEYS = {}
        for path_file_survey in dataloader.choices['load_data_parameters']['load_path_file_survey']:
            print ('###'*3)
            bs = BIRDSNACK(loader={'path_file_survey':path_file_survey}, configname='loader_config.yaml', edit_dict={'load_data_parameters':{'load_mode':'preproc'}},dfmeta=dataloader.dfmeta)
            bs.trim_sample()
            SURVEYS[path_file_survey[-1]] = {'retained_lcs':bs.lcs,'trimmed_lcs':bs.trimmed_lcs, 'reasons':bs.reasons}
        #Use SURVEYS dictionary to get common sample of SNS, with those SNe appearing in multiple surveys selected according to Pecking order
        SNSsnpy_combined = dataloader.get_SNSsnpy_combined(SURVEYS)#Save SNSsnpy_combined.pkl


    bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_combined}, configname='loader_config.yaml', edit_dict={'load_data_parameters':{'load_mode':'preproc'}}, dfmeta=dataloader.dfmeta)
    bs.trim_sample()
    bs.get_peak_mags(savekey='combined')
    bs.additional_cuts()

    #Remaining SNe are fiducial sample
    fiducial_sns     = [sn for sn in list(bs.DF_M[0].index)]
    SNSsnpy_fiducial = {sn:SNSsnpy_combined[sn] for sn in fiducial_sns}
    with open(f'{dataloader.SNSpath}SNSsnpy_fiducial.pkl','wb') as f:
        pickle.dump(SNSsnpy_fiducial,f)

print ('Fiducial dictionary of SNe from multiple literature surveys compiled')
