#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *
dataloader = LOAD_DATA()

#'''Get Sample Split of Surveys
with open(f'{dataloader.SNSpath}SNSsnpy_fiducial.pkl','rb') as f:
    SNSsnpy_fiducial = pickle.load(f)

surveys = {}
for sn in SNSsnpy_fiducial:
    survey = SNSsnpy_fiducial[sn]['survey']
    if survey in surveys:
        surveys[survey].append(sn)
    else:
        surveys[survey] = [sn]

print ('###'*5)
for survey in surveys:
    print (survey, len(surveys[survey]))
print ('###'*5)
#'''

#'''#Get Breakdown of reasons for why SNe were cut from sample
import copy
def get_sn_reasons(sn,reasons):
    sn_reasons = []
    for reason in reasons:
        if sn in reasons[reason]:
            sn_reasons.append(reason)
    return sn_reasons


with open(f'{dataloader.SNSpath}SURVEYS.pkl','rb') as f:
    SURVEYS = pickle.load(f)

Unique_sne = []
for survey in SURVEYS:
    combined_lcs = {**SURVEYS[survey]['trimmed_lcs'],**SURVEYS[survey]['retained_lcs']}
    Unique_sne.extend(list(combined_lcs.keys()))
Unique_sne = list(set(Unique_sne))
Unique_sne.sort()
print (len(Unique_sne),"Remember also that 14yw doesnt load so add +1; 09dc appears here (it is just the CfA file that is ignored, so no additional SN)")
print (f"2009dc in Unique_sne: {'2009dc' in Unique_sne}")
print (f"Therefore total unique SNe is {len(Unique_sne)+1}")

Pecking_order_of_reasons = [
['Has_pblist'],
['snpy_cut'],
['Points_Before_After_RefBand_Max','None_Tmax_GP_restframe','Large_Tmax_GP_restframe_std'],
['B_0','V_0','r_0','i_0'],
['J_0','H_0'],
]

def get_Number_of_SNe_dropped():
    Remaining_SNe = copy.deepcopy(Unique_sne)
    for set_of_reasons in Pecking_order_of_reasons:
        print ('***'*10)
        Ndropped_at_this_stage = 0 ; sne_to_remove = {}
        for sn in Remaining_SNe:
            #print ('###'*5)
            keep = False ; surveys_this_sn_appeared_in = []
            for survey in dataloader.choices['load_data_parameters']['Pecking_Order']:
                combined_lcs = {**SURVEYS[survey]['trimmed_lcs'],**SURVEYS[survey]['retained_lcs']}
                if sn in list(combined_lcs.keys()):
                    surveys_this_sn_appeared_in.append(survey)
                    if sn in list(SURVEYS[survey]['retained_lcs'].keys()):
                        keep = True
            #print (f'SN {sn} appeared in {surveys_this_sn_appeared_in} surveys')
            if keep:
                #print (f'SN {sn} was kept in survey {surveys_this_sn_appeared_in[0]}')
                pass
            if not keep:
                lowest_survey_in_pecking_order = surveys_this_sn_appeared_in[-1]
                reasons_for_dropping = get_sn_reasons(sn,SURVEYS[lowest_survey_in_pecking_order]['reasons'])
                #print (f'SN {sn} dropped because in survey {lowest_survey_in_pecking_order}: {reasons_for_dropping}')
                for rr in set_of_reasons:
                    if rr in reasons_for_dropping:
                        sne_to_remove[sn] = lowest_survey_in_pecking_order
                        Ndropped_at_this_stage += 1
                        break

        Remaining_SNe = [sn for sn in Remaining_SNe if sn not in sne_to_remove]
        if set_of_reasons==['Points_Before_After_RefBand_Max','None_Tmax_GP_restframe','Large_Tmax_GP_restframe_std']:
            #print (f'For {set_of_reasons}, SNe and survey removed are')
            #print (sne_to_remove)
            TO_INSPECT_VISUALLY = sne_to_remove
        print ('~~~'*5)
        print (f'For {set_of_reasons} there are {Ndropped_at_this_stage} removed')
    return TO_INSPECT_VISUALLY

TO_INSPECT_VISUALLY = get_Number_of_SNe_dropped()
#'''

'''#Inspect the B-band-cut SNe visually
SNSsnpy_full = {}
for path_file_survey in dataloader.choices['load_data_parameters']['load_path_file_survey']:
    with open(f"{dataloader.SNSpath}{path_file_survey[1]}",'rb') as f:
            SNSsnpy_survey = pickle.load(f)
    SNSsnpy_full[path_file_survey[-1]] = SNSsnpy_survey

SNSsnpy_inspect = {}
for sn in TO_INSPECT_VISUALLY:
    survey = TO_INSPECT_VISUALLY[sn]
    SNSsnpy_inspect[sn] = SNSsnpy_full[survey][sn]

from birdsnack_model import BIRDSNACK
bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_inspect}, configname='loader_config.yaml', edit_dict={'plotting_parameters':{'show':True}}, dfmeta=dataloader.dfmeta)
bs.plot_lcs()
'''


'''#Show final trims to fiducial sample
with open(f"{dataloader.SNSpath}SNSsnpy_combined.pkl",'rb') as f:
        SNSsnpy_combined = pickle.load(f)
from birdsnack_model import BIRDSNACK
bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_combined}, configname='loader_config.yaml', edit_dict={'load_data_parameters':{'load_mode':'preproc'}}, dfmeta=dataloader.dfmeta)
bs.get_peak_mags(savekey='combined')
bs.additional_cuts()
#'''
