#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from load_raw_data import *
from birdsnack_model import BIRDSNACK

#Load up light curves of fiducial sample of SNe, and metadata
dataloader = LOAD_DATA()
SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')
dfmeta           = dataloader.dfmeta
edit_dict        = {}
###################
##u-band Analysis
#edit_dict = {  'preproc_parameters'    :   {'DF_savekey':'uBVriJH','pblist':[s for s in 'uBVriJH']},
#                'analysis_parameters'   :   {'HBM_savekey':'uBVriJH','lam_choice':'central'}}
#edit_dict = {  'preproc_parameters'    :   {'DF_savekey':'uBVriJH','pblist':[s for s in 'uBVriJH']},
#                'analysis_parameters'   :   {'HBM_savekey':'uBVriJH_Cens1.0','lam_choice':'central','CensoredData':True,'CensoredCut':1.0}}
#edit_dict = {  'preproc_parameters'    :   {'DF_savekey':'uBVriJH','pblist':[s for s in 'uBVriJH']},
#                'analysis_parameters'   :   {'HBM_savekey':'uBVriJH_LowRVs_Cens1.0','lam_choice':'central','CensoredData':True,'CensoredCut':1.0},
#                'additional_cut_parameters' : {'extra_drop_SNe':{sn:'Low RVs' for sn in ['2009ds','16abc']} } }
#edit_dict = {  'preproc_parameters'    :   {'DF_savekey':'uBVriJH','pblist':[s for s in 'uBVriJH']},
#                'analysis_parameters'   :   {'HBM_savekey':'uBVriJH_Phasemax4_LowRVs_Cens1.0','lam_choice':'central','CensoredData':True,'CensoredCut':1.0},
#                'additional_cut_parameters' : {'extra_drop_SNe':{sn:'Low RVs' for sn in ['2009ds','16abc']} ,'phase_max':4} }
###################
##CSP-Only Analysis
#SNSsnpy_fiducial = dataloader.load_SNSsnpy('snpy_SNS_CSP.pkl')
#edit_dict = {  'preproc_parameters'    :   {'DF_savekey':'CSPOnly'},
#                'analysis_parameters'   :   {'HBM_savekey':'CSPOnly_Cens1.0','CensoredData':True,'CensoredCut':1.0}}
###################
#Skewed RV
#edit_dict = { 'analysis_parameters'   :   {'HBM_savekey':'SkewRV_Cens1.0','CensoredData':True,'CensoredCut':1.0,'skew_RV':True,'n_sampling':6000}}
#edit_dict = { 'analysis_parameters'   :   {'HBM_savekey':'SkewRV_LowRVs_Cens1.0','CensoredData':True,'CensoredCut':1.0,'skew_RV':True,'n_sampling':12000},
#              'additional_cut_parameters' : {'extra_drop_SNe':{sn:'Low RVs' for sn in ['2009ds','16abc']}}}
###################
#Student-T RV
#edit_dict = { 'analysis_parameters'   :   {'HBM_savekey':'RVstudT_Cens1.0','CensoredData':True,'CensoredCut':1.0,'RVprior':'StudentT','n_sampling':12000}}
#edit_dict = { 'analysis_parameters'   :   {'HBM_savekey':'RVstudT_LowRVs_Cens1.0','CensoredData':True,'CensoredCut':1.0,'RVprior':'StudentT','n_sampling':12000},
#              'additional_cut_parameters' : {'extra_drop_SNe':{sn:'Low RVs' for sn in ['2009ds','16abc']}}}
###################
#Skewed Intrinsic
#edit_dict = { 'analysis_parameters'   :   {'HBM_savekey':'SkewInt_Cens1.0','CensoredData':True,'CensoredCut':1.0,'skew_int':True,'n_sampling':12000}}
edit_dict = { 'analysis_parameters'   :   {'HBM_savekey':'SkewInt_LowRVs_Cens1.0','CensoredData':True,'CensoredCut':1.0,'skew_int':True,'n_sampling':12000},
              'additional_cut_parameters' : {'extra_drop_SNe':{sn:'Low RVs' for sn in ['2009ds','16abc']}}}
###################
#Skewed RV and Intrinsic
#edit_dict = { 'analysis_parameters'   :   {'HBM_savekey':'SkewRV_SkewInt_Cens1.0','CensoredData':True,'CensoredCut':1.0,'skew_RV':True,'skew_int':True,'n_sampling':12000}}
###################
#Load into Bird-Snack
bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dfmeta, edit_dict=edit_dict)

#Trim Sample [Turn on for u-band or CSP-only analysis]
#bs.trim_sample()

#Get peak magnitudes
bs.get_peak_mags()

#Perform additional sample cuts
bs.additional_cuts()

#Plot Mag Deviations
#bs.plot_mag_deviations()

#Plot Colours
#bs.plot_colour_corner()

#Plot LCs
#bs.plot_lcs()

#Fit HBM to data
bs.fit_stan_model()#Rhat_threshold=1.05)

#Plot posterior samples
bs.plot_posterior_samples()
#summary = bs.plot_posterior_samples(returner=True)
#line = [summary['NSNe(NCens)']]+[value for key,value in summary.items() if key!='NSNe(NCens)']
#line = ' & '.join(line) + ' \\\\ '
#print (line)
