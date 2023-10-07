"""
Perform Simulation-Based Calibration
--------------------
Simple script that uses SBC class to:
	simulate a bunch of SN datasets
	Fit them using Bird-Snack model
	Assess Recovery of Input (Hyper)parameters

Functions:
	update_edit_dict_for_ppc(sbc_choices, edit_dict)
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

import sys
sys.path.append('model_files/')
from SBC import *
import yaml
from sbc_plot_functions import update_edit_dict_for_ppc

additional_pars = ['AVs','RVs']#Save these as well as dust hyperparameters, optional

#Directory to periodically_delete stan build, conserve memory
stan_build_dir      = '/Users/samward/Library/Caches/httpstan/'
periodically_delete = 3

#Choose to fit with Exp or Gamma
BIRDSNACK_EDIT_DICT = {'analysis_parameters':{'HBM_savekey':'PPC_FitdoubleAVdist_wExp','AVprior':'Exp','n_warmup':1000,'n_sampling':1000,'n_thin':1000}}
#BIRDSNACK_EDIT_DICT = {'analysis_parameters':{'HBM_savekey':'PPC_FitdoubleAVdist_wGamma','AVprior':'Gamma','n_warmup':1000,'n_sampling':3000,'n_thin':1000}}

#Simulate x4 High Red SNe
edit_dict1 = {'simulate_parameters':{'Nsims':100,'pre_defined_hyps':{'load_file':'HighRedSNe_AVGauss_FixedInt','load_int_file':'AVExp_BVcut1.0','load_ext_file':'HighRedSNe_AVGauss_FixedInt'}}}
with open('ppc.yaml') as f:
	sbc_choices = yaml.load(f, Loader=yaml.FullLoader)
edit_dict1 = update_edit_dict_for_ppc(sbc_choices,edit_dict1)
sbc1 = SBC_CLASS(sbc_choices,edit_dict1)
sbc1.simulate_truths()
TRUTHS_DICT1 = sbc1.get_truths()

print ('###'*20)
print ('Done simulating x4 high-red SNe, now simulate x65 low-to-moderate-reddening SNe')
print ('###'*20)

#Simulate x65 Low-to-moderate Reddening SNe
edit_dict2 = {'simulate_parameters':{'Nsims':100,'pre_defined_hyps':{'load_file':'AVExp_BVcut1.0'}}}
with open('ppc.yaml') as f:
	sbc_choices = yaml.load(f, Loader=yaml.FullLoader)
edit_dict2 = update_edit_dict_for_ppc(sbc_choices,edit_dict2)
sbc2 = SBC_CLASS(sbc_choices,edit_dict2)
sbc2.simulate_truths()
TRUTHS_DICT2 = sbc2.get_truths()

#Add on the measurements of the 4 high reddening SNe to the remaining x65
for i,truths in TRUTHS_DICT2.items():
	truths_highred  = TRUTHS_DICT1[i]
	truths.mobs     = list(np.concatenate((truths.mobs,truths_highred.mobs)))
	truths.errors   = list(np.concatenate((truths.errors,truths_highred.errors)))
	TRUTHS_DICT2[i] = truths

#Fit combined dataset of x69 SNe
sbc2.TRUTHS_DICT = TRUTHS_DICT2
sbc2.fit_truths(edit_dict=BIRDSNACK_EDIT_DICT,periodically_delete=periodically_delete,outputdir=stan_build_dir,additional_pars=additional_pars)
