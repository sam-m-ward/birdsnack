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
##########################################################################################
'''
#PPC Recovery of nuA=1; Applying Gamma dist. to fit Exp-simulated Data
BIRDSNACK_EDIT_DICT = {'analysis_parameters':
						{'HBM_savekey':'PPC_CensoredCut1.0_ExpFitGamma',
						'CensoredData':True,'CensoredCut':1.0,
						'AVprior':'Gamma','n_warmup':1000,'n_sampling':3000,'n_thin':1000}}

#Choices for simulating data based on previous stan fit with BirdSnack
edit_dict = {'simulate_parameters':{'Nsims':20,'S':250,'pre_defined_hyps':{'load_file':'AVExp_Cens1.0'}}}
#'''
##########################################################################################
#'''
#For CensoredData Simulations
BIRDSNACK_EDIT_DICT = {'analysis_parameters':
						{'HBM_savekey':'PPC_muRV2.5sigRV0.1tauA0.5_LowBVNoCens',
						'CensoredData':True,'CensoredCut':0.3,
						'AVprior':'Exp','n_warmup':1000,'n_sampling':2000,'n_thin':1000}}
#BIRDSNACK_EDIT_DICT = {'analysis_parameters':
#						{'HBM_savekey':'PPC_muRV2.5sigRV0.1tauA0.5_WithCensoredData',
#						'CensoredData':True,'CensoredCut':'inf',
#						'AVprior':'Exp','n_warmup':1000,'n_sampling':2000,'n_thin':1000}}
#edit_dict = {'simulate_parameters':{'Nsims':100,'S':100,'tauA':0.5,'muRV':2.5,'sigRV':0.5,'PredefinedExtrinsicHyps':False,
#			'pre_defined_hyps':{'load_file':'AVExp_Cens1.0'}}}
edit_dict = {'simulate_parameters':{'Nsims':100,'S':100,'tauA':0.5,'muRV':2.5,'sigRV':0.1,'PredefinedExtrinsicHyps':False,
			'pre_defined_hyps':{'load_file':'AVExp_Cens1.0'}}}
#'''
##########################################################################################
'''
#For X-H Sims (X-H model applied to Deviations-simulated SNe)
BIRDSNACK_EDIT_DICT = {'analysis_parameters':
						{'HBM_savekey':'PPC_CensoredCut1.0_DevFitXH',
						'CensoredData':True,'CensoredCut':1.0,
						'IntrinsicModel':'X-H','DataTransformation' : 'X-H','n_warmup':1000,'n_sampling':2000,'n_thin':1000}}

#Choices for simulating data based on previous stan fit with BirdSnack
edit_dict = {'simulate_parameters':{'Nsims':20,'pre_defined_hyps':{'load_file':'AVExp_Cens1.0'}}}
#edit_dict = {'simulate_parameters':{'Nsims':20,'S':250,'pre_defined_hyps':{'load_file':'AVExp_Cens1.0'}}}
#'''
##########################################################################################
#Directory to periodically_delete stan build, conserve memory
stan_build_dir      = '/Users/samward/Library/Caches/httpstan/'
periodically_delete = 3

if __name__ == "__main__":
	with open('ppc.yaml') as f:
		sbc_choices = yaml.load(f, Loader=yaml.FullLoader)

	#Get Edit Dictionary for Posterior Predictive Checks
	edit_dict = update_edit_dict_for_ppc(sbc_choices,edit_dict)

	#Get SBC_CLASS
	sbc = SBC_CLASS(sbc_choices,edit_dict)

	#Simulate SNe Datasets
	print ('Simulating SN Datasets')
	sbc.simulate_truths()

	#Load up the simulated datasets
	print ('Loading up SN Datasets')
	TRUTHS_DICT = sbc.get_truths()
	print ('####'*10)

	#Loop through Simulated Datasets and fit BirdSnack HBM
	#BIRDSNACK_EDIT_DICT implements e.g. fit with AVprior='Gamma'
	print ('Begin Fitting SN Datasets')
	sbc.fit_truths(edit_dict=BIRDSNACK_EDIT_DICT,periodically_delete=periodically_delete,outputdir=stan_build_dir)
