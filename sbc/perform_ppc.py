"""
Perform Simulation-Based Calibration
--------------------
Simple script that uses SBC class to:
	simulate a bunch of SN datasets
	Fit them using Bird-Snack model
	Assess Recovery of Input (Hyper)parameters
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""
def update_edit_dict_for_ppc(sbc_choices, edit_dict):
	"""
	Update edit_dictionary for Posterior Predictive Checks

	Parameters
	----------
	sbc_choices : dict
		dicionary of choices from .yaml regarding how datasets are simulated
	edit_dict : dict
		edits to above .yaml choices, suited for ppc

	Returns
	----------
	edit_dict : dict
		edited to have all important data to perform ppc
	"""
	#Open original BirdSnack fit to real data
	with open(f"{sbc_choices['load_parameters']['path_to_birdsnack_rootpath']}products/stan_fits/FITS/FIT{edit_dict['simulate_parameters']['pre_defined_hyps']['load_file']}.pkl",'rb') as f:
		FIT = pickle.load(f)
	#Get number of SNe (including censored SNe)
	S = FIT['stan_data']['S']
	#Get reference band (for building full vector mu_int where ref-band has element entry=0)
	zero_index      = FIT['choices']['analysis_parameters']['zero_index']
	#Whether using model for colours or deviations
	IntrinsicModel  = FIT['choices']['analysis_parameters']['IntrinsicModel']
	#Assign folder name using these values
	tauA,muRV,sigRV = FIT['df'].median(axis=0).round(4)[['tauA','mu_RV','sig_RV']].values
	dust_hyps = dict(zip(['tauA','muRV','sigRV'],[tauA,muRV,sigRV]))
	#Upload above to edit_dict
	edit_dict['simulate_parameters'] = {**edit_dict['simulate_parameters'],**dust_hyps,**{'S':S}}
	edit_dict['simulate_parameters']['pre_defined_hyps'] = {**edit_dict['simulate_parameters']['pre_defined_hyps'],**{'ppc_zero_index':zero_index,'ppc_IntrinsicModel':IntrinsicModel}}
	return edit_dict

import sys
sys.path.append('model_files/')
from SBC import *
import yaml

#Choices for applying HBM to simulated data
BIRDSNACK_EDIT_DICT = {'analysis_parameters':
						{'HBM_savekey':'PPC_CensoredCut1.0_ExpFitGamma',
						'CensoredData':True,'CensoredCut':1.0,
						'AVprior':'Gamma','n_warmup':1000,'n_sampling':1000}}

#Choices for simulating data based on previous stan fit with BirdSnack
edit_dict = {'simulate_parameters':{'Nsims':20,'pre_defined_hyps':{'load_file':'Fiducial_CensoredCut1.0'}}}

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
	sbc.fit_truths(edit_dict=BIRDSNACK_EDIT_DICT)
