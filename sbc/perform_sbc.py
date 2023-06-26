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

import sys
sys.path.append('model_files/')
from SBC import *
import argparse,yaml

parser = argparse.ArgumentParser(description="SBC Input Dust Hyperparameters")
parser.add_argument("--tauA",default=0.5,help='Float for AV Dust Hyp')
parser.add_argument("--muRV",default=2.5,help='Float for RV Pop. Mean')
parser.add_argument("--sigRV",default=0.5,help='Float for RV Gaussian Dist.')

args = parser.parse_args().__dict__
dust_hyps = {x:float(args[x]) for x in ['tauA','muRV','sigRV']}
print (f"Performing SBC using input dust hyperparameters: {dust_hyps}")


if __name__ == "__main__":
	with open('sbc.yaml') as f:
		sbc_choices = yaml.load(f, Loader=yaml.FullLoader)

	#Draw dust hyperparameters from simulation distributions; create 1000 datasets
	edit_dict   = {'simulate_parameters':{**{'Nsims':200},**dust_hyps}}

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
	print ('Begin Fitting SN Datasets')
	sbc.fit_truths()
