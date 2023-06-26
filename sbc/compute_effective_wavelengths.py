"""
Compute Effective Wavelengths
--------------------
Simple script that uses SBC class to:
	simulate a bunch of SN datasets
	get effective wavelengths
	save as DLAM_NEW.txt
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk
"""

import sys
sys.path.append('model_files/')
from SBC import *
import yaml


if __name__ == "__main__":
	with open('sbc.yaml') as f:
		sbc_choices = yaml.load(f, Loader=yaml.FullLoader)

	#Draw dust hyperparameters from simulation distributions; create 1000 datasets
	edit_dict   = {'simulate_parameters':{'Nsims':1000,'tauA':None,'muRV':None,'sigRV':None}}

	#Get SBC_CLASS
	sbc = SBC_CLASS(sbc_choices,edit_dict)

	#Simulate SNe Datasets
	print ('Simulating SNe Datasets')
	sbc.simulate_truths()

	#Load up the simulated datasets
	print ('Loading up SNe Datasets')
	TRUTHS_DICT = sbc.get_truths()

	#Compute lam_effs; simply loops through each SN
	print ('Computing effective wavelengths')
	df = sbc.get_leffs_df()

	#Print Effective Wavelengths
	print ('Effective Wavelengths are:')
	lam_eff = df.copy().median(axis=0)
	print (lam_eff)

	#Save list of lam_effs
	print (f"Saving as DLAM_NEW.txt at {sbc.rootpath}DLAMS/DLAM_NEW.txt")
	df.to_csv(f"{sbc.rootpath}DLAMS/DLAM_NEW.txt",index=False)
