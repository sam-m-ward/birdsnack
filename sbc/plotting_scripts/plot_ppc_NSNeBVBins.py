import matplotlib.pyplot as pl
import numpy as np
import yaml

path_to_rootpath = '../'
path_to_birdsnack_rootpath = '../../'
import sys
sys.path.append(path_to_birdsnack_rootpath+'sbc/model_files/')
from SBC import *
from sbc_plot_functions import update_edit_dict_for_ppc

if __name__ == "__main__":
	Nsims = 100

	#'''
	FITS = [
	'AVExp_Full',
	'AVExp_CensInf',
	'AVExp_Cens1.0',
	'AVGamma_Full',
	'AVGamma_CensInf',
	'AVGamma_Cens1.0'
	]
	titles = [	'Exp. $A_V^s$ Fit to Full Samp. (69 SNe)',
				'Exp. $A_V^s$ Fit to Low $B-V$ Samp. (62 SNe) \nw/ Censored Data (7 SNe)',
				'Exp. $A_V^s$ Fit to Low $B-V$ Samp. (62 SNe) \nw/ Censored Data (3 SNe)',
				'Gamma $A_V^s$ Fit to Full Samp. (69 SNe)',
				'Gamma $A_V^s$ Fit to Low $B-V$ Samp. (62 SNe) \nw/ Censored Data (7 SNe)',
				'Gamma $A_V^s$ Fit to Low $B-V$ Samp. (62 SNe) \nw/ Censored Data (3 SNe)']
	#'''
	'''
	FITS = [
	'AVExp_Full_AVRVBeta',
	'AVExp_AVRVBeta_BVcut1.0_Cens',
	'AVExp_AVRVBeta_BVcut0.3_Cens',
	'AVGamma_Full_AVRVBeta',
	'AVGamma_AVRVBeta_BVcut1.0_Cens',
	'AVGamma_AVRVBeta_BVcut0.3_Cens',
	]
	titles = [	'Exp. $A_V^s$ Fit to Full Samp. (69 SNe) \nw/ $\\mathbf{\\beta}$-Model for $\\mu^s_{R_V}|A_V^s$',
				'Exp. $A_V^s$ Fit to $B-V<1.0\,$mag Samp. (65 SNe) \nw/ $\\mathbf{\\beta}$-Model for $\\mu^s_{R_V}|A_V^s$ & 4 Cens. SNe',
				'Exp. $A_V^s$ Fit to $B-V<0.3\,$mag Samp. (62 SNe) \nw/ $\\mathbf{\\beta}$-Model for $\\mu^s_{R_V}|A_V^s$ & 7 Cens. SNe',
				'Gamma $A_V^s$ Fit to Full Samp. (69 SNe) \nw/ $\\mathbf{\\beta}$-Model for $\\mu^s_{R_V}|A_V^s$',
				'Gamma $A_V^s$ Fit to $B-V<1.0\,$mag Samp. (65 SNe) \nw/ $\\mathbf{\\beta}$-Model for $\\mu^s_{R_V}|A_V^s$ & 4 Cens. SNe',
				'Gamma $A_V^s$ Fit to $B-V<0.3\,$mag Samp. (62 SNe) \nw/ $\\mathbf{\\beta}$-Model for $\\mu^s_{R_V}|A_V^s$ & 7 Cens. SNe',
				]
	#'''
	mapper = dict(zip(FITS,titles))

	load_samples_path = path_to_birdsnack_rootpath+'/products/stan_fits/FITS/'

	minBVs = []
	for file in FITS:
		#Get Posterior Median Hyperparameters from this file
		edit_dict = {   'load_parameters':{'path_to_rootpath':path_to_rootpath,'path_to_birdsnack_rootpath':path_to_birdsnack_rootpath},
						'simulate_parameters':{'Nsims':100,'pre_defined_hyps':{'load_file':file}}}
		with open(f"{path_to_rootpath}ppc.yaml") as f:
			sbc_choices = yaml.load(f, Loader=yaml.FullLoader)

		#Get Edit Dictionary for Posterior Predictive Simulations
		edit_dict = update_edit_dict_for_ppc(sbc_choices,edit_dict)

		#Get SBC_CLASS
		sbc = SBC_CLASS(sbc_choices,edit_dict)
		print (sbc.simfolder)
		#Simulate SNe Datasets
		print ('Simulating SN Datasets')
		sbc.simulate_truths()

		#Load up the simulated datasets
		print ('Loading up SN Datasets')
		TRUTHS_DICT = sbc.get_truths()

		Ncens = {'low':[],'mid':[],'high':[]}
		for ISIM, truths in TRUTHS_DICT.items():
			BVs = [truths.mexts[_][0]-truths.mexts[_][1] for _ in range(truths.S)]

			lowBVs  = [BV for BV in BVs if BV<=-0.3]
			midBVs  = [BV for BV in BVs if BV>=0.3 and BV<=1.0]
			highBVs = [BV for BV in BVs if BV>=1.0]

			Ncens['low'].append(len(lowBVs))
			Ncens['mid'].append(len(midBVs))
			Ncens['high'].append(len(highBVs))

			minBVs.append(min(BVs))

		FS=18
		from matplotlib.ticker import MultipleLocator
		pl.figure()
		pl.title(f'Predicted No. of SNe in High $B-V$ Bins from \n{mapper[file]}',fontsize=FS)
		pl.xlabel(r'$N_{SNe}$ (0.3$<B-V<$1.0 mag)',fontsize=FS)
		pl.ylabel(r'$N_{SNe}$ ($B-V>$1.0 mag)',fontsize=FS)
		try:
			if choices['CensoredData'] and choices['CensoredCut']==1.0:pl.plot(3,0,marker='o',c='r',markersize=20,fillstyle='none',linestyle='None',label='Real Data')
			elif choices['CensoredData'] and choices['CensoredCut']=='inf':pl.plot(3,4,marker='o',c='r',markersize=20,fillstyle='none',linestyle='None',label='Real Data')
			else: raise Exception()
		except:
			if 'Full' in file or 'CensInf' in file or 'BVcut1.0_Cens' in file or 'BVcut0.3_Cens' in file:
				pl.plot(3,4,marker='o',c='r',markersize=20,fillstyle='none',linestyle='None',label='Real Data')
			elif 'Cens1.0' in file: pl.plot(3,0,marker='o',c='r',markersize=20,fillstyle='none',linestyle='None',label='Real Data')

		pl.scatter(Ncens['mid'],Ncens['high'],alpha=0.2,label='Simulations')

		X = pl.gca().get_xlim()
		xticks = np.arange(int(X[0]),int(X[1]+1),int((X[1]+1)/min([int(X[1]+1),5])))
		pl.xticks(xticks)
		pl.gca().xaxis.set_minor_locator(MultipleLocator(1))

		Y = pl.gca().get_ylim()
		yticks = np.arange(int(Y[0]),int(Y[1]+1),int((Y[1]+1)/min([int(Y[1]+1),5])))
		pl.yticks(yticks)
		pl.gca().yaxis.set_minor_locator(MultipleLocator(1))

		mmin = max([min(xticks),min(yticks)])
		mmax = min([max(xticks),max(yticks)])
		pl.plot([mmin,mmax],[mmin,mmax],c='black',alpha=0.2)
		if file=='AVGamma_Cens1.0':
			pl.legend(fontsize=FS-4,ncol=1,loc='center right')
		else:
			pl.legend(fontsize=FS-4,ncol=1)

		pl.tick_params(labelsize=FS)
		pl.tight_layout()
		pl.savefig(f"{sbc.plotpath}{file}_NSNBV.pdf",bbox_inches='tight')
		#pl.show()

		print ('Done')
		print ('####'*10)
