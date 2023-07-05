def exponential(lam,x):
    return np.array([0 if xi<=0 else xi for xi in lam*np.exp(-lam*x)])

def gamma(alpha,beta,x):
    num  = (x**(alpha-1))*np.exp(-beta*x)*(beta**alpha)
    den  = math.gamma(alpha)
    func = num/den
    return func

#Append model_files to path directory
import math,sys
import numpy as np
sys.path.append('../../model_files/')
from load_raw_data import *
from birdsnack_model import BIRDSNACK
from HBM_preprocessing import HBM_preprocessor

#Load up light curves of fiducial sample of SNe, and metadata
dataloader = LOAD_DATA(rootpath='../')
SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')
dfmeta           = dataloader.dfmeta

#Get the indices of SNe in low-reddening sample
edit_dict = {'rootpath':'../../','analysis_parameters':{'CensoredData':True}}
bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='../loader_config.yaml', dfmeta=dfmeta, edit_dict=edit_dict)
bs.get_peak_mags()
modelloader = HBM_preprocessor(bs.choices, bs.DF_M)
modelloader.get_censored_data()
lowinds = [i for i,sn in enumerate(list(bs.DF_M[0].index)) if sn in modelloader.RetainedSNe]

import matplotlib.pyplot as pl
#Files for inspection
FITS = ['AVExp_Full',
        'AVExp_LowBV',
        'AVGamma_Full',
        'AVGamma_LowBV']

FITS = FITS + [ff.replace('Full','CensInf').replace('LowBV','Cens1.0') for ff in FITS ]
modes = ['Exponential','Gamma','Exp. w/ Cens. SNe', 'Gamma w/ Cens. SNe']
Modemap = dict(zip(FITS,['Science ' for _ in range(4)] + ['CensoredData' for _ in range(4)]))
AVpriors = dict(zip(FITS,['Exp','Exp','Gamma','Gamma','Exp','Exp','Gamma','Gamma']))


FS=18 ; MS = 8 ; STR = '' ; CS = 4
#'''
pl.figure(figsize=(8,6))
pl.title('Mean $A_V^s$ in Low-Reddening Sample',fontsize=FS)
counter = -1
for file in FITS:
    Mode = Modemap[file]
    counter += 1
    if counter%CS==0: str = f'{modes[counter//CS]} & '
    with open(f"{bs.FITSpath}FIT{file}.pkl",'rb') as f:
        FIT = pickle.load(f)
    df     = FIT['df']
    AVcols = [col for col in df.columns if 'AVs.' in col]
    if 'Full' in file: AVcols = [col for col in AVcols if float(col.split('.')[1])-1 in lowinds]
    if Mode == 'CensoredData':
        SC = FIT['stan_data']['SC'] ; S = FIT['stan_data']['S']
        AVcols = [col for col in AVcols if float(col.split('.')[1])-1 < (S-SC)]
    AVs = df[AVcols].mean(axis=1)
    str += f'${round(AVs.median(),2)}\\pm{round(AVs.std(),2)}$ & '
    if (counter+1)%CS==0:
        str = str[:-2] +  ' \\\\ '
        STR += str
    xerr=AVs.std(); yerr=df['mu_RV'].std()
    xerr = [[AVs.quantile(0.5)-AVs.quantile(0.16)],[AVs.quantile(0.84)-AVs.quantile(0.5)]]
    yerr = [[df['mu_RV'].quantile(0.5)-df['mu_RV'].quantile(0.16)],[df['mu_RV'].quantile(0.84)-df['mu_RV'].quantile(0.5)]]
    pl.errorbar(AVs.median(),df['mu_RV'].median(),xerr=xerr,yerr=yerr,
    marker=['o','x','s','d'][2*(counter//4)+counter%2],
    c=f'C{int(counter//2)}',
    label=f'{modes[counter//2]}'*(counter%2==0),
    markersize=MS,alpha=0.8,capsize=2)
###########
X = pl.gca().get_xlim()
Y = pl.gca().get_ylim()
for cntr in np.arange(CS):
    pl.errorbar(-1,-1,marker=['o','x','s','d'][cntr],c='black',markersize=MS,alpha=0.7,capsize=2,label=['69 [Full Samp.]','62 [Low-Red. Samp.]','62(+7)','62(+3)'][cntr])
pl.xlim(X)
pl.ylim(Y)
pl.ylim([1.5,5.0])
pl.xlabel(r'$\overline{A_V^s}$ (mag)',fontsize=FS)
pl.ylabel(r'$\mu_{R_V}$',fontsize=FS)
pl.legend(fontsize=FS-6,title=r'$A^s_V$ Prior Distribution          N$_{\rm{SNe}}$ (N$_{\rm{Censored}}$)       ',
title_fontsize=FS-4,ncol=2,loc='upper left',framealpha=0.5)
pl.tick_params(labelsize=FS)
pl.tight_layout()
pl.savefig(f"{bs.plotpath}AVmeansNoCensAndCens.pdf",bbox_inches='tight')
pl.show()
print (STR)
#'''

#'''
fig,axs = pl.subplots(2,1,figsize=(8,4*2),sharex=True)
ax = fig.axes
ax[0].set_title('$A_V^s$ Population Distributions',fontsize=FS)
modes = ['Exponential','Gamma','Exponential', 'Gamma']
AVmax = 2
x     = np.linspace(0,AVmax,1000)
counter=-1
for file in FITS:
    Mode = Modemap[file]
    counter += 1
    if counter%CS==0: str = f'{modes[counter//CS]} & '
    with open(f"{bs.FITSpath}FIT{file}.pkl",'rb') as f:
        FIT = pickle.load(f)
    df     = FIT['df']
    AVsprior_params = [par for par in ['tauA','nu'] if 'RV' not in par]
    y     = np.zeros(len(x))
    tauAs = df['tauA'].values
    try: nus = df['nu'].values
    except: nus= np.zeros(len(tauAs))
    #########
    Nthinsize = 4000
    samples = [tauAs,nus]
    for p in range(len(samples)):
        di         = int((len(samples[p]))/Nthinsize)
        samples[p] = np.array([samples[p][_] for _ in np.arange(0,len(samples[p]),di) if di<=len(samples[p]-1)])
        samples[p] = [np.median(samples[p])]
    tauAs,nus = samples[:]
    #########
    if AVpriors[file]=='Exp':
        Title = r'$A_V^s \sim$Exp$(\tau_A)$'
        for tauA in tauAs:
            y  += exponential(1/tauA,x)
    elif AVpriors[file]=='Gamma':
        Title = r'$A_V^s \sim$Gamma$(\nu_A,\tau_A)$'
        for tauA,nu in zip(tauAs,nus):
            y  += gamma(nu,1/tauA,x)

    ir = counter//4
    ax[ir].plot(x,y,c=f'C{counter//2}',linestyle={0:'-',1:'--'}[counter%2],label=modes[counter//2]*(counter%2==0))
    ax[ir].set_yticks([])
    ax[ir].set_xlim([0,AVmax])
    ax[ir].set_ylim([0,3])
    ax[ir].set_yticklabels([])
    ax[ir].tick_params(labelsize=FS)
    if (counter//2)%2==1:
        ax[ir].plot([0,0.1],[-0.1,-0.1],color='black',linestyle={0:'-',1:'--'}[counter%2],label=['69 [Full Samp.]','62 [Low-Red. Samp.]','62(+7)','62(+3)'][(counter//2)-1+counter%2])

ax[0].legend(fontsize=FS-2,ncol=2,title=r'$A_V^s$-Prior'+'\t\t\t'+r'$N_{\rm{SNe}}$[Sample]', title_fontsize=FS)
ax[1].legend(fontsize=FS-2,ncol=2,title=r'       $A_V^s$-Prior'+'\t      '+r'$N_{\rm{SNe}}(N_{\rm{Cens}})$', title_fontsize=FS)
ax[0].annotate(r'No Censored SNe',        fontsize=FS+1,xy=(0.5,0.375),xycoords='axes fraction',weight='bold')
ax[1].annotate(r'With Censored SNe',        fontsize=FS+1,xy=(0.5,0.375),xycoords='axes fraction',weight='bold')
ax[-1].set_xlabel(r'$A_V^s$',fontsize=FS)

fig.add_subplot(111, frameon=False)#For common Y label
pl.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
pl.ylabel(r'Posterior-Median $A_V^s$ Distributions',fontsize=FS)
pl.tight_layout()
pl.savefig(f"{bs.plotpath}AVsdistposterioraveragedNoCensCens.pdf",bbox_inches='tight')
pl.show()
#'''
