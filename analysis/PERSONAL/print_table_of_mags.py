#Append model_files to path directory
import sys
sys.path.append('../model_files/')
#from load_raw_data import *
from birdsnack_model import BIRDSNACK

#Load up light curves of fiducial sample of SNe, and metadata
#dataloader = LOAD_DATA()
#SNSsnpy_fiducial = dataloader.load_SNSsnpy('SNSsnpy_fiducial.pkl')
#dfmeta           = dataloader.dfmeta
#edit_dict        = {}

#Load into Bird-Snack
#bs = BIRDSNACK(loader={'SNSsnpy':SNSsnpy_fiducial}, configname='loader_config.yaml', dfmeta=dfmeta, edit_dict=edit_dict)
bs = BIRDSNACK(configname='loader_config.yaml')
#Get peak magnitudes
bs.get_peak_mags()

dfpeak = bs.DF_M[0]
pblist = bs.choices['preproc_parameters']['pblist']
errstr = bs.choices['preproc_parameters']['errstr']


errcols = [pb+errstr for pb in pblist]

STRS = []
for sn in list(dfpeak.index):
    magrow = dfpeak.loc[sn][pblist].round(3).values
    errrow = dfpeak.loc[sn][errcols].round(3).values


    STRS.append(f"{sn} & " +' '.join([f"{m} & {merr} &" for m,merr in zip(magrow,errrow)])[:-1]+ ' \\\\')

header = 'SN &  '+' '.join([f"${pb}$ & $\\sigma_{pb}$ &" for pb in pblist])[:-1] + ' \\\\'
print (header)
for s in STRS[:5]:
    print (s)
