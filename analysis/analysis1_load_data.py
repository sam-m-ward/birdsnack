#Append model_files to path directory
import sys
sys.path.append('../model_files/')
from model1_load_data import *
#Load up and save raw data as snpyfile dictionaries, load up metadata
dataloader = LOAD_DATA()

from birdsnack_model import BIRDSNACK

for path_file_survey in dataloader.choices['load_data_parameters']['load_path_file_survey']:
    print ('###'*3)
    #if path_file_survey[2]=='Misc':
    bs = BIRDSNACK(loader={'path_file_survey':path_file_survey}, configname='loader_config.yaml', dfmeta=dataloader.dfmeta)
