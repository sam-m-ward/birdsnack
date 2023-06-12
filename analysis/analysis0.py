import sys
sys.path.append('../model_files/')

#For analysis1
from model1_load_data import *

#Load up metadata including SN spectroscopic subclassification and host galaxy stellar masses
#Also load up all lcs into snpy, and save survey name
dataloader = LOAD_DATA()
