"""
LOAD_DATA class

Module to load SN metadata as pandas df, and load up SNooPy .txt files with SNooPy, then save pre-loaded snsnpy files

Contains
--------------------
LOAD_DATA class
    inputs: configname='loader_config.yaml'

    Methods are:
        load_all_SNSsnpy()
        get_SNSsnpy(datafolder,SNSfile,surveyname)
        load_SNSsnpy(filename)
        load_meta_data()
        get_SNSsnpy_combined(SURVEYS,overwrite=False):
--------------------

Written by Sam M. Ward: smw92@cam.ac.uk

"""

from glob import glob
import os,pickle,re,yaml
import pandas as pd
import snpy
from miscellaneous import ensure_folders_to_file_exist

class LOAD_DATA:
    """
	LOAD_DATA Class Object

    Class object that takes snpytxtfiles and meta data from the datapath
    Then, stores snpyfiles, ready for pre-processing cuts


	Parameters
	----------
	configname: str
		name of config.yaml file used for analysis

	Methods
	----------
	load_all_SNSsnpy()
		returns DF_C, DF_M, each of which is a {tc:df} dict of GP apparent colours or magnitudes, respectively

	get_SNSsnpy(datafolder,SNSfile,surveyname)
        loads up snpytxtfiles from a specific survey and saves files (reduces time complexity by performing this operation only once)

    load_meta_data()
    	loads up SN metadata, such as eg host masses, spectroscopic classification, using data compiled from literature
	"""

    def __init__(self,configname='loader_config.yaml',load_all_SNSsnpy=False):
        """
        Initialisation

        Set pathnames
        Load up config.yaml choices
        Load up metadata
        Load up snpytxtfiles and save as dictionary
        Convert snpyfiles to snana lcs

        Parameters
        ----------
        configname : str (optional; default='loader_config.yaml')
            .yaml file that sets analysis choices

        load_all_SNSsnpy : bool (optional; default=True)
            if True, automatically loads up txt files into snpy when class is called (saves on time complexity as this only needs to be done once)
        """
        #Set Configname
        self.configname   = configname

        #Load up config.yaml choices
        with open(self.configname) as f:
            self.choices = yaml.load(f, Loader=yaml.FullLoader)

        #Set Pathname
        self.rootpath     = self.choices['rootpath']
        self.analysispath = self.rootpath+'analysis/'
        self.datapath     = self.rootpath+'data/'
        self.productpath  = self.rootpath+'products/'
        self.SNSpath      = self.productpath+'snpy_SNS/'
        self.snanapath    = self.SNSpath+'snana_copies/'
        for path in [self.analysispath,self.datapath,self.productpath,self.SNSpath,self.snanapath]: ensure_folders_to_file_exist(path)

        #Load up metadata
        self.load_meta_data()

        if load_all_SNSsnpy:
            #Load up snpytxtfiles into snpy (time complex)
            self.load_all_SNSsnpy()

    def load_all_SNSsnpy(self):
        """
        Load All SNS snpy files

        Method to load up snpy lc files from various surveys, and save dictionaries of:
        SNSsnpy = {SN:{'lc':snpyfile,'survey':surveyname}} in self.SNSpath

        End Product(s)
        ----------
        for each survey saves SNSsnpy dictionary
        """
        for x in self.choices['load_data_parameters']['load_path_file_survey']:
            path,file,survey = x[:]
            self.get_SNSsnpy(path,file,survey)

    def get_SNSsnpy(self,datafolder,SNSfile,surveyname,returner=False):
        """
        Get SNSsnpy

        Method to load up a dictionary of snpy files (important because snpy.get_sn is time consuming)
        Can put in new entries, checks to see each file in lcspath is in SNS by matching SNS keys with lcspaths
        saves SNSsnpy = {SN:{'lc':snpyfile,'survey':surveyname}} in self.SNSpath

        Parameters
        ----------
        datafolder: str
            folder where snpytxtfiles are

        SNSfile: str
            name of SNS.pkl

        surveyname: str
            name of survey

        returner: bool (optional; default=False)
            if True, return SNSsnpy

        End Product(s)
        ----------
        SNS: dict
            SNSsnpy = {SN:{'lc':snpyfile,'survey':surveyname}} in self.SNSpath
        """

        lcspath       = self.datapath + datafolder
        SNS_path      = self.SNSpath  + SNSfile
        exclude_paths = [f'{self.datapath}{path}' for path in self.choices['load_data_parameters']['exclude_paths']]
        try:
            snspath = glob(lcspath+'*_snpy.txt')
            snspath.sort()

            SNS = None
            with open(SNS_path,'rb') as f:
                SNS = pickle.load(f)

            not_loaded = []
            for snpath in snspath:
                for isn,sn in enumerate(SNS):
                    if sn in snpath:        break
                    elif isn==len(SNS)-1:   not_loaded.append(snpath)
            if not_loaded!=[]:
                not_loaded = [path for path in not_loaded if path not in exclude_paths]
                if not_loaded!=[]:
                    print (f'SNe Not Loaded: {not_loaded}')
                    snspath = not_loaded
                    raise Exception()
        except Exception as e:
            print (e)
            if SNS is None: SNS = {}
            faillist = []
            for snfile in snspath:
                try:
                    snsnpy            = snpy.get_sn(snfile)
                    snsnpy.name       = snsnpy.name.replace('SN','')
                    SNS[snsnpy.name]  = {'lc':snsnpy,'survey':surveyname}
                    with open(SNS_path,'wb') as f:
                        pickle.dump(SNS,f)
                except Exception as e:
                    print (e)
                    print ('Failed:',snfile)
                    faillist.append(snfile)
            print ('Final SNS:',SNS)
            print (f'faillist (Len: {len(faillist)}):',faillist)
            with open(SNS_path,'wb') as f:
                pickle.dump(SNS,f)
        self.SNSsnpy = SNS
        if returner: return self.SNSsnpy

    def load_SNSsnpy(self,filename):
        """
        Load SNSsnpy using filename

        Parameters
        ----------
        filename : str
            e.g. 'SNSsnpy_fiducial.pkl'

        Returns
        ----------
        SNSsnpy : dict
            {sn:{'lc':lc,'survey':survey}}
        """
        with open(self.SNSpath+filename,'rb') as f:
            SNSsnpy = pickle.load(f)
        return SNSsnpy


    def load_meta_data(self):
        """
        Load Meta Data

        Method to load up meta data of Spectroscopic Sub-class, Host galaxy mass, from various sources

        End Product(s)
        ----------
        self.df_combined: pandas df
            contains SN name, host galaxy stellar masses, spectroscopic sub-classification, and references
        """

        #CSPDR3 Dataframes
        df_spectral_types_cspdr3		   = pd.read_csv(f'{self.datapath}cspdr3/Table2Krisciunas17/table2_edited.txt',sep='\t')
        df_host_masses_cspdr3    		   = pd.read_csv(f'{self.datapath}cspdr3/TableC1Uddin20/TableC1Uddin20.txt', sep='\s*&\s*', names=['SN','Mlow','Mbest','Mhigh','ProjDist'])
        df_host_masses_cspdr3['ProjDist']  = df_host_masses_cspdr3['ProjDist'].apply(lambda x: x.replace('\\',''))
        df_host_masses_cspdr3['HostMass_Source'] = 'cspdr3' ; df_spectral_types_cspdr3['SpecType_Source'] = 'cspdr3'

        #CfA Dataframes
        df_host_masses_CfA			   = pd.read_csv(f'{self.datapath}CfA/Kelly10HostMasses.csv',sep='\s+&\s+')
        df_spectral_types_CfA    	   = pd.read_csv(f'{self.datapath}CfA/Friedmann15SpecClasses.csv', sep=r'\s*\\ & \\\s*', names=['SN','RA','Dec','HostName','HostMorphology','zhelio','sigmazhelio','zref','Discoveryref','Discoverer','Typeref','Spectral_Type','\\'])
        def specmapper(x): return 'normal' if x=='Ia' else 'peculiar'
        df_spectral_types_CfA['SubtypeIa']     = df_spectral_types_CfA['Spectral_Type'].apply(specmapper)
        df_spectral_types_CfA['SN']            = df_spectral_types_CfA['SN'].apply(lambda x: x.replace('\sn{}',''))
        df_host_masses_CfA['HostMass_Source'] = 'CfA' ; df_spectral_types_CfA['SpecType_Source'] = 'CfA'

        #Ponder21 Masses
        #No mass entries in the Ponder21 table for ['00ca','01bt','01cz','01el','11by'] unfortunately
        df_host_masses_Ponder21 = pd.read_csv(f'{self.datapath}CfA/Ponder21masses.csv',sep='\s*&\s*')
        df_host_masses_Ponder21['SN'] = df_host_masses_Ponder21['SN'].apply(lambda x: x.replace('SN~',''))
        df_host_masses_Ponder21['HostMass_Source'] = 'P21'

        #Uddin17 Masses
        df_host_masses_Uddin17  = pd.read_csv(f'{self.datapath}CfA/Uddin17host.csv', sep='\s+')
        df_host_masses_Uddin17['SN'] = df_host_masses_Uddin17['SN'].apply(lambda x: x.replace('sn',''))
        df_host_masses_Uddin17['HostMass_Source'] = 'U17'
        #Uddin20 exclude: 04gc,07A,07if,07mm,08bf,08ff because of poor galaxy association
                        #06dd excluded because lack of peak-time photometry
                        #07sr excluded because of galaxy association

        #Johansson21 Masses
        df_host_masses_Johansson21 = pd.read_csv(f'{self.datapath}RATIR/Johansson21/HostMasses.csv',sep='\s+&\s+')
        df_host_masses_Johansson21['Mbest'] = df_host_masses_Johansson21['Mbest'].apply(lambda x: x.replace('\t',''))
        df_host_masses_Johansson21['Mbest'] = df_host_masses_Johansson21['Mbest'].apply(lambda x: x.replace('\\',''))
        df_host_masses_Johansson21['SN']    = df_host_masses_Johansson21['SN'].apply(lambda x: re.findall(r'(\w+)\s+\(SN\\,\w+\)',x)[0] if x[-1]==')' else x)
        df_host_masses_Johansson21['SN']    = df_host_masses_Johansson21['SN'].apply(lambda x: x.replace('iPTF',''))
        df_host_masses_Johansson21 = df_host_masses_Johansson21.append(pd.DataFrame(data = {'SN':['2014J'],'Mbest':[9]}))#Assume SN2014J is low mass; Greco 2012"Measurement of the mass and stellar population distribution in M82 with LBT"; dynamical mass approaches 10^10 (but this inevitably includes e.g. DM as well
        df_host_masses_Johansson21['HostMass_Source'] = 'J21'

        #RATIR Spec Types
        df_spectral_types_RATIR = df_host_masses_Johansson21[['SN']].copy()
        def get_ratir_spec(sn):
            #if sn in ['14atg','14apg']:#These ones excluded in J21
            #	df_spectral_types_RATIR[df_spectral_types_RATIR['SN']==sn]['SubtypeIa'] = 'peculiar'
            if sn in ['13abc','13ebh','14atg','14bdn','14apg']:#'13asv','13dge','16abc','17lf'
                return 'peculiar'
            elif sn in ['14ale']:
                return 'unknown'
            else:
                return 'normal'
        df_spectral_types_RATIR['SubtypeIa'] = df_spectral_types_RATIR['SN'].apply(get_ratir_spec)
        df_spectral_types_RATIR['SpecType_Source'] = 'RATIR'

        #Avelino Spectroscopic Types; Assign normal to M20 sample of SNe
        #Misc Spectroscopic Types; Assign normal to all these SNe
        AvelinoSNs                = [re.findall(rf'{self.datapath}BayeSNLCs/snpytxtfiles/(.*)_snpy.txt',f)[0] for f in glob(f'{self.datapath}BayeSNLCs/snpytxtfiles/*_snpy.txt')]
        MiscSNs                   = [re.findall(rf'{self.datapath}Misc/snpytxtfiles/(.*)_snpy.txt',f)[0] for f in glob(f'{self.datapath}Misc/snpytxtfiles/*_snpy.txt')]
        AvelinoM20MiscSNs = AvelinoSNs+MiscSNs
        df_spectral_types_AvelinoM20Misc = pd.DataFrame(data = dict(zip(['SN','SubtypeIa'],[AvelinoM20MiscSNs,['normal' for _ in range(len(AvelinoM20MiscSNs))]] ) ) )
        df_spectral_types_AvelinoM20Misc['SpecType_Source'] = 'AvelinoM20Misc'

        #Neill09 Host Masses
        df_SNhostmapper        = pd.read_csv(f'{self.datapath}Misc/Neill09/Table1.txt',sep='\s+&\s+',names=['SN','host','c','d','e','f'])[['SN','host']]
        df_Neill09_masses      = pd.read_csv(f'{self.datapath}Misc/Neill09/Table2.txt',sep='\s+&\s+',names=['host','typeNum','Agelow','Agebest','Agehigh','Mlow','Mbest','Mhigh','SFRlow','SFRbest','SFRhigh','EBHh'])[['host','Mbest']]
        df_host_masses_Neill09 = df_SNhostmapper.merge(df_Neill09_masses)
        df_host_masses_Neill09['SpecType_Source'] = 'Neill09'

        #Rose19Table11 "Think Global, Act Local: The Influence of Environment Age and Host Mass on Type Ia Supernova Light Curves"
        df_host_masses_Rose19 = pd.DataFrame(data={'SN':['2011by','2011fe'],'Mbest':[9.8,9.9]})
        df_host_masses_Rose19['SpecType_Source'] = 'Rose19'

        #Merge style may have double entries which conflict (i.e. same SN appears in multiple tables), where entries are the same only pick out once, else merge list of entries into one string
        df_host_masses    = pd.concat([df_host_masses_cspdr3,df_host_masses_CfA,df_host_masses_Ponder21, df_host_masses_Uddin17, df_host_masses_Johansson21, df_host_masses_Neill09],axis=0)
        df_spectral_types = pd.concat([df_spectral_types_cspdr3,df_spectral_types_CfA,df_spectral_types_RATIR, df_spectral_types_AvelinoM20Misc],axis=0)

        df_host_masses    = df_host_masses[['SN','Mbest','HostMass_Source']]
        df_spectral_types = df_spectral_types[['SN','SubtypeIa','SpecType_Source']]
        df_combined = df_host_masses.merge(df_spectral_types, on='SN',how='outer')

        self.dfmeta = df_combined

    def get_SNSsnpy_combined(self,SURVEYS,overwrite=False):
        """
        Get SNSsnpy_combined

        Method for getting SNSsnpy by combining multiple surveys according to pecking order

        Parameters
        ----------
        SURVEYS : dict
            key,values are {'retained_lcs':lcs,'trimmed_lcs':bs.trimmed_lcs, 'reasons':reasons}

        overwrite : bool (optional; default=False)
            if True, rewrite SNSsnpy_combined file even if if already exists

        End Product(s)
        ----------
        save files for:
            SURVEYS : dict
            retained_sns : dict with {sn:survey}
            SNSsnpy_combined : dict with {sn:{lc:snpyfile,survey:survey}}
        """
        #Save SURVEYS dictionarys
        if overwrite or not os.path.exists(f'{self.SNSpath}SURVEYS.pkl'):
            with open(f'{self.SNSpath}SURVEYS.pkl','wb') as f:
                pickle.dump(SURVEYS,f)

        #Create/load/save retained_sns dict
        if overwrite or not os.path.exists(f'{self.SNSpath}retained_sns.pkl'):
            #Determine SNe that are retained, keep those highest in pecking order
            retained_sns = {}
            for survey in self.choices['load_data_parameters']['Pecking_Order']:
                for sn in SURVEYS[survey]['retained_lcs']:
                    if sn not in list(retained_sns.keys()):
                        retained_sns[sn] = survey

            #Save retained_sns dicionary as .pkl
            with open(f'{self.SNSpath}retained_sns.pkl','wb') as f:
                pickle.dump(retained_sns,f)
        else:
            with open(f'{self.SNSpath}retained_sns.pkl','rb') as f:
                retained_sns = pickle.load(f)

        #Create/load/save SNSsnpy_combined
        if overwrite or not os.path.exists(f'{self.SNSpath}SNSsnpy_combined.pkl'):
            SNSsnpy_combined = {}
            for survey in self.choices['load_data_parameters']['Pecking_Order']:
                for x in self.choices['load_data_parameters']['load_path_file_survey']:
                    if x[-1]==survey:
                        path,file,survey = x[:]
                        SNSsnpy = self.get_SNSsnpy(path,file,survey,returner=True)

                for sn in retained_sns:
                    if sn in list(SNSsnpy.keys()) and retained_sns[sn]==SNSsnpy[sn]['survey']:
                        SNSsnpy_combined[sn] = {'lc':SNSsnpy[sn]['lc'],'survey':SNSsnpy[sn]['survey']}
            with open(f'{self.SNSpath}SNSsnpy_combined.pkl','wb') as f:
                pickle.dump(SNSsnpy_combined,f)
        else:
            with open(f'{self.SNSpath}SNSsnpy_combined.pkl','rb') as f:
                SNSsnpy_combined = pickle.load(f)

        return SNSsnpy_combined
