"""
LOAD_DATA class and methods

Module to load SN metadata, and SNooPy .txt files, then convert them to snana files, save files, and store/return as LCs

Contains
--------------------
LOAD_DATA class
    inputs: rootpath='../',configname='default_config.yaml'

    Methods are:
        load_all_SNSsnpy()
        get_SNSsnpy(datafolder,SNSfile,surveyname)
        convert_all_snpy_to_snana()
        convert_snpyfiles_to_snana(SNSfile)
        convert_snpyfile_to_snana(sn,snsnpy)
        load_meta_data()

--------------------
Functions use simple operations

Written by Sam M. Ward: smw92@cam.ac.uk

"""

from glob import glob
import snpy
import pickle,re,yaml,os
import pandas as pd

class LOAD_DATA:
    """
	LOAD_DATA Class Object

    Class object that takes snpytxtfiles and meta data from the datapath
    Then, stores snpyfiles and .snana files, ready for pre-processing cuts


	Parameters
	----------
	rootpath : str
        path/to/root/directory from which e.g. datapath, productspath are accessible

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

    def __init__(self,rootpath='../',configname='default_config.yaml'):
        """
        Initialisation

        Set pathnames
        Load up config.yaml choices
        Load up metadata
        Load up snpytxtfiles and save as dictionary
        Convert snpyfiles to snana lcs

        """
        #Set Pathname
        self.rootpath     = rootpath
        self.configname   = configname
        self.datapath     = self.rootpath+'data/'
        self.productpath  = self.rootpath+'products/'
        self.analysispath = self.rootpath+'analysis/'
        self.SNSpath      = self.productpath+'snpy_SNS/'

        #Load up config.yaml choices
        with open(self.configname) as f:
            self.choices = yaml.load(f, Loader=yaml.FullLoader)

        #Load up metadata
        self.load_meta_data()

        #Load up snpytxtfiles
        self.load_all_SNSsnpy()

        #Convert snpytxtfiles to snana files
        #self.convert_all_snpy_to_snana()

    def load_all_SNSsnpy(self):
        """
        Load All SNS snpy files

        Method to load up snpy lc files from various surveys, and save dictionaries of:
        SNSsnpy = {SN:{'lc':snpyfile,'survey':surveyname}} in self.SNSpath
        """

        self.get_SNSsnpy('cspdr3/snpytxtfiles/','snpy_SNS_CSP.pkl','CSP')
        self.get_SNSsnpy('CfA/snpytxtfiles/','snpy_SNS_CfA.pkl','CfA')
        self.get_SNSsnpy('RATIR/snpytxtfiles/','snpy_SNS_RATIR.pkl','RATIR')
        self.get_SNSsnpy('BayeSNLCs/snpytxtfiles/','snpy_SNS_Avelino.pkl','AvelinoM20')
        self.get_SNSsnpy('Misc/snpytxtfiles/','snpy_SNS_Misc.pkl','Misc')

    def convert_all_snpy_to_snana(self):
        """
        Load All SNSsnpy dictionaries and covert to snana files

        Method to load up SNSsnpy = {SN:{'lc':snpyfile,'survey':surveyname}} in self.SNSpath
        then convert each snpyfile to an snana file
        """

        self.convert_snpyfiles_to_snana('snpy_SNS_CSP.pkl')
        self.convert_snpyfiles_to_snana('snpy_SNS_CfA.pkl')
        self.convert_snpyfiles_to_snana('snpy_SNS_RATIR.pkl')
        self.convert_snpyfiles_to_snana('snpy_SNS_Avelino.pkl')
        self.convert_snpyfiles_to_snana('snpy_SNS_Misc.pkl')

    def get_SNSsnpy(self,datafolder,SNSfile,surveyname):
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

    def convert_snpyfiles_to_snana(self,SNSfile):
        """
        Get Light Curves from SNooPy

        Simple function to iterate over the collection of snoopy objects and collect the snana lcs

        Parameters
        ----------
        SNSsnpy: dict
            {snname:snpy.get_sn(snfile)} for each sn in lcspath

        outputdir: str
            path/to/dir/where_snana_lcs_are_saved

        Returns
        ----------
        lcs : dict
            {sn:lc} where sn is str of snname, and lc is :py:class:`astropy.table.Table` light curve object
        """
        SNS_path = self.SNSpath + SNSfile
        with open(SNS_path,'rb') as f:
            SNSsnpy = pickle.load(f)

        lcs = {}
        for sn in SNSsnpy:
            snsnpy  = SNSsnpy[sn]['lc']
            lcs[sn] = {'lc':self.convert_snpyfile_to_snana(sn,snsnpy),'survey':SNSsnpy[sn]['survey']}

        return lcs

    def convert_snpyfile_to_snana(self,sn,snsnpy):
        """
        Get Light Curve from snpy object

        Function loads up an snana lc given a sn name, and creates this file from a snpy object if it doesn't already exist

        Parameters
        ----------
        sn: str
            name of the sn

        snsnpy: snpy.sn.sn object
            snpyfile object

        Returns
        ----------
        lc: :py:class:`astropy.table.Table`
            light curve object
        """
        import os
        import snanaio as io

        path = self.SNSpath+'snana_copies/'+sn+'.snana.dat'
        if not os.path.exists(path):
            mjd,flt,mag,magerr = [],[],[],[]
            for filt in snsnpy.allbands():
                mjd.extend(list(snsnpy.data[filt].MJD))
                flt.extend([filt for _ in range(len(snsnpy.data[filt].MJD))])
                mag.extend(list(snsnpy.data[filt].mag))
                magerr.extend(list(snsnpy.data[filt].e_mag))
            snname  = sn; tmax=0; z_helio = snsnpy.z; z_cmb=snsnpy.get_zcmb(); z_cmb_err=0; ebv_mw = snsnpy.EBVgal
            io.write_snana_lcfile(self.SNSpath+'snana_copies/', snname, mjd, flt, mag, magerr, tmax, z_helio, z_cmb, z_cmb_err, ebv_mw,ra=snsnpy.ra,dec=snsnpy.decl)
        sn , lc = io.read_snana_lcfile(path)
        os.remove(path)#Delete the newly created file straight away
        return lc

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

        #Avelino Spectroscopic Types; #Assign normal to M20 sample of SNe
        AvelinoSNs                = [re.findall(rf'{self.datapath}BayeSNLCs/snpytxtfiles/(.*)_snpy.txt',f)[0] for f in glob(f'{self.datapath}BayeSNLCs/snpytxtfiles/*_snpy.txt')]
        df_spectral_types_Avelino = pd.DataFrame(data = dict(zip(['SN','SubtypeIa'],[AvelinoSNs,['normal' for _ in range(len(AvelinoSNs))]] ) ) )

        #Misc Spectroscopic Types; Assign normal to all these SNe
        MiscSNs                   = [re.findall(rf'{self.datapath}Misc/snpytxtfiles/(.*)_snpy.txt',f)[0] for f in glob(f'{self.datapath}Misc/snpytxtfiles/*_snpy.txt')]
        df_spectral_types_Avelino = pd.DataFrame(data = dict(zip(['SN','SubtypeIa'],[MiscSNs,['normal' for _ in range(len(MiscSNs))]] ) ) )
        df_spectral_types_Avelino['SpecType_Source'] = 'AvelinoM20Misc'

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
        df_spectral_types = pd.concat([df_spectral_types_cspdr3,df_spectral_types_CfA,df_spectral_types_RATIR, df_spectral_types_Avelino],axis=0)

        df_host_masses    = df_host_masses[['SN','Mbest','HostMass_Source']]
        df_spectral_types = df_spectral_types[['SN','SubtypeIa','SpecType_Source']]
        df_combined = df_host_masses.merge(df_spectral_types, on='SN')

        self.dfmeta = df_combined
