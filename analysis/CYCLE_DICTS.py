CYCLE_DICT_Science = {
                'COMMON_CHANGES' : {'newdict':{},'HBMappender':''},

                'RUNS' : {

				'AVexp_Full'                 :       {    'newdict': {'AVprior':'Exp'},
													 'label':'$A^s_V \sim \\rm{Exp}(\\tau_A)$'},
				'AVGamma_Full'               :       {    'newdict': {'AVprior':'Gamma','n_sampling':5000},
													 'label':'$A_V^s \sim $ Gamma$(\\nu_A,\\tau_A)$\\tnote{b}'},
                'AVexp_lowBV'           :       {    'newdict': {'AVprior':'Exp','BVcut':True},
													 'label':'$A^s_V \sim \\rm{Exp}(\\tau_A)$'},
				'AVGamma_lowBV'         :       {    'newdict': {'AVprior':'Gamma','n_sampling':5000,'BVcut':True},
													 'label':'$A_V^s \sim $ Gamma$(\\nu_A,\\tau_A)$\\tnote{b}'},
				}
            }

CYCLE_DICT_CensoredData = {
                'COMMON_CHANGES' : {'newdict':{'CensoredData':True,'CensoredCut':1.0},'HBMappender':'Cens1.0'},

                'RUNS' : {

				'AVexp'                 :       {    'newdict': {'AVprior':'Exp'},
													 'label':'$A^s_V \sim \\rm{Exp}(\\tau_A)$'},
				'AVGamma'               :       {    'newdict': {'AVprior':'Gamma','n_warmup':2000,'n_sampling':10000},
													 'label':'$A_V^s \sim $ Gamma$(\\nu_A,\\tau_A)$\\tnote{e}'},
                ###############
				'Central_Lam'				:       {   'newdict': {'lam_choice':'central','n_sampling':3000},
											'label':'Central-$\\lambda$'},
				'PreProc_Interpflts_BVriJH'			:       {   'newdict': {'interpflts':'BVriJH','n_sampling':3000,'DF_savekey':'interpfltsBVriJH'},
											'label':'$BVriJH$'},
				'PreProc_EBVMW_PreSubtract'			:       {   'newdict': {'apply_EBVMW_corrections':'Presubtract','n_sampling':3000,'DF_savekey':'EBVMWpresubtract'},
											'label':'Pre-subtract $E(B-V)_{\\rm{MW}}$'},
				'PreProc_1DGPTMax_KCorrNoMangle'   :       {   'newdict': {'mangle': False,'n_sampling':3000,'DF_savekey':'KcorrNoMangle'},
											'label':'No Mangling'},
				'PreProc_2DGPTmax'                 :       {   'newdict': {'Tmax_method':'2DGP','n_sampling':3000,'DF_savekey':'2DGPTmax'},
											'label':'2DGP'},
				'PreProc_SNPYTmax'                 :       {   'newdict': {'Tmaxchoicestr': 'Tmax_snpy_fitted','n_sampling':3000,'DF_savekey':'snpyTmax'},
											'label':'\\textsc{SNooPy}'},

				'PreProc_2Dfl'          :       {   'newdict': {'method' : '2DGP', 'bright_mode':'flux','n_sampling':3000,'DF_savekey':'2DGPfluxinterp'},
											'label':'2DGP Flux Interp.'},
				'PreProc_1Dmg'          :       {   'newdict': {'method' : '1DGP', 'bright_mode':'mag','n_sampling':3000,'DF_savekey':'1DGPmaginterp'},
											'label':'1DGP Mag Interp.'},
				'PreProc_1Dfl'          :       {   'newdict': {'method' : '1DGP', 'bright_mode':'flux','n_sampling':3000,'DF_savekey':'1DGPfluxinterp'},
											'label':'1DGP Flux Interp.'},
                ###############
				'HighMass'   			:		{   'newdict': {'mass_mode':'high_masses','n_sampling':3000},
											'label':'High-Mass'},
				'LowMass'    			:		{   'newdict': {'mass_mode':'low_masses','n_sampling':3000},
											'label':'Low-Mass'},
                ###############
				'AdjCols'      :       {    'newdict': {'DataTransformation':'Adjacent','IntrinsicModel':'Adjacent',
															'n_sampling':2000,'n_warmup':2000},
											'label':'Adjacent Colours'},
				'BXCols'          :       {    'newdict': {'DataTransformation':'B-X','IntrinsicModel':'B-X',
															'n_sampling':2000,'n_warmup':2000},
											'label':'$B-X$'},
				'XHCols'          :       {    'newdict': {'DataTransformation':'X-H','IntrinsicModel':'X-H',
													'n_sampling':2000,'n_warmup':2000},
											'label':'$X-H$'},
                ###############
                'LCShapeInc'         	:       {   'newdict': {'include_LCshape': True,'trim_on_extras':True,'DF_savekey':'IncLCShape',
                                                            'n_sampling':3000},
													'label':'With LC Shape'},

				'AdjCols_wLCshape'      :       {    'newdict': {'include_LCshape': True,'DataTransformation':'Adjacent','IntrinsicModel':'Adjacent','trim_on_extras':True,'DF_savekey':'IncLCShape',
															'n_sampling':2000,'n_warmup':2000},
											'label':'Adjacent Colours'},
				'BXCols_wLCshape'          :       {    'newdict': {'include_LCshape': True,'DataTransformation':'B-X','IntrinsicModel':'B-X','trim_on_extras':True,'DF_savekey':'IncLCShape',
															'n_sampling':2000,'n_warmup':2000},
											'label':'$B-X$'},
				'XHCols_wLCshape'          :       {    'newdict': {'include_LCshape': True,'DataTransformation':'X-H','IntrinsicModel':'X-H','trim_on_extras':True,'DF_savekey':'IncLCShape',
													'n_sampling':2000,'n_warmup':2000},
											'label':'$X-H$'},
                ###############
				'Deviations_NoIntVar'   :		 {    'newdict': {
															'n_sampling':1000,'n_warmup':1000,'include_residuals':False},
											'label':'Deviations'},
				'AdjCols_NoIntVar'      :       {    'newdict': {'DataTransformation':'Adjacent','IntrinsicModel':'Adjacent',
															'n_sampling':1000,'n_warmup':1000, 'include_residuals':False},
											'label':'Adjacent Colours'},
				'BXCols_NoIntVar'          :       {    'newdict': {'DataTransformation':'B-X','IntrinsicModel':'B-X',
															'n_sampling':1000,'n_warmup':1000, 'include_residuals':False},
											'label':'$B-X$'},
				'XHCols_NoIntVar'          :       {    'newdict': {'DataTransformation':'X-H','IntrinsicModel':'X-H',
															'n_sampling':1000,'n_warmup':1000, 'include_residuals':False},
											'label':'$X-H$'},

                }
            }
