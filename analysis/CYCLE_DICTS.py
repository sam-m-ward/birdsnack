CYCLE_DICT_Science = {
				'COMMON_CHANGES' : {'newdict':{},'HBMappender':''},

				'RUNS' : {

				'AVExp_Full'                 :       {    'newdict': {'AVprior':'Exp'},
													 'label':'$A^s_V \sim \\rm{Exp}(\\tau_A)$ Full'},
				'AVGamma_Full'               :       {    'newdict': {'AVprior':'Gamma','n_sampling':5000},
													 'label':'$A_V^s \sim $ Gamma$(\\nu_A,\\tau_A)$ Full'},
				'AVExp_lowBV'           :       {    'newdict': {'AVprior':'Exp','BVcut':True},
													 'label':'$A^s_V \sim \\rm{Exp}(\\tau_A)$ low BV'},
				'AVGamma_lowBV'         :       {    'newdict': {'AVprior':'Gamma','n_sampling':5000,'BVcut':True},
													 'label':'$A_V^s \sim $ Gamma$(\\nu_A,\\tau_A)$ low BV'},
				#'uBVriJH'				:		{	 'newdict':{'DF_savekey':'uBVriJH','pblist':[s for s in 'uBVriJH'],'lam_choice':'central'},
				#										'label':'$uBVriJH$'},
				}
			}

CYCLE_DICT_BVcut1 = {
				'COMMON_CHANGES' : {'newdict':{'CensoredData':False,'BVcut':True,'BVcutval':1.0},'HBMappender':'BVcut1.0'},

				'RUNS' : {

				'AVExp'                 :       {    'newdict': {'AVprior':'Exp','n_sampling':2000},
													 'label':'$A^s_V \sim \\rm{Exp}(\\tau_A)$'},
				'AVGamma'               :       {    'newdict': {'AVprior':'Gamma','n_sampling':10000},
													 'label':'$A_V^s \sim $ Gamma$(\\nu_A,\\tau_A)$'},
				###############
				'uBVriJH'				:		{	 'newdict':{'DF_savekey':'uBVriJH','pblist':[s for s in 'uBVriJH'],'lam_choice':'central'},
														'label':'$uBVriJH$'},
				###############
				'Central_Lam'				:       {   'newdict': {'lam_choice':'central','n_sampling':3000},
											'label':'Central-$\\lambda$'},
				'PreProc_Interpflts_BVriJH'			:       {   'newdict': {'interpflts':'BVriJH','n_sampling':3000,'DF_savekey':'interpfltsBVriJH'},
											'label':'Interp. Filts. $BVriJH$'},
				'PreProc_EBVMW_PreSubtract'			:       {   'newdict': {'apply_EBVMW_corrections':'Presubtract','n_sampling':3000,'DF_savekey':'EBVMWpresubtract'},
											'label':'Pre-subtract $E(B-V)_{\\rm{MW}}$'},
				'PreProc_KCorrNoMangle'   :       {   	'newdict': {'mangle': False,'n_sampling':3000,'DF_savekey':'KcorrNoMangle'},
											'label':'K-corrections No Mangling'},
				'PreProc_2DGPTmax'                 :       {   'newdict': {'Tmax_method':'2DGP','n_sampling':3000,'DF_savekey':'2DGPTmax'},
											'label':'2DGP $T_{B;\,\\rm{max}}$'},
				'PreProc_SNPYTmax'                 :       {   'newdict': {'Tmaxchoicestr': 'Tmax_snpy_fitted','n_sampling':3000,'DF_savekey':'snpyTmax'},
											'label':'\\textsc{SNooPy} $T_{B;\,\\rm{max}}$'},

				'PreProc_2Dfl'          :       {   'newdict': {'method' : '2DGP', 'bright_mode':'flux','n_sampling':3000,'DF_savekey':'2DGPfluxinterp'},
											'label':'2DGP Flux Interp.'},
				'PreProc_1Dmg'          :       {   'newdict': {'method' : '1DGP', 'bright_mode':'mag','n_sampling':3000,'DF_savekey':'1DGPmaginterp'},
											'label':'1DGP Mag Interp.'},
				'PreProc_1Dfl'          :       {   'newdict': {'method' : '1DGP', 'bright_mode':'flux','n_sampling':3000,'DF_savekey':'1DGPfluxinterp'},
											'label':'1DGP Flux Interp.'},
				###############
				'AdjCols'      :       {    'newdict': {'DataTransformation':'Adjacent','IntrinsicModel':'Adjacent',
															'n_sampling':5000,'n_warmup':2000},
											'label':'Adjacent Colours'},
				'BXCols'          :       {    'newdict': {'DataTransformation':'B-X','IntrinsicModel':'B-X',
															'n_sampling':5000,'n_warmup':2000},
											'label':'$B-X$ Colours'},
				'XHCols'          :       {    'newdict': {'DataTransformation':'X-H','IntrinsicModel':'X-H',
													'n_sampling':10000,'n_warmup':2000},
											'label':'$X-H$ Colours'},
				###############
				'LCShapeInc'         	:       {   'newdict': {'include_LCshape': True,'trim_on_extras':True,'DF_savekey':'IncLCShape',
															'n_sampling':3000},
													'label':'Deviations w/ LC Shape'},

				'AdjCols_wLCshape'      :       {    'newdict': {'include_LCshape': True,'DataTransformation':'Adjacent','IntrinsicModel':'Adjacent','trim_on_extras':True,'DF_savekey':'IncLCShape',
															'n_sampling':2000,'n_warmup':2000},
											'label':'Adjacent Colours w/ LC Shape'},
				'BXCols_wLCshape'          :       {    'newdict': {'include_LCshape': True,'DataTransformation':'B-X','IntrinsicModel':'B-X','trim_on_extras':True,'DF_savekey':'IncLCShape',
															'n_sampling':2000,'n_warmup':2000},
											'label':'$B-X$ Colours w/ LC Shape'},
				'XHCols_wLCshape'          :       {    'newdict': {'include_LCshape': True,'DataTransformation':'X-H','IntrinsicModel':'X-H','trim_on_extras':True,'DF_savekey':'IncLCShape',
													'n_sampling':10000,'n_warmup':2000},
											'label':'$X-H$ Colours w/ LC Shape'},
				###############
				'HighMass'   			:		{   'newdict': {'mass_mode':'high_masses','n_sampling':3000},
											'label':'High-Mass'},
				'LowMass'    			:		{   'newdict': {'mass_mode':'low_masses','n_sampling':3000},
											'label':'Low-Mass'},
				###############
				'Phasemax4'             :       {   'newdict': {'phase_max':4},
											'label':'Data within 4 days of Peak'},
				'Phasemax3'             :       {   'newdict': {'phase_max':3},
											'label':'Data within 3 days of Peak'},
				'Phasemax2'             :       {   'newdict': {'phase_max':2},
											'label':'Data within 2 days of Peak'},
				###############
				'Deviations_NoIntVar'   :		 {    'newdict': {
															'n_sampling':1000,'n_warmup':1000,'include_residuals':False},
											'label':'Deviations no Int. Var.'},
				'AdjCols_NoIntVar'      :       {    'newdict': {'DataTransformation':'Adjacent','IntrinsicModel':'Adjacent',
															'n_sampling':1000,'n_warmup':1000, 'include_residuals':False},
											'label':'Adjacent Colours no Int. Var.'},
				'BXCols_NoIntVar'          :       {    'newdict': {'DataTransformation':'B-X','IntrinsicModel':'B-X',
															'n_sampling':1000,'n_warmup':1000, 'include_residuals':False},
											'label':'$B-X$ Colours no Int. Var.'},
				'XHCols_NoIntVar'          :       {    'newdict': {'DataTransformation':'X-H','IntrinsicModel':'X-H',
															'n_sampling':1000,'n_warmup':1000, 'include_residuals':False},
											'label':'$X-H$ Colours no Int. Var.'},
				################
				}
			}

CYCLE_DICT_RVBins = {'COMMON_CHANGES' : {'newdict':{'BinnedRVFit':True,'n_sampling':2000},'HBMappender':''},

						'RUNS' : {

						'AVExp_Full_AVRVBeta'			:       {    'newdict': {'BVbinboundaries':[],'RVstyles':['AVRVBeta']},
															 		'label':'$\\mathbf{\\beta}$-Model'},
						'AVExp_AVRVBeta_BVcut1.0_Cens'	:       {    'newdict': {'BVbinboundaries':[],'RVstyles':['AVRVBeta'],'CensoredData':True,'BVcutval':1.0},
															 		'label':'$\\mathbf{\\beta}$-Model Cut $B-V>1.0$~mag'},
						'AVExp_AVRVBeta_BVcut0.3_Cens'	:       {    'newdict': {'BVbinboundaries':[],'RVstyles':['AVRVBeta'],'CensoredData':True,'BVcutval':0.3},
															 		'label':'$\\mathbf{\\beta}$-Model Cut $B-V>0.3$~mag'},
						'AVGamma_Full_AVRVBeta'			:       {    'newdict': {'BVbinboundaries':[],'RVstyles':['AVRVBeta'],'AVprior':'Gamma'},
															 		'label':'$\\mathbf{\\beta}$-Model'},
						'AVGamma_AVRVBeta_BVcut1.0_Cens'	:       {    'newdict': {'BVbinboundaries':[],'RVstyles':['AVRVBeta'],'CensoredData':True,'BVcutval':1.0,'AVprior':'Gamma'},
															 		'label':'$\\mathbf{\\beta}$-Model Cut $B-V>1.0$~mag'},
						'AVGamma_AVRVBeta_BVcut0.3_Cens'	:       {    'newdict': {'BVbinboundaries':[],'RVstyles':['AVRVBeta'],'CensoredData':True,'BVcutval':0.3,'AVprior':'Gamma'},
															 		'label':'$\\mathbf{\\beta}$-Model Cut $B-V>0.3$~mag'},
						}
					}
