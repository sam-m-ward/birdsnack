class HBM_preprocessor:

    def __init__(self,choices):
        self.choices = choices


    def get_leff_rest():

        l_eff_rest  = np.array([fset[pb].ave_wave for pb in pblist])

		if 'effective' in lam_choice:
			if lam_choice=='effective_thetaON':
				LAM = pd.read_csv('sbc/DLAMS/DLAM_thetaON.txt',index_col=False)[fclist].median(axis=0).values
				print (f'Changing effective wavelength (thetaON): {LAM}')
			if lam_choice=='effective_thetaOFF':
				LAM = pd.read_csv('sbc/DLAMS/DLAM_thetaOFF.txt',index_col=False)[fclist].median(axis=0).values
				print (f'Changing effective wavelength (thetaOFF): {LAM}')
			if lam_choice=='effective_thetaOFFepsOFF':
				LAM = pd.read_csv('sbc/DLAMS/DLAM_thetaOFFepsOFF.txt',index_col=False)[fclist].median(axis=0).values
				print (f'Changing effective wavelength (thetaOFFepsOFF): {LAM}')
			if lam_choice=='effective_taustar':
				if BVcut is False or not CensoredData:
					LAM = pd.read_csv('sbc/DLAMS/DLAM_tau0.45.txt',index_col=False)[fclist].median(axis=0).values
					print (f'Changing effective wavelength (tau*=0.45): {LAM}')
				elif BVcut is True or CensoredData:
					LAM = pd.read_csv('sbc/DLAMS/DLAM_tau0.32.txt',index_col=False)[fclist].median(axis=0).values
					print (f'Changing effective wavelength (tau*=0.32): {LAM}')
			l_eff_rest = LAM

        return l_eff_rest

    def get_dustlaw():
        if dustlaw=='fitzpatrick99':
			xk   = np.array([0.0, 1e4/26500., 1e4/12200., 1e4/6000., 1e4/5470., 1e4/4670., 1e4/4110., 1e4/2700., 1e4/2600.])
			KD_x = spline_utils.invKD_irr(xk)
			M_fitz_block  = spline_utils.spline_coeffs_irr(1e4/l_eff_rest, xk, KD_x)
			if ColourChoice=='Adjacent':#Adjacent Colours
				dM_fitz_block   = np.array([ M_fitz_block[i,:]-M_fitz_block[i+1,:] for i in range(M_fitz_block.shape[0]-1) ])
			elif ColourChoice=='B-X':#B-X colours
				dM_fitz_block   = np.array([ M_fitz_block[0,:]-M_fitz_block[i+1,:] for i in range(M_fitz_block.shape[0]-1) ])
			elif ColourChoice=='X-H':#X-H colours
				dM_fitz_block   = np.array([ M_fitz_block[i,:]-M_fitz_block[-1,:]  for i in range(M_fitz_block.shape[0]-1) ])
		else:
			import extinction
			RV1,RV2 = 2,3#Dummy RV values used to back out the a,b values
			_DRV    = (1/RV1 - 1/RV2)
			if dustlaw=='ccm89':
				b_vec   = (extinction.ccm89(l_eff_rest,1,RV1)-extinction.ccm89(l_eff_rest,1,RV2))/_DRV
				a_vec   =  extinction.ccm89(l_eff_rest,1,RV1)-b_vec/RV1
			elif dustlaw=='ccm89_o94':
				b_vec   = (extinction.odonnell94(l_eff_rest,1,RV1)-extinction.odonnell94(l_eff_rest,1,RV2))/_DRV
				a_vec   =  extinction.odonnell94(l_eff_rest,1,RV1)-b_vec/RV1
			Del_a_vec = a_vec[:-1]-a_vec[1:]
			Del_b_vec = b_vec[:-1]-b_vec[1:]
