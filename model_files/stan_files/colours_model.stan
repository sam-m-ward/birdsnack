data {
	int S;                         //Total Number of SNe (Censored + Uncensored)
	int SC;                        //Number of Censored SNe
	int<lower=1> Nc;               //Number of colours per SN

	vector[Nc*(S-SC)] capps;               //Vector of apparent colours
	array[S-SC] matrix[Nc,Nc] capps_errs;  //Vector of matrices, each matrix is covariance of colour measurement errors
	vector[SC] BVerrs_Cens;        		     //Vector of B-V measurements errors on censored SNe

	vector[9] xk;                  //Fitzpatrick99 knot wavelengths
	matrix [Nc, 9] DelM;        //Fitzpatrick99 knots->lambda

	real muRVmin;                  //Lower bound on muRV prior
	real muRVmax;                  //Upper bound on muRV prior
	real RVsmin;                   //Lower Truncation on RVs parameters
	real disp_sigmaRV;             //Dispersion on Half-Normal prior on sigma_RV, default is 2

	real a_sigma_cint;             //Dispersion of Half-Cauchy Hyperprior on Chromatic Intrinsic Dispersion Hyperparameters

	vector[S] dm15Bs;              //Light curve shape point estimates,    if gamma_shape=0, then these are dummy values
	vector[S] dm15B_errs;          //Light curve shape measurement errors, if gamma_shape=0, then these are dummy values
	real gamma_shape;              //If 1, include LC shape term, if 0, exclude LC shape term
	real gamma_res;                //If 1, include intrinsic colour residuals, if 0, exclude residuals

	matrix[Nc,Nc] CC;              //Matrix to transform from 1 set of reference colours to another
	matrix[Nc,Nc] CC_to_adj;       //Matrix to transform from 1 set of reference colours to another
}

transformed data {
	//Fitzpatrick 99 constants
	real f99_x0 = 4.596;
	real f99_gamma = 0.99;
	real f99_c3 = 3.23;
	real f99_c4 = 0.41;
	real f99_c5 = 5.9;
	vector[2] f99_d = square(xk[8:9]) ./ (square(square(xk[8:9]) - square(f99_x0)) + square(f99_gamma*xk[8:9]));

	int sp = S-SC;
}

parameters {
	//Parameters
	vector[Nc*SC] eta_cint_Cens;                       //Re-parameterised residual intrinsic deviations of Censored SNe (only BV saves on No. of parameters)
	vector<lower=0.3>[SC] BVs_Cens;                   //The censored data which have B-V>0.3
	vector<lower=0>[S] AVs;                           //Dust extinction parameters
	//Re-parameterised Dust Parameters
	vector<lower=0,upper=1>[S] nuRVs;                 //Re-parameterised Individual Dust-Law Shape Parameters
	vector[S] dm15B_latent;                           //Light curve shape parameters

	//Hyperparameters
	real<lower=muRVmin,upper=muRVmax> mu_RV;          //RV population mean hyperparameter
	real<lower=0> eta_sig_RV;                         //Transform of RV population dispersion hyperparameter
	real<lower=0, upper=pi()/2> tauA_tform;           //Transform of tau_A, AV expontential dist. hyperparameter
	//real<lower=0,upper=pi()/2> nu_tform;            //Gamma AVs distribution shape hyperparameter
	vector[Nc] FPC0;                                  //The zeroth intrinsic colour functional principal component (same for each SN)
	vector[Nc] FPC1;                                  //The first intrinsic colour functional principal component (same for each SN)
	cholesky_factor_corr[Nc] L_cint_eta;              //Cholesky factor of unscaled intrinsic colour covariance matrix
	vector<lower=0,upper=pi()/2>[Nc] sigma_cint_eta;  //Scaling of covariance matrix for each intrinsic colour

}

transformed parameters {
	matrix[Nc,S-SC] capps_matrix;      //Matrix of reddenings and intrinsic colours
	vector[Nc*(S-SC)] capps_latent;    //Vector of latent apparent colours
	vector[Nc*SC] c_latent_Cens;       //Vector of modelling colours for censored SNe
	vector[SC] BVslatent_Cens;         //Vector of latent BV colours of censored SNe


	//Rescaled Parameters
	vector<lower=RVsmin>[S] RVs;       //Individual Dust-Law Shape Parameters

	//Rescaled Hyperparameters
	real<lower=0> tauA;                //Dust extinction Hyperparameter
	real<lower=0> sig_RV;              //RV population dispersion Hyperparameter

	vector<lower=0>[Nc] sigma_cint;    //Vector of intrinsic colour dispersions
	matrix[Nc,Nc] L_cint;              //Cholesky decomposition of intrinsic colours

	real<upper=0> alpha;
	//real<lower=0> beta;
	//real<lower=0> nu;

	//Fitzpatrick99 parameters
	matrix[S,9] yk;
	real f99_c1;
	real f99_c2;
	real RV;


	//Rescaled Dust Extinction Population Distribution Hyperparameters
	tauA    = tan(tauA_tform);
	sig_RV  = eta_sig_RV*disp_sigmaRV;

	//Tranform nuRVs->RVs via truncated Gaussian (RVsmin=RVsmin, RVsmax=inf)
	alpha  = (RVsmin - mu_RV)/sig_RV;
	RVs    = mu_RV - sig_RV * ( inv_Phi ( nuRVs * Phi (-alpha) ) );
	//Upper-Truncated-Normal
	//beta   = (RVsTrunc - mu_RV)/sig_RV;
	//RVs  = mu_RV + sig_RV*(inv_Phi(nuRVs*(Phi(beta)-Phi(alpha)) + Phi(alpha)));

	//Compute Fitzpatrick99 parameters
	for (s in 1:S) {
		RV = RVs[s];
		f99_c2 = -0.824 + 4.717/RV;
		f99_c1 = 2.030 - 3.007*f99_c2;

		yk[s,1] = -RV;
		yk[s,2] =	0.26469*RV/3.1 - RV;
		yk[s,3] = 0.82925*RV/3.1 - RV;
		yk[s,4] = -0.422809 + 1.00270*RV + 2.13572e-4*square(RV) - RV;
		yk[s,5] = -5.13540e-2 + 1.00216*RV - 7.35778e-5*square(RV) - RV;
		yk[s,6] = 0.700127 + 1.00184*RV - 3.32598e-5*square(RV) - RV;
		yk[s,7] = 1.19456 + 1.01707*RV - 5.46959e-3*square(RV) + 7.97809e-4*pow(RV,3) - 4.45636e-5*pow(RV,4) - RV;
		yk[s,8:9] = to_row_vector(f99_c1 + f99_c2*xk[8:9] + f99_c3*f99_d);
	}

	//Re-scaled Covariance Matrix for Intrinsic Colours
	sigma_cint   = a_sigma_cint*tan(sigma_cint_eta);
	L_cint       = diag_pre_multiply(sigma_cint,L_cint_eta);

	//Dust Reddening is Fitzpatrick99 evaluated at central wavelength
	for (s in 1:S-SC) {
		capps_matrix[:,s]  = AVs[s] * DelM * to_vector(yk[s,:]) / RVs[s] + FPC0 + gamma_shape*(dm15B_latent[s]-1.05)*FPC1;
	}
	//Apparent Colours are combination of Intrinsic Colour and Dust Reddening
  capps_latent = to_vector(capps_matrix);

	if (SC>0) {
		for (s in S-SC+1:S) {
			c_latent_Cens[Nc*(s-sp-1)+1:Nc*(s-sp)] = AVs[s]*DelM*to_vector(yk[s,:])/RVs[s] + FPC0 + gamma_res*L_cint*eta_cint_Cens[Nc*(s-sp-1)+1:Nc*(s-sp)] + gamma_shape*(dm15B_latent[s]-1.05)*FPC1; //Tripp98 offset is 1.05
			c_latent_Cens[Nc*(s-sp-1)+1:Nc*(s-sp)] = CC_to_adj*c_latent_Cens[Nc*(s-sp-1)+1:Nc*(s-sp)];
			BVslatent_Cens[s-sp] = c_latent_Cens[Nc*(s-sp-1)+1]; //Whatever the colours being modelled, transform to adjacent colours and pick out first entry which is B-V
		}
	}
}

model {
	//Priors on Parameters
	eta_cint_Cens ~ std_normal();
	nuRVs         ~ uniform(0,1);
	AVs           ~ exponential(inv(tauA));
	//AVs         ~ gamma(nu,1/tauA);

	//Priors on Hyperparameters
	mu_RV          ~ uniform(muRVmin,muRVmax);
	eta_sig_RV     ~ std_normal();
	tauA_tform     ~ uniform(0,pi()/2);
	//nu_tform     ~ uniform(0,pi()/2);
	FPC0           ~ normal(0,10);
	FPC1           ~ normal(0,10);
	L_cint_eta     ~ lkj_corr_cholesky(1);
	sigma_cint_eta ~ uniform(0,pi()/2);

	//Likelihood
	dm15Bs   ~ normal(dm15B_latent,dm15B_errs);
	for (s in 1:S-SC) {
		capps[Nc*(s-1)+1:Nc*s] ~ multi_normal(capps_latent[Nc*(s-1)+1:Nc*s], capps_errs[s]+gamma_res*CC*L_cint*L_cint'*CC');
	}
	BVs_Cens ~ normal(BVslatent_Cens, BVerrs_Cens);
}
