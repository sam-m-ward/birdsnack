data {
	int S;                         //Total Number of SNe (Censored + Uncensored)
	int SC;                        //Number of Censored SNe
	int<lower=1> Nm;               //Number of mags per SN

	vector[(S-SC)*Nm] mags;        //Vector of apparent magnitude point estimates
	vector[(S-SC)*Nm] mags_errs;   //Vector of apparent magnitude measurement errors
	vector[SC] BVerrs_Cens;        //Vector of B-V measurements errors on censored SNe

	vector[9] xk;                  //Fitzpatrick99 knot wavelengths
	matrix [Nm, 9] Mmatrix;        //Fitzpatrick99 knots->lambda

	real muRVmin;                  //Lower bound on muRV prior
	real muRVmax;                  //Upper bound on muRV prior
	real RVsmin;                   //Lower Truncation on RVs parameters
	real disp_sigmaRV;             //Dispersion on Half-Normal prior on sigma_RV, default is 2
	real disp_sigmaAV;             //Dispersion on Half-Normal prior on sigma_AV, default is 1

	real a_sigma_mint;             //Dispersion of Half-Cauchy Hyperprior on Chromatic Intrinsic Dispersion Hyperparameters

	vector[S] dm15Bs;              //Light curve shape point estimates,    if gamma_shape=0, then these are dummy values
	vector[S] dm15B_errs;          //Light curve shape measurement errors, if gamma_shape=0, then these are dummy values
	real gamma_shape;              //If 1, include LC shape term, if 0, exclude LC shape term
	real gamma_res;                //If 1, include intrinsic colour residuals, if 0, exclude residuals

	real zero_index;               //Index denoting reference band that is set to zero in mu_int
	real BVcutval;                 //The upper BV value used for cutting sample, default is often 0.3~mag, this value only matters when SC>0

	matrix[Nm,Nm] L_mint;              //Cholesky decomposition of intrinsic magnitudes
	vector[Nm] FPC0;                  //The zeroth intrinsic mag functional principal component (same for each SN)
	vector[Nm] FPC1;                  //The first intrinsic mag functional principal component (same for each SN)
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
	vector<lower=0,upper=100>[S-SC] mus;              //Distance Parameters
	vector[(S-SC)*Nm] eta_mint;                       //Re-parameterised residual intrinsic deviations
	vector[SC*2] eta_mint_Cens;                       //Re-parameterised residual intrinsic deviations of Censored SNe (only BV saves on No. of parameters)
	vector<lower=BVcutval>[SC] BVs_Cens;              //The censored data which have B-V>0.3
	vector<lower=0,upper=1>[S] nuAVs;                 //Dust extinction parameters
	//Re-parameterised Dust Parameters
	vector<lower=0,upper=1>[S] nuRVs;                 //Re-parameterised Individual Dust-Law Shape Parameters
	vector[S] dm15B_latent;                           //Light curve shape parameters

	//Hyperparameters
	real<lower=muRVmin,upper=muRVmax> mu_RV;          //RV population mean hyperparameter
	real<lower=0> eta_sig_RV;                         //Transform of RV population dispersion hyperparameter
	real<lower=0> eta_sig_AV;                         //Transform of AV population dispersion hyperparameter
	real<lower=0, upper=pi()/2> tauA_tform;           //Transform of tau_A, AV expontential dist. hyperparameter
	//Transformed Gamma AVs distribution shape hyperparameter
	//real<lower=0,upper=pi()/2> nu_tform;
}

transformed parameters {
	matrix[Nm,S-SC] mags_matrix;       //Matrix of latent apparent magnitudes
	vector[(S-SC)*Nm] mags_latent;     //Vector of latent apparent magnitudes
	matrix[2,SC] mags_matrix_Cens;     //Matrix of latent apparent magnitudes for censored SNe
	vector[SC] BVslatent_Cens;         //Vector of latent BV colours of censored SNe

	//Rescaled Parameters
	vector<lower=0>[S] AVs;            //Individual Dust Extinction Parameters
	vector<lower=RVsmin>[S] RVs;       //Individual Dust-Law Shape Parameters
	vector[S*Nm] xivec;               //Dust law deviations

	//Rescaled Hyperparameters
	real<lower=0> tauA;                //Dust extinction Hyperparameter
	real<lower=0> sig_AV;              //Dust extinction population dispersion Hyperparameter
	real<lower=0> sig_RV;              //RV population dispersion Hyperparameter

	matrix[2,2] L_mint_Cens;           //Cholesky decomposition of BV intrinsic magnitudes

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
	sig_AV  = eta_sig_AV*disp_sigmaAV;
	sig_RV  = eta_sig_RV*disp_sigmaRV;
	//nu      = tan(nu_tform);

	//Tranform nuRVs->RVs via truncated Gaussian (RVsmin=RVsmin, RVsmax=inf)
	alpha  = (RVsmin - mu_RV)/sig_RV;
	RVs    = mu_RV - sig_RV * ( inv_Phi ( nuRVs * Phi (-alpha) ) );
	AVs    = tauA - sig_AV * ( inv_Phi ( nuAVs * Phi (tauA/sig_AV) ) );
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

	//Dust Reddening is Fitzpatrick99 evaluated at central wavelength
	for (s in 1:S-SC) {
		xivec[Nm*(s-1)+1:Nm*s] = Mmatrix * to_vector(yk[s,:])/RVs[s]; //No Mean Subtraction
		mags_matrix[:,s]       = mus[s] + AVs[s]*xivec[Nm*(s-1)+1:Nm*s] + FPC0 + gamma_res * L_mint * eta_mint[Nm*(s-1)+1:Nm*s] + gamma_shape*(dm15B_latent[s]-1.05)*FPC1; //Tripp98 offset is 1.05
	}
	//Apparent Magnitudes are combination of Intrinsic Abs Mags and Dust Extinction and Distance
	mags_latent = to_vector(mags_matrix);

	if (SC>0) {
		L_mint_Cens = L_mint[1:2,1:2];
		for (s in S-SC+1:S) {
			xivec[Nm*(s-1)+1:Nm*s]   = (Mmatrix * to_vector(yk[s,:])) / RVs[s];
			mags_matrix_Cens[:,s-sp] = AVs[s]*xivec[Nm*(s-1)+1:Nm*(s-1)+2] + FPC0[1:2] + gamma_res*L_mint_Cens*eta_mint_Cens[2*(s-sp-1)+1:2*(s-sp-1)+2] + gamma_shape*(dm15B_latent[s]-1.05)*FPC1[1:2]; //Tripp98 offset is 1.05
			BVslatent_Cens[s-sp]     = mags_matrix_Cens[1,s-sp]-mags_matrix_Cens[2,s-sp];
		}
	}
}

model {
	//Priors on Parameters
	mus           ~ uniform(0,100);
	eta_mint      ~ std_normal();
	eta_mint_Cens ~ std_normal();
	nuRVs         ~ uniform(0,1);
	nuAVs         ~ uniform(0,1);
	//AVs         ~ gamma(nu,1/tauA);

	//Priors on Hyperparameters
	mu_RV          ~ uniform(muRVmin,muRVmax);
	eta_sig_RV     ~ std_normal();
	eta_sig_AV		 ~ std_normal();
	tauA_tform     ~ uniform(0,pi()/2);
	//nu_tform     ~ uniform(0,pi()/2);

	//Likelihood
	dm15Bs   ~ normal(dm15B_latent,dm15B_errs);
	mags     ~ normal(mags_latent, mags_errs);
	BVs_Cens ~ normal(BVslatent_Cens, BVerrs_Cens);
}
