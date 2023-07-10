import pickle

path = '../../../rvgps/products/stan_fits/FITS/'
files = ['AVGamma_False','AVGamma_True','AVGamma_Cens_inf','AVGamma_Cens_1.0']


for f in files:
    file = f"{path}FIT{f}.pkl"
    print('###'*10)
    print (f)
    with open(file,'rb') as f:
        FIT = pickle.load(f)
    df = FIT['df']
    print ('tauA',df['tauA'].median().round(2),df['tauA'].std().round(2))
    df['nutauA'] = df['nu']*df['tauA']
    print ('nuA*tauA')
    print (f"{df['nutauA'].median().round(2)}\\pm{df['nutauA'].std().round(2)}^{{\,\,{df['nutauA'].quantile(0.95).round(2)}}}_{{\,\,{df['nutauA'].quantile(0.05).round(2)}}}")
