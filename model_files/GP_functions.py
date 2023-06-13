"""
Module contains functions to get GP interpolation using (x,y,yerr) data

Contains:
--------------------
GPFitObj:
    class that stores important GP fit information

get_x_pred:
    function to return grid of x-values where predictions are made based off some xmin and xmax

GP_1D_squaredexp:
    function to interpolate x,y,yerr data with GP 1D squared exponential kernel, outputs GPFitObj class

GP_2D_Matern:
    function to interpolate x,y,yerr and wavelengths with GP 2D Matern kernel, outputs GPFitObj class

--------------------
Classes/Functions take in various arguments

GPFitObj: def __init__(self,x,y,yerr,hyperparameters,gp)

get_x_pred(xmin=None,xmax=None,x=None,Ngrid=10000)

GP_1D_squaredexp(x,y,yerr,x_pred=None,xmin=None,xmax=None,tau_guess=10, return_cov=False)

GP_2D_Matern(x,y,yerr,lambdaC,wavelengths,x_pred=None,xmin=None,xmax=None,tau_guess=10, return_cov=False)

--------------------
Functions use simple operations + george gp. object

Written by Sam M. Ward: smw92@cam.ac.uk
Credits Boone19 https://github.com/kboone/avocado/tree/master
"""

import george
from george import kernels
import numpy as np
from scipy.optimize import minimize
import pandas as pd
import copy

class GPFitObj:
    """
    GPfit class object

    A convenient way of storing GP fit for a single object's light curves

    Parameters
    ----------
    x,y,yerr: each np.array
        the GP interpolations

    hyperparameters: np.array
        the MLE hyperparameters

    gp: the george gp object

    x2: (optional; default=None)
        Nothing yet implemented, may include wavelengths

    Returns
    ---------
    self.df: pandas dataframe
        columns are x,y,yerr
    """
    def __init__(self,x,y,yerr,hyperparameters,gp,x2=None):
        self.x    = x
        self.y    = y
        self.yerr = yerr
        self.df   = pd.DataFrame({'x':x,'y':y,'yerr':yerr})
        self.hyperparameters = hyperparameters
        self.gp   = gp
        if x2 is not None:
            self.x2   = x2


def get_x_pred(xmin=None,xmax=None,x=None,Ngrid=10000):
    """
    Get Predicted x values

    Get grid of x values where GP will predict y values

    Parameters
    ----------
    xmin: (float or int) or None (optional; default=None)
        minimum possible x value to interpolate to (e.g. as defined by model phase range)

    xmax: (float or int) or None (optional; default=None)
        maximum possible x value to interpolate to

    x: (list or array) or None (optional; default=None)
        x-values of data

    Ngrid: int (optional; default=None)
        Optional Parameter (Default: 10,000)
        Number of x grid points to evaluate at

    Returns
    ----------
        x_pred: np.array
            grid of x where GP will interpolate at
    """
    if xmin is None and x is not None:
        xmin = min(x)
    if xmax is None and x is not None:
        xmax = max(x)

    if xmin is not None and xmax is not None:
        x_pred = np.linspace(xmin, xmax, Ngrid)
        return x_pred
    else:
        raise Exception('Need more x information to compute grid of x-values where GP predictions will be made')


def GP_1D_squaredexp(x,y,yerr,x_pred=None,xmin=None,xmax=None,tau_guess=10, return_cov=False):
    """
    Gaussian Process 1-D Squared Exponential

    Compute a Gaussian Process interpolation for (x,y,yerr) data by taking MLE hyperparameters:
    amplitude and length/time-scale

    Parameters
    ----------
    x: list or array
        x-values of data (should have zero or negligible error)

    y: list or array
        y-values of data that we wish to interpolate

    yerr: list or array
        68% uncertainties on y-values

    x_pred: list or array (optional; default=None)
        Grid of x-values where GP will interpolate to

    xmin: float or int or None (optional; default=None)
        minimimum allowed x-value

    xmax: float or int or None (optional; default=None)
        maximum allowed x-value

    tau_guess: float or int (optional; default=10; e.g. 10~d for typical variations)
        Used to initialise time-scale (amplitude is initialised at Std.Dev. of y-data)

    return_cov: bool (optional; default=False)
        if True, return covariance matrix of GP predictions

    Returns
    ----------
    GPfit: GPFitObj class
        contains interpolated y,yerr at x_pred grid, also the george gp object and the MLE hyperparameters
    """

    def neg_ln_like(p):
        #Calculates and returns GP negative log likelihood using current vector of hyperparameters
        gp.set_parameter_vector(p)

        return -gp.log_likelihood(y)

    def grad_neg_ln_like(p):
        #Calculates and returns gradient of GP negative log likelihood using current vector of hyperparameters
        gp.set_parameter_vector(p)

        return -gp.grad_log_likelihood(y)


    #initialise kernel
    kernel = np.mean(y)**2+np.var(y) * kernels.ExpSquaredKernel(tau_guess**2)

    #create GP george object
    gp = george.GP(kernel)
    gp.compute(x, yerr)

    #optimise log likelihood by varying hyperparameters;
    try:                                                                                                                          #e^(1/2)==1.649day minimum for tau
        result = minimize(neg_ln_like, gp.get_parameter_vector(), jac=grad_neg_ln_like, bounds = [(-np.inf,np.inf),(-np.inf,np.inf),(1,np.inf)])
        gp.set_parameter_vector(result.x)

        #k1:log_constant --------> The ln(mean^2)--->Transform via e(Hyp/2)
        #k2:k1:log_constant------> The ln(A^2)------>Transform via e(Hyp/2)
        #k2:k2:metric:log_M_0_0--> The ln(tau^2)---->Transform via e(Hyp/2)

        if x_pred is None:
            #get grid of x-values where GP predictions will be made
            x_pred = get_x_pred(xmin=xmin,xmax=xmax,x=x)
        else:
            pass

        #use george to predict y conditioned on data
        if return_cov:
            y_pred, y_pred_disp = gp.predict(y, x_pred, return_cov=True)
        else:
            y_pred, y_pred_var = gp.predict(y, x_pred, return_var=True)
            y_pred_disp        = y_pred_var**0.5

        GPfit = GPFitObj(np.asarray(x_pred), np.asarray(y_pred), np.asarray(y_pred_disp), np.asarray(result.x), gp)
        return GPfit

    except Exception as e:
        print (f'1DGP Fit Failed because: {e}')
        return None

def GP_2D_Matern(x,y,yerr,lambdaC,wavelengths,x_pred=None,xmin=None,xmax=None,tau_guess=10, return_cov=False):
    """
    Gaussian Process 2-D Matern

    Compute a Gaussian Process interpolation for (x,y,yerr) data by taking MLE hyperparameters:
    const, amplitude and length/time-scale (fix lambda scale at 6000A following Boone19 )

    This function is based on code from https://github.com/kboone/avocado/blob/master/avocado/astronomical_object.py

    Parameters
    ----------
    x: list or array
        x-values of data (should have zero or negligible error)

    y: list or array
        y-values of data that we wish to interpolate

    yerr: list or array
        68% uncertainties on y-values

    lambdaC: list or array
        central wavelengths of each passband, length N_unique_passbands

    wavelengths: list or array
        length of data vectors, contains the central wavelength value

    x_pred: list or array (optional; default=None)
        Grid of x-values where GP will interpolate to

    xmin: float or int or None (optional; default=None)
        minimimum allowed x-value

    xmax: float or int or None (optional; default=None)
        maximum allowed x-value

    tau_guess: float or int (optional; default=10; e.g. 10~d for typical variations)
        Used to initialise time-scale (amplitude is initialised at Std.Dev. of y-data)

    return_cov: bool (optional; default=False)
        if True, return covariance matrix of GP predictions


    Returns
    ----------
    mini_FIT: dict
        each key is the central wavelength, and each item is a GPfit object for GP interpolations at that wavelength;

    So each value is:
    GPfit: GPFitObj class
        contains interpolated y,yerr at x_pred grid, also the george gp object and the MLE hyperparameters
    """

    def neg_ln_like(p):
        #Calculates and returns GP negative log likelihood using current vector of hyperparameters
        gp.set_parameter_vector(p)

        return -gp.log_likelihood(y)

    def grad_neg_ln_like(p):
        #Calculates and returns gradient of GP negative log likelihood using current vector of hyperparameters
        gp.set_parameter_vector(p)

        return -gp.grad_log_likelihood(y)

    if x_pred is None:
        #get grid of x-values where GP predictions will be made
        x_pred = get_x_pred(xmin=xmin,xmax=xmax,x=x)
    else:
        pass

    #Number of grid points in GP time axis
    Ngridt      = copy.deepcopy(len(x_pred))

    #For each time point in grid, have a central_wavelength grid point
    lam_pred    = np.array([ll for ll in lambdaC for _ in range(len(x_pred))])

    #Combine time and wavelength to have a 2D grid
    x_pred      = np.hstack((x_pred for _ in range(len(lambdaC))))
    new_xpred   = np.vstack([x_pred,lam_pred]).T

    #Initialise Kernel using Mean + Matern32Kernel
    #Hyperparameters are:
             #Constant      #Amplitude of variability        #Time Scale Hyperparam   #Lambda Length scale at 6000A
    kernel = np.mean(y)**2+ np.var(y)*kernels.Matern32Kernel([tau_guess ** 2,          6000 ** 2], ndim=2)

    #Freeze Lambda Length Scale Hyperparameter
    kernel.freeze_parameter("k2:k2:metric:log_M_1_1")

    #Create George GP instance and load in data
    gp = george.GP(kernel)
    guess_parameters = gp.get_parameter_vector()
    x_data           = np.vstack([x, wavelengths]).T
    gp.compute(x_data, yerr)


    #Set Bounds on Hyperparameters
    bounds = [  (np.log( (np.mean(y)-np.std(y)*10)**2),  np.log( (np.mean(y)+np.std(y)*10)**2) ),   #Mean between within 10 Std. Devs. of Data mean
                (guess_parameters[1] - 10, guess_parameters[1] + 10),   #Amplitude Param +/-10 in 2log(A) space
                (0, np.log(1000 ** 2))  #Time scale somewhere between 0 and 1000 days
                ]

    #MLE estimate of Hyperparameters
    try:
        result = minimize(neg_ln_like, gp.get_parameter_vector(), jac=grad_neg_ln_like, bounds=bounds)
        gp.set_parameter_vector(result.x)

        #Mean;                     Hyperparameter is ln(mean^2)--> transformation is e^(Hyp/2); k1:log_constant
        #Amplitude of Variability; Hyperparmeter is ln(A^2)    --> tranformation is e^(Hyp/2);  k2:k1:log_constant
        #Timescale of Variability; Hyperparameter is ln(tau)   --> transformation is e^(Hyp/2); k2:k2:metric:log_M_0_0

        #predict y conditioned on data
        if return_cov:
            y_pred, y_pred_disp = gp.predict(y, new_xpred, return_cov=True)
        else:
            y_pred, y_pred_var = gp.predict(y, new_xpred, return_var=True)
            y_pred_disp        = y_pred_var**0.5

        #Split up total fit into respective passbands
        mini_FIT = {}
        for il,ll in enumerate(lambdaC):
            xp,yp,yperr  = x_pred[il*Ngridt:(il+1)*Ngridt],y_pred[il*Ngridt:(il+1)*Ngridt],y_pred_disp[il*Ngridt:(il+1)*Ngridt]
            GPfit        = GPFitObj(np.asarray(xp), np.asarray(yp), np.asarray(yperr), np.asarray(result.x), gp)
            mini_FIT[ll] = GPfit
        return mini_FIT

    except Exception as e:
        print (f'2DGP Fit Failed because: {e}')
        mini_FIT = {}
        for il,ll in enumerate(lambdaC):
            mini_FIT[ll] = None
        return mini_FIT
