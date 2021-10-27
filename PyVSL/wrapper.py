import os
import numpy as np
try:
    # rpy2
    from rpy2.robjects.packages import importr
    import rpy2.robjects.numpy2ri
    import rpy2.robjects as ro
    rpy2.robjects.numpy2ri.activate()
except ImportError:
    pass

import oct2py


def VSL_R(syear, eyear, phi, T, P, T1=8, T2=23, M1=0.01, M2=0.05, Mmax=0.76, Mmin=0.01,
          alph=0.093, m_th=4.886, mu_th=5.8, rootd=1000, M0=0.2, substep=0, I_0=1, I_f=12,
          hydroclim='P', Rlib_path=None):

    ''' VS-Lite tree-ring PSM

    Porting the VSLiteR code (https://github.com/fzhu2e/VSLiteR),
    which is a fork of the original version by Suz Tolwinski-Ward (https://github.com/suztolwinskiward/VSLiteR)
    NOTE: Need to install the R library following the below steps:
        1. install.packages("devtools")
        2. library(devtools)
        3. install_github("fzhu2e/VSLiteR")

    Args:
        syear (int): Start year of simulation.
        eyear (int): End year of simulation.
        phi (float): Latitude of site (in degrees N).
        T (array): temperature timeseries in deg C, with length of 12*Nyrs
        P (array): precipitation timeseries in mm/month, with length of 12*Nyrs
        # original param. in R package: T (12 x Nyrs) Matrix of ordered mean monthly temperatures (in degEes C).
        # original param. in R package: P (12 x Nyrs) Matrix of ordered accumulated monthly precipitation (in mm).
        T1: Lower temperature threshold for growth to begin (scalar, deg. C).
        T2: Upper temperature threshold for growth sensitivity to temp (scalar, deg. C).
        M1: Lower moisture threshold for growth to begin (scalar, v.v).
        M2: Upper moisture threshold for growth sensitivity to moisture (scalar, v/v).
        Mmax: Scalar maximum soil moisture held by the soil (in v/v).
        Mmin: Scalar minimum soil moisture (for error-catching) (in v/v).
        alph: Scalar runoff parameter 1 (in inverse months).
        m_th: Scalar runoff parameter 3 (unitless).
        mu_th: Scalar runoff parameter 2 (unitless).
        rootd: Scalar root/"bucket" depth (in mm).
        M0: Initial value for previous month's soil moisture at t = 1 (in v/v).
        substep: Use leaky bucket code with sub-monthly time-stepping? (TRUE/FALSE)
        I_0: lower bound of integration window (months before January in NH)
        I_f: upper bound of integration window (months after January in NH)
        hydroclim: Switch; value is either "P" (default) or "M" depending on whether the second input climate variable
            is precipitation, in which case soil moisture isestimated using the Leaky Bucket model of the CPC,
            or soil moisture, in which case the inputs are used directly to compute the growth response.

    Returns:
        res (dict): a dictionary with several lists, keys including:
            trw: tree ring width
            gT, gM, gE, M, potEv, sample.mean.width, sample.std.width

    Reference:
        + The original VSLiteR code by Suz Tolwinski-Ward (https://github.com/suztolwinskiward/VSLiteR)
        + Tolwinski-Ward, S.E., M.N. Evans, M.K. Hughes, K.J. Anchukaitis, An efficient forward model of the climate controls
            on interannual variation in tree-ring width, Climate Dynamics, doi:10.1007/s00382-010-0945-5, (2011).
        + Tolwinski-Ward, S.E., K.J. Anchukaitis, M.N. Evans, Bayesian parameter estimation and interpretation for an
            intermediate model of tree-ring width , Climate of the Past, doi:10.5194/cpd-9-615-2013, (2013).
        + Tolwinski-Ward, S.E., M.P. Tingley, M.N. Evans, M.K. Hughes, D.W. Nychka, Probabilistic reconstructions of localtemperature and
            soil moisture from tree-ring data with potentially time-varying climatic response, Climate Dynamics, doi:10.1007/s00382-014-2139-z, (2014).
    '''
    if Rlib_path is None:
        Rlib_path = '/Library/Frameworks/R.framework/Versions/Current/Resources/library'  # on macOS

    ro.r('.libPaths("{}")'.format(Rlib_path))
    VSLiteR = importr('VSLiteR')
    nyr = eyear - syear + 1
    T_model = T.reshape((nyr, 12)).T
    P_model = P.reshape((nyr, 12)).T

    res = VSLiteR.VSLite(syear, eyear, phi, T_model, P_model, T1=T1, T2=T2, M1=M1, M2=M2, Mmax=Mmax, Mmin=Mmin,
                         alph=alph, m_th=m_th, mu_th=mu_th, rootd=rootd, M0=M0, substep=substep,
                         I_0=I_0, I_f=I_f, hydroclim=hydroclim)

    res_dict = dict(zip(res.names, map(list, list(res))))
    res_dict['trw'] = res_dict['trw'][0]


    return res_dict

def VSL_M(syear, eyear, phi, T, P, T1=8, T2=23, M1=0.01, M2=0.05, Mmax=0.76, Mmin=0.01,
          alph=0.093, m_th=4.886, mu_th=5.8, rootd=1000, M0=0.2, substep=0, I_0=1, I_f=12,
          hydroclim="P"):

    ''' VS-Lite tree-ring PSM

    Porting the VSLite Matlab code by Suz Tolwinski-Ward (https://github.com/suztolwinskiward/VSLite/blob/master/VSLite_v2_3.m)

    NOTE: Need to install Octave as a substitute of Matlab

    Args:
        syear (int): Start year of simulation.
        eyear (int): End year of simulation.
        phi (float): Latitude of site (in degrees N).
        T (array): temperature timeseries in deg C, with length of 12*Nyrs
        P (array): precipitation timeseries in mm/month, with length of 12*Nyrs
        # original param. in R package: T (12 x Nyrs) Matrix of ordered mean monthly temperatures (in degEes C).
        # original param. in R package: P (12 x Nyrs) Matrix of ordered accumulated monthly precipitation (in mm).
        T1: Lower temperature threshold for growth to begin (scalar, deg. C).
        T2: Upper temperature threshold for growth sensitivity to temp (scalar, deg. C).
        M1: Lower moisture threshold for growth to begin (scalar, v.v).
        M2: Upper moisture threshold for growth sensitivity to moisture (scalar, v/v).
        Mmax: Scalar maximum soil moisture held by the soil (in v/v).
        Mmin: Scalar minimum soil moisture (for error-catching) (in v/v).
        alph: Scalar runoff parameter 1 (in inverse months).
        m_th: Scalar runoff parameter 3 (unitless).
        mu_th: Scalar runoff parameter 2 (unitless).
        rootd: Scalar root/"bucket" depth (in mm).
        M0: Initial value for previous month's soil moisture at t = 1 (in v/v).
        substep: Use leaky bucket code with sub-monthly time-stepping? (TRUE/FALSE)
        I_0: lower bound of integration window (months before January in NH)
        I_f: upper bound of integration window (months after January in NH)
        hydroclim: Switch; value is either "P" (default) or "M" depending on whether the second input climate variable
            is precipitation, in which case soil moisture isestimated using the Leaky Bucket model of the CPC,
            or soil moisture, in which case the inputs are used directly to compute the growth response.

    Returns:
        res (dict): a dictionary with several lists, keys including:
            trw: tree ring width
            gT, gM, gE, M, potEv, sample.mean.width, sample.std.width

    Reference:
        + The original VSLite Matlab code by Suz Tolwinski-Ward (https://github.com/suztolwinskiward/VSLite/blob/master/VSLite_v2_3.m)
        + Tolwinski-Ward, S.E., M.N. Evans, M.K. Hughes, K.J. Anchukaitis, An efficient forward model of the climate controls
            on interannual variation in tree-ring width, Climate Dynamics, doi:10.1007/s00382-010-0945-5, (2011).
        + Tolwinski-Ward, S.E., K.J. Anchukaitis, M.N. Evans, Bayesian parameter estimation and interpretation for an
            intermediate model of tree-ring width , Climate of the Past, doi:10.5194/cpd-9-615-2013, (2013).
        + Tolwinski-Ward, S.E., M.P. Tingley, M.N. Evans, M.K. Hughes, D.W. Nychka, Probabilistic reconstructions of localtemperature and
            soil moisture from tree-ring data with potentially time-varying climatic response, Climate Dynamics, doi:10.1007/s00382-014-2139-z, (2014).
    '''
    dirpath = os.path.dirname(__file__)
    oc = oct2py.Oct2Py(temp_dir=dirpath)
    oc.addpath(dirpath)

    nyr = eyear - syear + 1
    T_model = T.reshape((nyr, 12)).T
    P_model = P.reshape((nyr, 12)).T

    trw, gT, gM, gE, M, potEv, width, width_mean, width_std = oc.feval('VSLite_v2_3',
        syear, eyear, phi,  T1, T2, M1, M2, T_model, P_model, Mmax=Mmax, Mmin=Mmin,
        alph=alph, m_th=m_th, mu_th=mu_th, rootd=rootd, M0=M0, substep=substep,
        I_0=I_0, I_f=I_f, hydroclim=hydroclim, nout=9
    )

    res_dict = {
        'trw': trw[0],
        'gT': gT,
        'gM': gM,
        'gE': gE,
        'M': M,
        'potEv': potEv,
        'width': width,
        'width_mean': width_mean,
        'width_std': width_std,
    }

    return res_dict

def leakybucket_monthly(syear, eyear, phi, T, P, Mmax=0.76, Mmin=0.01, alph=0.093, m_th=4.886, mu_th=5.8, rootd=1000, M0=0.2):
    ''' Leaky bucket model: simulate soil moisture with coarse monthly time step, outputs simulated soil moisture and potential evapotranspiration.

    Args:
        syear (int): Start year of simulation.
        eyear (int): End year of simulation.
        phi (float): Latitude of site (in degrees N).
        T (array): temperature timeseries in deg C, with length of 12*Nyrs
        P (array): precipitation timeseries in mm/month, with length of 12*Nyrs
        # original param. in R package: T (12 x Nyrs) Matrix of ordered mean monthly temperatures (in degEes C).
        # original param. in R package: P (12 x Nyrs) Matrix of ordered accumulated monthly precipitation (in mm).
        Mmax: Scalar maximum soil moisture held by the soil (in v/v).
        Mmin: Scalar minimum soil moisture (for error-catching) (in v/v).
        alph: Scalar runoff parameter 1 (in inverse months).
        m_th: Scalar runoff parameter 3 (unitless).
        mu_th: Scalar runoff parameter 2 (unitless).
        rootd: Scalar root/"bucket" depth (in mm).
        M0: Initial value for previous month's soil moisture at t = 1 (in v/v).

    Returns:
        M (array): soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs)
        potEv (array): potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm)
    '''

    dirpath = os.path.dirname(__file__)
    oc = oct2py.Oct2Py(temp_dir=dirpath)
    oc.addpath(dirpath)

    nyr = eyear - syear + 1
    T_model = T.reshape((nyr, 12)).T
    P_model = P.reshape((nyr, 12)).T

    M, potEv, ndl, cdays = oc.feval('leakybucket_monthly', syear, eyear, phi,  T_model, P_model, Mmax, Mmin, alph, m_th, mu_th, rootd, M0, nout=4)

    res_dict = {
        'M': M,
        'potEv': potEv,
        'ndl': ndl,
        'cdays': cdays,
    }

    return res_dict


def est_params(
    T, P, lat, TRW, nyr=None,
    nsamp=1000, errormod=0, gparpriors='fourbet',
    pt_ests='med', seed=0,
    beta_params=np.matrix([
        [9, 5, 0, 9],
        [3.5, 3.5, 10, 24],
        [1.5, 2.8, 0, 0.1],
        [1.5, 2.5, 0.1, 0.5],
    ])):
    ''' Run the VSL parameter estimatino Matlab precedure
    Args:
        T (1-D array): monthly surface air temperature [degC]
        P (1-D array): monthly accumulated precipitation [mm]
        lat (float): latitude of the site
        pt_ests (str): 'med' or 'mle'
        nsamp (int): the number of MCMC iterations
        errmod (int): 0: white noise, 1: AR(1) noise
        gparpriors (str): 'fourbet': beta distribution, 'uniform': uniform distribution
        beta_params (matrix): the beta distribution parameters for T1, T2, M1, M2
        seed (int): random seed
    '''
    dirpath = os.path.dirname(__file__)
    oc = oct2py.Oct2Py(temp_dir=dirpath)
    oc.addpath(dirpath)

    if nyr is None:
        nyr = np.size(T) // 12

    T_model = T.reshape((nyr, 12)).T
    P_model = P.reshape((nyr, 12)).T

    # T1, T2, M1, M2 = oct2py.octave.feval(
    T1, T2, M1, M2 = oc.feval(
        'estimate_vslite_params_v2_3',
        T_model, P_model, lat, TRW,
        'seed', seed, 'nsamp', nsamp, 'errormod', errormod,
        'gparpriors', gparpriors, 'fourbetparams', beta_params,
        'pt_ests', pt_ests, nargout=4, nout=4,
    )

    res_dict = {
        'T1': T1,
        'T2': T2,
        'M1': M1,
        'M2': M2,
    }

    return res_dict
