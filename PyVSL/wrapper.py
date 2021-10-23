# rpy2
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
import rpy2.robjects as ro
rpy2.robjects.numpy2ri.activate()


def VSL_R(syear, eyear, phi, T, P, T1=8, T2=23, M1=0.01, M2=0.05, Mmax=0.76, Mmin=0.01,
          alph=0.093, m_th=4.886, mu_th=5.8, rootd=1000, M0=0.2, substep=0, I_0=1, I_f=12,
          hydroclim="P", Rlib_path=None):

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
        T (array): temperature timeseries in K, with length of 12*Nyrs
        P (array): precipitation timeseries in kg/m2/s, with length of 12*Nyrs
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
    T_model = T.reshape((nyr, 12)).T - 273.15  # in monthly temperatures (degC)
    P_model = P.reshape((nyr, 12)).T * 3600*24*30  # in accumulated monthly precipitation (mm)

    res = VSLiteR.VSLite(syear, eyear, phi, T_model, P_model, T1=T1, T2=T2, M1=M1, M2=M2, Mmax=Mmax, Mmin=Mmin,
                         alph=alph, m_th=m_th, mu_th=mu_th, rootd=rootd, M0=M0, substep=substep,
                         I_0=I_0, I_f=I_f, hydroclim=hydroclim)

    res = dict(zip(res.names, map(list, list(res))))

    return res