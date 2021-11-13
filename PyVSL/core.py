'''
VS-Lite in pure Python
@author: Feng Zhu (fzhu@nuist.edu.cn)
'''
import numpy as np

def VSL(syear, eyear, phi, T, P, T1=8, T2=23, M1=0.01, M2=0.05, Mmax=0.76, Mmin=0.01, sensitivity='TM',
        alph=0.093, m_th=4.886, mu_th=5.8, rootd=1000, M0=0.2, substep=0, I_0=1, I_f=12, hydroclim='P'):
    '''
    Translated from VSLite_v2_3.m - Simulate tree ring width index given monthly climate inputs.

    Basic Usage:
       trw = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P)
       gives just simulated tree ring as ouput.
    
      [trw,gT,gM,gE,M] = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P))
       also includes growth response to temperature, growth response to soil
       moisture, scaled insolation index, and soil moisture estimate in outputs.
    
    Basic Inputs:
      syear = start year of simulation.
      eyear = end year of simulation.
      phi = latitude of site (in degrees N)
      T1 = scalar temperature threshold below which temp. growth response is zero (in deg. C)
      T2 = scalar temperature threshold above which temp. growth response is one (in deg. C)
      M1 = scalar soil moisture threshold below which moist. growth response is zero (in v/v)
      M2 = scalar soil moisture threshold above which moist. growth response is one (in v/v)
        (Note that optimal growth response parameters T1, T2, M1, M2 may be estimated
         using code estimate_vslite_params_v2_3.m also freely available at
         the NOAA NCDC Paleoclimatology software library.)
      T = (12 x Nyrs) matrix of ordered mean monthly temperatures (in degEes C)
      P = (12 x Nyrs) matrix of ordered accumulated monthly precipitation (in mm)
    
    Advanced Inputs (must be specified as property/value pairs):
        'lbparams':  Parameters of the Leaky Bucket model of soil moisture.
                   These may be specified in an 8 x 1 vector in the following
                   order (otherwise the default values are read in):
                      Mmax: scalar maximum soil moisture content (in v/v),
                        default value is 0.76
                      Mmin: scalar minimum soil moisture (in v/v), default
                        value is 0.01
                      alph: scalar runoff parameter 1 (in inverse months),
                        default value is 0.093
                      m_th: scalar runoff parameter 3 (unitless), default
                        value is 4.886
                      mu_th: scalar runoff parameter 2 (unitless), default
                        value is 5.80
                      rootd: scalar root/"bucket" depth (in mm), default
                        value is 1000
                      M0: initial value for previous month's soil moisture at
                        t = 1 (in v/v), default value is 0.2
                      substep: logical 1 or 0; perform monthly substepping in
                        leaky bucket (1) or not (0)? Default value is 0.
        'intwindow': Integration window. Which months' growth responses should
                     be intregrated to compute the annual ring-width index?
                     Specified as a 2 x 1 vector of integer values. Both
                     elements are given in integer number of months since January
                     (July) 1st of the current year in the Northern (Southern)
                     hemisphere, and specify the beginning and end of the integration
                     window, respectively. Defaults is [1 ; 12] (eg. integrate
                     response to climate over the corresponding calendar year,
                     assuming location is in the northern hemisphere).
        'hydroclim': Value is a single character either taking value ['P'] or ['M'].
                     If ['M'], then 9th input is interpreted as an estimate of
                     soil moisture content (in v/v) rather than as precipitation.
                     Model default is to read in precipitation and use the CPC's
                     Leaky Bucket model of hydrology to estimate soil moisture,
                     however if soil moisture observations or a more sophisticated
                     estimate of moisture accounting for snow-related processes
                     is available, then using these data directly are recommended
                     (and will also speed up code).

        For more detailed documentation, see:
        1) Tolwinski-Ward et al., An efficient forward model of the climate
        controls on interannual variation in tree-ring width, Climate Dynamics (2011)
        DOI: 10.1007/s00382-010-0945-5

        2) Tolwinski-Ward et al., Erratum to: An efficient forward model of the climate
        controls on interannual variation in tree-ring width, Climate Dynamics (2011)
        DOI: 10.1007/s00382-011-1062-9

        3) Tolwinski-Ward et al., Bayesian parameter estimation and
        interpretation for an intermediate model of tree-ring width, Clim. Past
        (2013), DOI: 10.5194/cp-9-1-2013

        4) Documentation available with the model at http://www.ncdc.noaa.gov/paleo/softlib/softlib.html
    '''
    iyear = np.arange(syear, eyear+1)
    nyrs = np.size(iyear)
    T = T.reshape((nyrs, 12)).T
    P = P.reshape((nyrs, 12)).T

    # output vars
    Gr = np.ndarray((12, nyrs))
    gT = np.ndarray((12, nyrs))
    gM = np.ndarray((12, nyrs))
    M = np.ndarray((12, nyrs))
    potEv = np.ndarray((12, nyrs))

    if hydroclim == 'M':
        M = P
    elif substep == 1:
        M = leakybucket_submonthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0)['M']
    elif substep == 0:
        M = leakybucket_monthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0)['M']
    else:
        raise ValueError('substep must either be set to 1 or 0')

    gE = Compute_gE(phi)
    gT = (T-T1)/(T2-T1)
    gT[T<T1] = 0
    gT[T>T2] = 1
    gM = (M-M1)/(M2-M1)
    gM[M<M1] = 0
    gM[M>M2] = 1

    # Compute overall growth rate:
    if sensitivity == 'T':
        Gr = np.transpose(gT.T*gE)
    elif sensitivity == 'M':
        Gr = np.transpose(gM.T*gE)
    else:
        Gr = np.transpose(np.amin(np.array([gT, gM]), axis=0).T * gE)

    width = np.ones(nyrs)
    width[:] = np.nan
    if phi>0: # if site is in the Northern Hemisphere:
        if I_0<0: # if we include part of the previous year in each year's modeled growth:
            startmo = 13+I_0
            endmo = I_f
            # use average of growth data across modeled years to estimate first year's growth due
            # to previous year:
            width[0] = np.sum(Gr[0:endmo,0]) + np.sum(np.mean(Gr[startmo-1:12,:], axis=1))
            for cyear in range(1, nyrs):
                width[cyear] = np.sum(Gr[startmo-1:12,cyear-1]) + np.sum(Gr[0:endmo,cyear])
        else: # no inclusion of last year's growth conditions in estimates of this year's growth:
            startmo = I_0+1
            endmo = I_f
            for cyear in range(nyrs):
                width[cyear] = np.sum(Gr[startmo-1:endmo,cyear])
    elif phi<0: # if site is in the Southern Hemisphere:
        # (Note: in the Southern Hemisphere, ring widths are dated to the year in which growth began!)
        startmo = 7+I_0 # (eg. I_0 = -4 in SH corresponds to starting integration in March of cyear)
        endmo = I_f-6 # (eg. I_f = 12 in SH corresponds to ending integraion in June of next year)
        for cyear in range(nyrs-1):
            width[cyear] = np.sum(Gr[startmo-1:12,cyear]) + np.sum(Gr[0:endmo,cyear+1])
        # use average of growth data across modeled years to estimate last year's growth due to the next year:
        width[nyrs-1] = np.sum(Gr[startmo-1:12,nyrs-1])+ np.sum(np.mean(Gr[0:endmo,:], axis=1))

    trw = (width-np.mean(width))/np.std(width) # proxy series is standardized width.
    width_mean = np.mean(width)
    width_std = np.std(width)

    res_dict = {
        'trw': trw,
        'gT': gT,
        'gM': gM,
        'gE': gE,
        'Gr': Gr,
        'M': M,
        'potEv': potEv,
        'width': width,
        'width_mean': width_mean,
        'width_std': width_std,
    }

    return res_dict


def leakybucket_monthly(syear,eyear,phi,T,P,Mmax=0.76, Mmin=0.01, alph=0.093, m_th=4.886, mu_th=5.8, rootd=1000, M0=0.2):
    '''
    Translated from leackybucket_monthly.m - Simulate soil moisture with coarse monthly time step.

    Usage: [M,potEv,ndl,cdays] = leakybucket_monthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0)
       outputs simulated soil moisture and potential evapotranspiration.

    Inputs:
      syear = start year of simulation.
      eyear = end year of simulation.
      phi = latitude of site (in degrees N)
      T = (12 x Nyrs) matrix of ordered mean monthly temperatures (in degEes C)
      P = (12 x Nyrs) matrix of ordered accumulated monthly precipitation (in mm)
      Mmax = scalar maximum soil moisture held by the soil (in v/v)
      Mmin = scalar minimum soil moisture (for error-catching) (in v/v)
      alph = scalar runoff parameter 1 (in inverse months)
      m_th = scalar runoff parameter 3 (unitless)
      mu_th = scalar runoff parameter 2 (unitless)
      rootd = scalar root/"bucket" depth (in mm)
      M0 = initial value for previous month's soil moisture at t = 1 (in v/v)

    Outputs:
      M = soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs)
      potEv = potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm)
    '''
    iyear = np.arange(syear, eyear+1)
    nyrs = np.size(iyear)
    if len(np.shape(T)) == 1:
        T = T.reshape((nyrs, 12)).T
    if len(np.shape(P)) == 1:
        P = P.reshape((nyrs, 12)).T

    # output vars
    M = np.ndarray((12, nyrs))
    potEv = np.ndarray((12, nyrs))

    # Compute normalized daylength (neglecting small difference in calculation for leap-years)
    latr = phi*np.pi/180;  # change to radians
    ndays = np.array([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    cdays = np.cumsum(ndays)
    sd = np.arcsin(np.sin(np.pi*23.5/180) * np.sin(np.pi * ((np.arange(1,366) - 80)/180)))   # solar declination
    y = -np.tan(np.ones(365)* latr) * np.tan(sd)
    y[y>=1] = 1
    y[y<=-1] = -1
    hdl = np.arccos(y)
    dtsi = (hdl* np.sin(np.ones(365)*latr)*np.sin(sd))+(np.cos(np.ones(365)*latr)*np.cos(sd)*np.sin(hdl))
    ndl=dtsi/np.max(dtsi) # normalized day length

    # calculate mean monthly daylength (used for evapotranspiration in soil moisture calcs)
    jday = cdays[0:12] + .5*ndays[1:13]
    m_star = 1-np.tan(latr)*np.tan(23.439*np.pi/180*np.cos(jday*np.pi/182.625))
    mmm = np.empty(12)
    mmm[:] = np.nan

    for mo in range(12):
        if m_star[mo] < 0:
            mmm[mo] = 0
        elif m_star[mo] >0 and m_star[mo] < 2:
            mmm[mo] = m_star[mo]
        elif m_star[mo] > 2:
            mmm[mo] = 2


    nhrs = 24*np.arccos(1-mmm)/np.pi # the number of hours in the day in the middle of the month
    L = (ndays[1:13]/30)*(nhrs/12) # mean daylength in each month.


    for cyear in range(nyrs):    # begin cycling over years
        for t in range(12):      # begin cycling over months in a year
            # Compute potential evapotranspiration for current month after Thornthwaite:
            if T[t,cyear] < 0:
                Ep = 0
            elif T[t,cyear] >= 0 and T[t,cyear] < 26.5:
                istar = T[:,cyear]/5
                istar[istar<0] = 0

                I = np.sum(istar**1.514)
                a = (6.75e-7)*I**3 - (7.71e-5)*I**2 + (1.79e-2)*I + .49
                Ep = 16*L[t]*(10*T[t,cyear]/I)**a
            elif T[t,cyear] >= 26.5:
                Ep = -415.85 + 32.25*T[t,cyear] - .43* T[t,cyear]**2

            potEv[t,cyear] = Ep
            # Now calculate soil moisture according to the CPC Leaky Bucket model (see J. Huang et al, 1996).
            if t > 0:
                # evapotranspiration:
                Etrans = Ep*M[t-1,cyear]*rootd/(Mmax*rootd)
                # groundwater loss via percolation:
                G = mu_th*alph/(1+mu_th)*M[t-1,cyear]*rootd
                # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
                R = P[t,cyear]*(M[t-1,cyear]*rootd/(Mmax*rootd))**m_th + (alph/(1+mu_th))*M[t-1,cyear]*rootd
                dWdt = P[t,cyear] - Etrans - R - G
                M[t,cyear] = M[t-1,cyear] + dWdt/rootd
            elif t == 0 and cyear > 0:
                # evapotranspiration:
                Etrans = Ep*M[11,cyear-1]*rootd/(Mmax*rootd)
                # groundwater loss via percolation:
                G = mu_th*alph/(1+mu_th)*M[11,cyear-1]*rootd
                # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
                R = P[t,cyear]*(M[11,cyear-1]*rootd/(Mmax*rootd))**m_th + (alph/(1+mu_th))*M[11,cyear-1]*rootd
                dWdt = P[t,cyear] - Etrans - R - G
                M[t,cyear] = M[11,cyear-1] + dWdt/rootd
            elif t == 0 and cyear == 0:
                if M0 < 0:
                    M0 = .20
                # evapotranspiration (take initial soil moisture value to be 200 mm)
                Etrans = Ep*M0*rootd/(Mmax*rootd)
                # groundwater loss via percolation:
                G = mu_th*alph/(1+mu_th)*(M0*rootd)
                # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
                R = P[t,cyear]*(M0*rootd/(Mmax*rootd))**m_th + (alph/(1+mu_th))*M0*rootd
                dWdt = P[t,cyear] - Etrans - R - G
                M[t,cyear] = M0 + dWdt/rootd

            # error-catching:
            if M[t,cyear] <= Mmin:
                M[t,cyear] = Mmin
            if M[t,cyear] >= Mmax:
                M[t,cyear] = Mmax
            if np.isnan(M[t,cyear])==1:
                M[t,cyear] = Mmin

    res_dict = {
        'M': M,
        'potEv': potEv,
        'ndl': ndl,
        'cdays': cdays,
    }
    return res_dict



def leakybucket_submonthly(syear,eyear,phi,T,P,Mmax=0.76, Mmin=0.01, alph=0.093, m_th=4.886, mu_th=5.8, rootd=1000, M0=0.2):
    '''
    NOTE: This function is not tested yet!

    Translated from leackybucket_submonthly.m - Simulate soil moisture; substeps within monthly timesteps
    to better capture nonlinearities and improve moisture estimates.

    Usage: [M,potEv,ndl,cdays] = leakybucket_monthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0)
       outputs simulated soil moisture and potential evapotranspiration.

    Inputs:
      syear = start year of simulation.
      eyear = end year of simulation.
      phi = latitude of site (in degrees N)
      T = (12 x Nyrs) matrix of ordered mean monthly temperatures (in degEes C)
      P = (12 x Nyrs) matrix of ordered accumulated monthly precipitation (in mm)
      Mmax = scalar maximum soil moisture held by the soil (in v/v)
      Mmin = scalar minimum soil moisture (for error-catching) (in v/v)
      alph = scalar runoff parameter 1 (in inverse months)
      m_th = scalar runoff parameter 3 (unitless)
      mu_th = scalar runoff parameter 2 (unitless)
      rootd = scalar root/"bucket" depth (in mm)
      M0 = initial value for previous month's soil moisture at t = 1 (in v/v)

    Outputs:
      M = soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs)
      potEv = potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm)
    '''
    iyear = np.arange(syear, eyear+1)
    nyrs = np.size(iyear)
    if len(np.shape(T)) == 1:
        T = T.reshape((nyrs, 12)).T
    if len(np.shape(P)) == 1:
        P = P.reshape((nyrs, 12)).T

    # output vars
    M = np.ndarray((12, nyrs))
    potEv = np.ndarray((12, nyrs))

    # Compute normalized daylength (neglecting small difference in calculation for leap-years)
    latr = phi*np.pi/180;  # change to radians
    ndays = np.array([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    cdays = np.cumsum(ndays)
    sd = np.arcsin(np.sin(np.pi*23.5/180) * np.sin(np.pi * ((np.arange(1,366) - 80)/180)))   # solar declination
    y = -np.tan(np.ones(365)* latr) * np.tan(sd)
    y[y>=1] = 1
    y[y<=-1] = -1
    hdl = np.arccos(y)
    dtsi = (hdl* np.sin(np.ones(365)*latr)*np.sin(sd))+(np.cos(np.ones(365)*latr)*np.cos(sd)*np.sin(hdl))
    ndl=dtsi/np.max(dtsi) # normalized day length

    # calculate mean monthly daylength (used for evapotranspiration in soil moisture calcs)
    jday = cdays[0:12] + .5*ndays[1:13]
    m_star = 1-np.tan(latr)*np.tan(23.439*np.pi/180*np.cos(jday*np.pi/182.625))
    mmm = np.empty(12)
    mmm[:] = np.nan

    for mo in range(12):
        if m_star[mo] < 0:
            mmm[mo] = 0
        elif m_star[mo] >0 and m_star[mo] < 2:
            mmm[mo] = m_star[mo]
        elif m_star[mo] > 2:
            mmm[mo] = 2


    nhrs = 24*np.arccos(1-mmm)/np.pi # the number of hours in the day in the middle of the month
    L = (ndays[1:13]/30)*(nhrs/12) # mean daylength in each month.


    for cyear in range(nyrs):    # begin cycling over years
        for t in range(12):      # begin cycling over months in a year
            # Compute potential evapotranspiration for current month after Thornthwaite:
            if T[t,cyear] < 0:
                Ep = 0
            elif T[t,cyear] >= 0 and T[t,cyear] < 26.5:
                istar = T[:,cyear]/5
                istar[istar<0] = 0

                I = np.sum(istar**1.514)
                a = (6.75e-7)*I**3 - (7.71e-5)*I**2 + (1.79e-2)*I + .49
                Ep = 16*L[t]*(10*T[t,cyear]/I)**a
            elif T[t,cyear] >= 26.5:
                Ep = -415.85 + 32.25*T[t,cyear] - .43* T[t,cyear]**2

            potEv[t,cyear] = Ep
            # Now calculate soil moisture according to the CPC Leaky Bucket model (see J. Huang et al, 1996).
            # (see J. Huang et al, 1996). Set n-steps according to 2 mm increments
            # have to update alpha and Ep as well - 2 mm increments came from
            # testing by K. Georgakakos, but one could use 5 or more, with less "accurate" results.
            # Stepping is necessary because the parametization is linearized around init condition.

            dp = 2.0 # mm of precip per increment
            nstep = np.floor(P[t,cyear]/dp)+1 # number of sub-monthly substeps
            Pinc = P[t,cyear]/nstep # precip per substep
            alphinc = alph/nstep  # runoff rate per substep time interval
            Epinc = Ep/nstep # potential evapotrans per substep.

            # handling for sm_init
            if t > 0:
                M0=M[t-1,cyear]
            elif (t == 0) and (cyear > 0):
                M0=M[11,cyear-1]
            sm0 = M0

            for istep in range(nstep):
                # evapotranspiration:
                Etrans = Epinc*sm0*rootd/(Mmax*rootd)
                # groundwater loss via percolation:
                G = mu_th*alphinc/(1+mu_th)*sm0*rootd
                # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
                R = Pinc*(sm0*rootd/(Mmax*rootd))**m_th + (alphinc/(1+mu_th))*sm0*rootd
                dWdt = Pinc - Etrans - R - G
                sm1 = sm0 + dWdt/rootd
                sm0=np.max([sm1,Mmin])
                sm0=np.min([sm0,Mmax])

            M[t,cyear] = sm0

            # error-catching:
            if M[t,cyear] <= Mmin:
                M[t,cyear] = Mmin
            if M[t,cyear] >= Mmax:
                M[t,cyear] = Mmax
            if np.isnan(M[t,cyear])==1:
                M[t,cyear] = Mmin

    res_dict = {
        'M': M,
        'potEv': potEv,
        'ndl': ndl,
        'cdays': cdays,
    }
    return res_dict


def Compute_gE(phi):
    gE = np.ones(12)
    gE[:] = np.nan

    # Compute normalized daylength (neglecting small difference in calculation for leap-years)
    latr = phi*np.pi/180;  # change to radians
    ndays = np.array([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    cdays = np.cumsum(ndays)
    sd = np.arcsin(np.sin(np.pi*23.5/180) * np.sin(np.pi * ((np.arange(1,366) - 80)/180)))   # solar declination
    y = -np.tan(np.ones(365)* latr) * np.tan(sd)
    y[y>=1] = 1
    y[y<=-1] = -1
    hdl = np.arccos(y)
    dtsi = (hdl* np.sin(np.ones(365)*latr)*np.sin(sd))+(np.cos(np.ones(365)*latr)*np.cos(sd)*np.sin(hdl))
    ndl=dtsi/np.max(dtsi) # normalized day length

    # calculate mean monthly daylength (used for evapotranspiration in soil moisture calcs)
    jday = cdays[0:12] + .5*ndays[1:13]
    m_star = 1-np.tan(latr)*np.tan(23.439*np.pi/180*np.cos(jday*np.pi/182.625))
    mmm = np.empty(12)
    mmm[:] = np.nan

    for mo in range(12):
        if m_star[mo] < 0:
            mmm[mo] = 0
        elif m_star[mo] >0 and m_star[mo] < 2:
            mmm[mo] = m_star[mo]
        elif m_star[mo] > 2:
            mmm[mo] = 2

    for t in range(12):
        gE[t] = np.mean(ndl[cdays[t]:cdays[t+1]])

    return gE