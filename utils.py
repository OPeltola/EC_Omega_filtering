import pandas as pd
import numpy as np
import datetime as dt
from cmcrameri import cm as cmp
import pysolar

import constants as c


def NEON_to_ameriflux(neon_id):

    ID_conversion = {'BARR':'US-xBA'
                     ,'BART':'US-xBR'
                     ,'CPER':'US-xCP'
                     ,'DCFS':'US-xDC'
                     ,'DELA':'US-xDL'
                     ,'DSNY':'US-xDS'
                     ,'GUAN':'PR-xGU'
                     ,'HARV':'US-xHA'
                     ,'HEAL':'US-xHE'
                     ,'JERC':'US-xJE'
                     ,'JORN':'US-xJR'
                     ,'KONA':'US-xKA'
                     ,'KONZ':'US-xKZ'
                     ,'LAJA':'PR-xLA'
                     ,'MOAB':'US-xMB'
                     ,'NIWO':'US-xNW'
                     ,'NOGP':'US-xNG'
                     ,'OAES':'US-xAE'
                     ,'ONAQ':'US-xNQ'
                     ,'ORNL':'US-xRN'
                     ,'OSBS':'US-xSB'
                     ,'SCBI':'US-xSC'
                     ,'SERC':'US-xSE'
                     ,'SRER':'US-xSR'
                     ,'STEI':'US-xST'
                     ,'STER':'US-xSL'
                     ,'TOOL':'US-xTL'
                     ,'UNDE':'US-xUN'
                     ,'WOOD':'US-xWD'}

    id_out = neon_id
    if neon_id in ID_conversion.keys():
        id_out = ID_conversion[neon_id]
    else:
        print('Error in converting NEON ID to Ameriflux ID')

    return id_out

def piecewise_linear(x, x1, y0, y1):
    return np.piecewise(x, [x < x1], [lambda x:(y1-y0)/x1*x + y0, lambda x:y1])


def calculate_MCC(obs,pred):
    indx = (~obs.isna()) & (~pred.isna())
    obs = obs[indx]
    pred = pred[indx]

    tp = ((obs==1) & (pred==1)).sum()
    fp = ((obs==0) & (pred==1)).sum()
    fn = ((obs==1) & (pred==0)).sum()
    tn = ((obs==0) & (pred==0)).sum()


    x = (tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)
    if x==0:x=1 # see https://en.wikipedia.org/wiki/Phi_coefficient
    return ((tp * tn) - (fp * fn))/x**0.5


def quality_filter(df,site):
    dfout = df.copy()
    flxvars = ['H','TAU','LE','FC','SH','SLE','SC','NEE','USTAR','W_SIGMA','ZL','MO_LENGTH','WS','WD']

    dfout.loc[df['W_SIGMA']/df['USTAR']>2,flxvars] = np.nan
    dfout.loc[df['W_SIGMA']/df['USTAR']<0.7,flxvars] = np.nan
    dfout.loc[df['USTAR']>2,flxvars] = np.nan
    dfout.loc[df['USTAR'].isna(),flxvars] = np.nan
    dfout.loc[df['W_SIGMA']>2,flxvars] = np.nan
    dfout.loc[df['WS']>20,flxvars] = np.nan
    dfout.loc[df['FC']>20,flxvars] = np.nan
    dfout.loc[df['FC']<-40,flxvars] = np.nan
    dfout.loc[(df['FC']<-10) & (df['SW_IN']<10),flxvars] = np.nan


    dfout.loc[df['NEE']>20,'NEE'] = np.nan
    dfout.loc[df['NEE']<-40,'NEE'] = np.nan
    dfout.loc[(df['NEE']<-10) & (df['SW_IN']<10),'NEE'] = np.nan


    ndflag = pd.Series(index=df.index,data=0)
    ndflag[df['SW_IN']>10] = 1
    dfout['FC'] = despike_Papale_2006(dfout['FC'],ndflag,lim=4)
    dfout['NEE'] = despike_Papale_2006(dfout['NEE'],ndflag,lim=4)

    if site=='SE-Deg':
        dfout.loc[(dfout.index>dt.datetime(2020,2,1)) & (dfout.index<dt.datetime(2020,5,13)),'TA_5'] = np.nan

    return dfout

def despike_Papale_2006(x,ndflag,lim=None):


    if lim is None:
        lim = 7

    dat = x.copy()

    di = dat.diff()
    di.iloc[0] = 0
    di2 = dat.diff()
    di2[di2.index[0:-1]] = di2.iloc[1:].values
    di2.iloc[-1] = 0

    di3 = di - di2

    nd = pd.Series(index=di3.index)
    nd[ndflag==0] = di3[ndflag==0]
    dd = pd.Series(index=di3.index)
    dd[ndflag==1] = di3[ndflag==1]

    nMd = nd.rolling(window=dt.timedelta(days=13)).median()
    dMd = dd.rolling(window=dt.timedelta(days=13)).median()

    nMAD = (nd-nMd).abs().rolling(window=dt.timedelta(days=13)).median()
    dMAD = (dd-dMd).abs().rolling(window=dt.timedelta(days=13)).median()

    spike = pd.Series(index=dat.index,data=False)
    spike[(nd<nMd-lim*nMAD/0.6745) | (nd>nMd+lim*nMAD/0.6745)] = True
    spike[(dd<dMd-lim*dMAD/0.6745) | (dd>dMd+lim*dMAD/0.6745)] = True

    dat[spike] = np.nan


    return dat


def calculate_reference_flux(fluxin,SWIN,USTARin,SWINlm=None,window=None,min_periods=None):
    """Calculating night time reference flux.

    Args:
        fluxin (Series): gas flux time series.
        SWIN (Series): Incoming short wave radiation in W m-2.
        USTARin (Series): Friction velocity in m s-1

    Returns:
        _type_: _description_
    """    
    if SWINlm is None:SWINlm = 10
    if window is None:window = dt.timedelta(days=10)
    if min_periods is None:min_periods = 5
    USTlm = .75

    flux = fluxin.copy()
    UST = USTARin.copy()

    flux[(SWIN>SWINlm) | (SWIN.isna())] = np.nan
    UST[(SWIN>SWINlm) | (SWIN.isna())] = np.nan
    flux[(UST<UST.quantile(USTlm)) | (UST.isna())] = np.nan
    flux[flux<0] = np.nan

    # running median within window
    refflux = flux.rolling(window = window, center = True, min_periods = min_periods).quantile(.5)
    refflux = refflux.interpolate(method='slinear',limit=int(48*5),limit_direction='both')
    return refflux



def fix_temperature_bias(TPOT,H,USTAR,z):

    TPOT_fixed = pd.DataFrame(index=TPOT.index,columns=TPOT.columns)
    bias = pd.Series(index=TPOT.columns)

    if len(TPOT.columns)>2:
        NN = (H.abs()<5) & (USTAR>0.5)
        if NN.sum()<100:
            NN = (H.abs()<10) & (USTAR>0.3)

        for col in TPOT.columns:
            TPOT_tmp = TPOT.copy()
            TPOT_tmp[col] = np.nan
            TPOT_tmp = TPOT_tmp.rename(columns=z).T
            TPOT_tmp = TPOT_tmp.sort_index()
            TPOT_tmp = TPOT_tmp.interpolate(method="linear", order=1, limit_direction="both")
            TPOT_tmp = TPOT_tmp.T
            TPOT_tmp = TPOT_tmp.rename(columns={v: k for k, v in z.items()})
            bias[col] = (TPOT_tmp.loc[NN,col]-TPOT.loc[NN,col]).median()
            TPOT_fixed[col] = TPOT[col] + bias[col]
    else:
        TPOT_fixed = TPOT
        
    return TPOT_fixed,bias


    # bias = TPOT[(H.abs()<5) & (USTAR>0.5)].median()

def convert_to_potential_temperature(Tmat,z):
    Tpot = pd.DataFrame(index=Tmat.index)
    for col in Tmat.columns:
        Tpot[col] = Tmat[col]+z[col]*c.g/c.cp
    return Tpot


def calculate_omega(W_SIGMA,Tmat,Uh,PAI,canz,z,zT,temporal_interpolation=None,verbose=None):
    """Function for calculating decoupling parameter Omega. 
    All the Series and DataFrames need to have the same indices.

    Args:
        W_SIGMA (Series): time series of vertical wind component std, unit m s-1
        Tmat (DataFrame): time series of air temperature profile, unit K
        Uh (Series): time series of wind speed at canopy height, unit m s-1
        PAI (Series): time series of plant area index, unit m2 m-2
        canz (float): canopy height in meters, unit m
        z (float): Eddy covariance measurement height, unit m
        zT (dict): measurement heights (unit m) for air temperatures. Keys much match Tmat columns
        temporal_interpolation (boolean): Flag for interpolating temperature data over time domain. Defaults to True
        verbose (boolean): Flag for printing output on the screen. Defaults to False

    Returns:
        Series: Time series of dimensionless decoupling parameter Omega
        Series: Time series of downward speed needed for the air parcel to reach the ground, unit m s-1
        Series: Time series of temperature difference between the air parcel and air column below EC, unit K
    """    
    NT = 273.15
    cd = 0.2

    if temporal_interpolation is None:
        temporal_interpolation = True
    if verbose is None:
        verbose = False

    for col in Tmat.columns:
        if Tmat[col].median()<100:
            print(col + ' appears to be in Celsius. Converting to Kelvin.')
            Tmat[col] = Tmat[col]+NT

    Tmat = Tmat.rename(columns=zT)
    # at least two temperature measurement heights are needed for Omega calculation
    indx = (~Tmat.isna()).sum(axis=1)<2

    Tmati = Tmat.copy()

    # making sure that there is a value at EC height
    if z not in Tmati.columns:
        Tmati[z] = np.nan
    # making sure that there is a value at ground level
    if 0 not in Tmati.columns:
        Tmati[0] = np.nan
    # making sure that there is a value at canopy height
    if (canz>1) & (canz not in Tmati.columns):
        Tmati[canz] = np.nan

    # spatial interpolation
    Tmati[~indx] = Tmati[~indx].interpolate(method='polynomial',order=1,axis=1,fill_value='extrapolate')

    Tmati = Tmati.reindex(sorted(Tmati.columns), axis=1)
    
    # temporal interpolation, at max 3 h 
    Tmati.loc[~indx,:] = Tmati.loc[~indx,:].interpolate(method='slinear',axis=0,limit_direction='both',limit_area='inside',limit=int(3*60/30))
    Tmati = Tmati.interpolate(method='index',axis=1,limit_direction='both',limit_area='inside')
    
    # averaging over duplicate measurements
    Tmati = Tmati.T.groupby(Tmati.T.index).mean().T


    if canz>1:
        alpha = 2.0
        beta = 1.5
        gamma = -(np.exp(-beta-alpha)-1)/(beta+alpha)
        # drag force term
        drg = gamma*cd*PAI*Uh

        # temperature of the air parcel
        T_ap = Tmati.loc[:,z]

        # mean below canopy height
        cols = Tmati.columns[(Tmati.columns<=canz)]
        reftheta = pd.Series(np.trapz(Tmati.loc[:,cols].T, x=cols, axis=0)/canz,index = W_SIGMA.index)


        reftheta1 = reftheta.copy()
        # not allowing negative buoyancy force, conservative approach
        # reftheta[reftheta>T_ap] = T_ap

        # wec at canopy height
        buoyancy = 2*c.g*canz*(T_ap-reftheta)/reftheta
        wec1 = -drg-(drg**2+buoyancy)**0.5
        
        # mean between canopy height and EC height
        cols = Tmati.columns[(Tmati.columns>=canz) & (Tmati.columns<=z)]
        reftheta = pd.Series(np.trapz(Tmati.loc[:,cols].T, x=cols, axis=0)/(z-canz),index = W_SIGMA.index)
        # not allowing negative buoyancy force, conservative approach
        # reftheta[reftheta>T_ap] = T_ap

        wec = -(wec1**2+2*c.g*(z-canz)*(T_ap-reftheta)/reftheta)**0.5

        wec[drg**2+buoyancy<0] = 0
    else:
        if verbose:
            print('Skipping canopy drag due to low canopy')
        # temperature of the air parcel
        T_ap = Tmati.loc[:,z]

        # mean below EC height
        cols = Tmati.columns[Tmati.columns<=z]
        reftheta = pd.Series(np.trapz(Tmati.loc[:,cols].T, x=cols, axis=0)/z,index = W_SIGMA.index)
        
        buoyancy = 2*c.g*z*(T_ap-reftheta)/reftheta
        wec = -(buoyancy)**0.5
        wec[buoyancy<0] = 0

        reftheta1 = reftheta.copy()


    omega = W_SIGMA/wec.abs()

    # at least two temperature measurements needed
    omega[indx] = np.nan
    wec[indx] = np.nan

    return omega,wec,T_ap-reftheta1


def calculate_Uh(WS,ZL,USTAR,z,canz,d):

    # stability correction function
    phi = 1+5*ZL
    phi[ZL<0] = 1
    
    dU = USTAR/(c.k*(z-d))*phi*(canz-z)

    # dU cannot be positive
    dU[dU>0] = 0
    dU[dU.isna()] = dU.mean()

    Uh = WS+dU
    Uh[Uh<0.5] = 0.5

    Uh[canz==0] = 0

    return Uh

def define_colors(z):

    val = np.random.rand(1,1)[0][0]*0.25
    if z<10:
        clr = list(cmp.batlow(0.05))
    elif (z>=10) & (z<30):
        clr = list(cmp.batlow(0.5))
    elif (z>=30):
        clr = list(cmp.batlow(0.8))

    clr = list(clr[0:3]+val)
    clr = [i if i<1 else 1.0 for i in clr]
    clr = [i if i>0 else 0.0 for i in clr]
    clr.append(1.0)
    return tuple(clr)

def interpolate_LAI(LAI,tvec):

    if len(LAI.index)>len(LAI.index.unique()):
        LAI = LAI.groupby(LAI.index).mean()

    if (~LAI.isna()).sum()>5:

        LAIi = pd.Series(index=tvec.union(LAI.index))
        LAIi = LAIi.sort_index()
        LAIi[LAI.index] = LAI
        LAIi = LAIi.resample('30min').mean()
        # gaps shorter than 3 months filled with linear interpolation
        LAIi = LAIi.interpolate(method='slinear',limit_area='inside',limit=int(48*30*3))
        # gaps longer than 3 months filled with mean annual pattern
        monthly_means = LAIi.groupby(LAIi.index.month).mean()
        for month in monthly_means.index:
            indx = (LAIi.index.month==month) & (LAIi.index.day==15) & (LAIi.index.hour==0) & (LAIi.index.minute==0) & (LAIi.isna())
            LAIi[indx] = monthly_means[month]
        LAIi = LAIi.interpolate(method='slinear',limit_area='inside')

        LAIi = LAIi.bfill()
        LAIi = LAIi.ffill()

        LAIout = LAIi[tvec]
    else:
        LAIout = pd.Series(index=tvec,data=LAI.median())
    
    return LAIout

    # plt.plot(LAI,'ko')
    # plt.plot(LAIi,'b-')
    # plt.show()

def calculate_Obukhov_length(H,UST,PA,T,H2O=None):
    """Calculation of Obukhov length

    Args:
        H (Series): Sensible heat flux in W m-2
        UST (Series): Friction velocity in m s-1
        PA (Series): Air pressure in kPa
        T (Series): Air temperature in C
        H2O (Series, optional): H2O dry mixing ratio in mmol mol-1. Defaults to None.

    Returns:
        Series: Obukhov length
    """

    if H2O is None:
        H2O = pd.Series(index=H.index,data=0)
    
    # air density
    rho = 1e3*PA*c.Md/(c.R*(T+c.NT))
    rho[rho.isna()] = rho.median()
    
    # virtual temperature
    Tv = (T+c.NT)*(1+0.61*H2O*1e-3)
    Tv[Tv.isna()] = T[Tv.isna()]+c.NT
    Tv[Tv.isna()] = Tv.median()

    TST = -H/(rho*c.cp)/UST
    Lmo = UST**2*Tv/(c.k*c.g*TST)

    return Lmo

def calculate_RH(H2O,T,PA):
    """Calculation of relative humidity (RH) from H2O mixing ratio, air temperature and atmospheric pressure

    Args:
        H2O (Series): H2O mixing ratio in mmol/mol
        T (Series): air temperature in Celsius degrees
        PA (Series): air pressure in kPa

    Returns:
        Series: relative humidity in %
    """    

    if np.median(T)>100:
        msg = 'Air temperature appears to be given in Kelvin. Converting to Celsius degrees.'                
        print(msg)
        T = T-c.NT
    if np.median(PA)>300:
        msg = 'Air pressure appears to be given in hPa. Converting to kPa.'                
        print(msg)
        PA = PA*1e-1
    
    esat = calculate_esat(T)

    # converting to mol/mol
    H2O = H2O*1e-3
    # H2O partial pressure
    ph2o = PA*H2O
    # relative humidity in %
    RH = ph2o/esat*100 
    return RH

def calculate_VPD(RH,T):
    """Calculation of vapor pressure deficit (VPD) from relative humidity and air temperature

    Args:
        RH (Series): relative humidity in %
        T (Series): air temperature in Celsius degrees

    Returns:
        Series: vapor pressure deficit in kPa
    """    

    if np.median(T)>100:
        msg = 'Air temperature appears to be given in Kelvin. Converting to Celsius degrees.'                
        print(msg)
        T = T-c.NT
    
    esat = calculate_esat(T)
    VPD = esat*(1-RH/100)
    return VPD


def calculate_esat(T):
    """Calculation of saturation vapor pressure (esat) based on August-Roche-Magnus formula

    Args:
        T (Series): Air temperature in Celsius degrees

    Returns:
        Series: saturation vapor pressure in kPa
    """    
    if np.median(T)>100:
        msg = 'Air temperature appears to be given in Kelvin. Converting to Celsius degrees.'                
        print(msg)
        T = T-c.NT

    # August-Roche-Magnus formula
    esat = 6.1094*np.exp(17.625*T/(T+243.04))
    # from hPa to kPa
    esat = esat*1e-1
    return esat



def NEON_site_name(siteID):

    namedict = {'BART':'BartlettExperimentalForest',
        'BLAN':'',
        'DSNY':'DisneyWildernessPreserve',
        'GUAN':'GuanicaForest',
        'GRSM':'',
        'DELA':'DeadLake',
        'DCFS':'',
        'CPER':'CentralPlainsExperimentalRange',
        'CLBJ':'',
        'ABBY':'',
        'BARR':'',
        'BONA':'',
        'DEJU':'',
        'HARV':'HarvardForest',
        'JERC':'JonesEcologicalResearchCenter',
        'JORN':'Jornada',
        'LAJA':'LajasExperimentalStation',
        'MOAB':'Moab',
        'NIWO':'NiwotRidgeMountainResearchStation',
        'STER':'NorthSterling',
        'ORNL':'OakRidge',
        'ONAQ':'OnaquiAult',
        'OSBS':'OrdwaySwisherBiologicalStation',
        'SRER':'SantaRita',
        'SCBI':'SmithsonianConservationBiologyInstitute',
        'SERC':'SmithsonianEnvironmentalResearchCenter',
        'STEI':'SteigerwaldtLandServices',
        'TALL':'TalladegaNationalForest',
        'UNDE':'Underc',
        'WOOD':'Woodworth'}
    
    if siteID in namedict.keys():
        name = namedict[siteID]
    else:
        name = ''

    return name


def MDS(df,vr,vrgf,drvvr,tol,Nmin=None,vrunc=None,timeisna=None,method=None,subsample=None):

    if Nmin is None:
        Nmin=2
    if Nmin is not round(Nmin):
        msg = 'MDS algorithm accepts only integer number for minimum amount of data. Rounding to integer.'
        print(msg)
        Nmin=round(Nmin)
    if Nmin<2:
        msg = 'Nmin must be larger than 1. Forcing to 2.'
        print(msg)
        Nmin=2
    if method is None:
        method = 'mean'
    if subsample is None:
        subsample = dict()
        for drv in drvvr:
            subsample[drv] = False

    # periods with gaps that are filled. If None, then filling all the gaps
    if timeisna is None:
        
        gap_flag = pd.Series(index=df.index,data=False)
        gap_flag[(df[vrgf].isnull())] = True
        N=48*60
        gap_flag = (gap_flag.groupby((gap_flag != gap_flag.shift()).cumsum()).transform(lambda x: len(x) < N)*gap_flag)
        timeisna = df.loc[gap_flag & ((~df[drvvr].isnull()).all(axis=1)),'time']
        print(str(len(timeisna)) +' (' + str(int(len(timeisna)/len(df)*100)) + ' %) gaps to be filled.')


    # MDS algorithm
    # step 1: Look-up table with window size +/- 7 days
    df = LUT(df,vr,vrgf,drvvr,tol,deltadays=7,Nmin=Nmin,vrunc=vrunc,timeisna=timeisna,method=method,subsample=subsample)
    intrsct =  list(set(timeisna) & set(df.loc[df[vrgf].isnull(),'time']))
    timeisna = pd.Series(index=intrsct,data=intrsct,name=timeisna.name,dtype='datetime64[ns]')
    # df[vrgf+'_1'] = df[vrgf]
    # step 2: Look-up table with window size +/- 14 days
    df = LUT(df,vr,vrgf,drvvr,tol,deltadays=14,Nmin=Nmin,vrunc=vrunc,timeisna=timeisna,method=method,subsample=subsample)
    intrsct =  list(set(timeisna) & set(df.loc[df[vrgf].isnull(),'time']))
    timeisna = pd.Series(index=intrsct,data=intrsct,name=timeisna.name,dtype='datetime64[ns]')
    # df[vrgf+'_2'] = df[vrgf]
    # step 3: Look-up table with main driver only, window size +/- 7 days
    df = LUT(df,vr,vrgf,[drvvr[0]],tol,deltadays=7,Nmin=Nmin,vrunc=vrunc,timeisna=timeisna,method=method,subsample=subsample)
    intrsct =  list(set(timeisna) & set(df.loc[df[vrgf].isnull(),'time']))
    timeisna = pd.Series(index=intrsct,data=intrsct,name=timeisna.name,dtype='datetime64[ns]')
    # df[vrgf+'_3'] = df[vrgf]
    # step 4: Mean diurnal course with window size of 0 (same day)
    df = MDC(df,vr,vrgf,deltadays=0,Nmin=Nmin,vrunc=vrunc,timeisna=timeisna,method=method)
    intrsct =  list(set(timeisna) & set(df.loc[df[vrgf].isnull(),'time']))
    timeisna = pd.Series(index=intrsct,data=intrsct,name=timeisna.name,dtype='datetime64[ns]')
    # df[vrgf+'_4'] = df[vrgf]
    # step 5: Mean diurnal course with window size of +/-1 day and +/-2 days
    df = MDC(df,vr,vrgf,deltadays=1,Nmin=Nmin,vrunc=vrunc,timeisna=timeisna,method=method)
    intrsct =  list(set(timeisna) & set(df.loc[df[vrgf].isnull(),'time']))
    timeisna = pd.Series(index=intrsct,data=intrsct,name=timeisna.name,dtype='datetime64[ns]')
    df = MDC(df,vr,vrgf,deltadays=2,Nmin=Nmin,vrunc=vrunc,timeisna=timeisna,method=method)
    intrsct =  list(set(timeisna) & set(df.loc[df[vrgf].isnull(),'time']))
    timeisna = pd.Series(index=intrsct,data=intrsct,name=timeisna.name,dtype='datetime64[ns]')
    # df[vrgf+'_5'] = df[vrgf]
    # step 6: Look-up table with window size of +/- 21 to +/- 70 days
    for tdays in range(21,70,7):
        df = LUT(df,vr,vrgf,drvvr,tol,deltadays=tdays,Nmin=Nmin,vrunc=vrunc,timeisna=timeisna,method=method,subsample=subsample)
        intrsct =  list(set(timeisna) & set(df.loc[df[vrgf].isnull(),'time']))
        timeisna = pd.Series(index=intrsct,data=intrsct,name=timeisna.name,dtype='datetime64[ns]')
    # df[vrgf+'_6'] = df[vrgf]
    # step 7: Look-up table with main driver only, window size of +/- 14 to +/- 70 days
    for tdays in range(14,70,7):
        df = LUT(df,vr,vrgf,[drvvr[0]],tol,deltadays=tdays,Nmin=Nmin,vrunc=vrunc,timeisna=timeisna,method=method,subsample=subsample)
        intrsct =  list(set(timeisna) & set(df.loc[df[vrgf].isnull(),'time']))
        timeisna = pd.Series(index=intrsct,data=intrsct,name=timeisna.name,dtype='datetime64[ns]')
    # df[vrgf+'_7'] = df[vrgf]
    # step 8: Mean diurnal course with window size of +/- 7 to +/- 210 days
    for tdays in range(7,210,7):
        df = MDC(df,vr,vrgf,deltadays=tdays,Nmin=Nmin,vrunc=vrunc,timeisna=timeisna,method=method)
        intrsct =  list(set(timeisna) & set(df.loc[df[vrgf].isnull(),'time']))
        timeisna = pd.Series(index=intrsct,data=intrsct,name=timeisna.name,dtype='datetime64[ns]')


    return df


def LUT(df,vr,vrgf,drvvr,tol,deltadays=None,Nmin=None,vrunc=None,timeisna=None,method=None,subsample=None):

    dfout = df.copy()

    if Nmin is None:
        Nmin=2
    if Nmin is not round(Nmin):
        msg = 'LUT algorithm accepts only integer number for minimum amount of data. Rounding to integer.'
        print(msg)
        Nmin=round(Nmin)
    if Nmin<2:
        msg = 'Nmin must be larger than 1. Forcing to 2.'
        print(msg)
        Nmin=2

    if deltadays is None:
        deltadays = 1
    if deltadays is not round(deltadays):
        msg = 'LUT algorithm accepts only integer number of days. Rounding to integer.'
        print(msg)
        deltadays = round(deltadays)
    if deltadays<1:
        msg = 'deltadays must be larger than 1. Forcing to 1.'
        print(msg)
        deltadays = 1

    if method is None:
        method = 'mean'

    if subsample is None:
        subsample = dict()
        for drv in drvvr:
            subsample[drv] = False

    # periods with gaps that are filled. If None, then filling all the gaps
    if timeisna is None:
        timeisna = df.loc[(df[vrgf].isnull()) & ((~df[drvvr].isnull()).all(axis=1)),'time']

    deltat = dt.timedelta(days=int(deltadays))

    df['bool'] = 1
    df.loc[df[vr].isnull(),'bool'] = 0
    df['group'] = ''

    # dftmp2 = df[vr]
    # asd = dict()
    # asd['time'] = pd.DataFrame(index=timeisna,columns=df.index)
    # asd['time'] = asd['time'].apply(df.index, axis=1)
    # asd['SW_IN'] = pd.DataFrame(index=timeisna,columns=df.index,data=df['SW_IN'])


    # looping over the gaps
    for t in timeisna:
        dftmp = df.loc[df['time'].between(t-deltat,t+deltat),drvvr+[vr,'time','group','bool']].copy()

        grp = 'a'
        for drv in drvvr:

            # value for the driver at this particular time
            drvval = dftmp.loc[dftmp['time']==t,drv].iloc[0]
            
            # tolerance around drvval
            tolval = tol[drv][0][1]
            for i in range(len(tol[drv])):
                if (drvval>tol[drv][i][0][0]) & (drvval<=tol[drv][i][0][1]):
                    tolval = tol[drv][i][1]

            grp = chr(ord(grp)+1)
            if subsample[drv]:
                lms = [drvval-tolval,drvval]
                dftmp.loc[dftmp[drv].between(lms[0],lms[1]),'group'] = dftmp.loc[dftmp[drv].between(lms[0],lms[1]),'group']+grp

                lms = [drvval,drvval+tolval]
                grp1 = chr(ord(grp)+1)
                dftmp.loc[(dftmp[drv].between(lms[0],lms[1])) & (~dftmp['group'].str.contains(grp)),'group'] = \
                    dftmp.loc[dftmp[drv].between(lms[0],lms[1]) & (~dftmp['group'].str.contains(grp)),'group']+grp1
                grp = grp1
            
            else:
                lms = [drvval-tolval,drvval+tolval]
                dftmp.loc[dftmp[drv].between(lms[0],lms[1]),'group'] = dftmp.loc[dftmp[drv].between(lms[0],lms[1]),'group']+grp
        dftmp = dftmp.loc[dftmp['group'].str.len()>=len(drvvr)]

        # amount of non-NaN values
        ndat = dftmp['bool'].sum()

        if ndat>=Nmin:
            # enough data for averaging            
            if method=='mean':
                dfout.loc[dfout['time']==t,vrgf] = dftmp.groupby('group').mean(numeric_only=True).mean(numeric_only=True)[vr]
            if method=='median':
                dfout.loc[dfout['time']==t,vrgf] = dftmp.groupby('group').median(numeric_only=True).median(numeric_only=True)[vr]
            if vrunc is not None:
                dfout.loc[dfout['time']==t,vrunc] = dftmp.groupby('group').std().mean(numeric_only=True)[vr]

        del dftmp

    return dfout


def MDC(df,vr,vrgf,deltadays=None,method=None,Nmin=None,vrunc=None,timeisna=None,deltat=None):

    dfout = df.copy()

    z4 = diurnal_pattern(df[vr],df['time'],method=method,deltadays=deltadays,Nmin=Nmin,deltat=deltat)
    if vrunc is not None:
        z5 = diurnal_pattern(df[vr],df['time'],method='std',deltadays=deltadays,Nmin=Nmin,deltat=deltat)



    # periods with gaps that are filled. If None, then filling all the gaps
    if timeisna is None:
        timeisna = df.loc[df[vrgf].isnull(),'time']
    
    for t in timeisna:
        dfout.loc[t,vrgf] = z4[t]
        if vrunc is not None:
            dfout.loc[t,vrunc] = z5[t]

    return dfout



def diurnal_pattern(dat,time,method=None,deltadays=None,Nmin=None,deltat=None):
    
    if Nmin is None:
        Nmin=1
    if Nmin is not round(Nmin):
        msg = 'Integer number for minimum amount of data should be given. Rounding to integer.'
        print(msg)
        Nmin=round(Nmin)
    if Nmin<1:
        msg = 'Nmin must be at least 1. Forcing to 1.'
        print(msg)
        Nmin=1

    if deltadays is None:
        deltadays = 1
    if deltadays is not round(deltadays):
        msg = 'Integer number of days should be given. Rounding to integer.'
        print(msg)
        deltadays = round(deltadays)
    if deltadays<1:
        msg = 'deltadays must be at least 1. Forcing to 1.'
        print(msg)
        deltadays = 1
    if Nmin>deltadays:
        msg = 'Nmin must be at least deltadays. Forcing to deltadays.'
        print(msg)
        Nmin = deltadays

    if method is None:
        method = 'mean'
    if (method!='mean') and (method!='median') and (method!='std'):
        msg = 'Unknown method given (%s) for calculating diurnal patterns. Using ''mean'' instead.'%(method,)
        print(msg)
        method = 'mean'
    if deltat is None:
        deltat = '30min'


    if (deltat!='30min') & (deltat!='1min'):
        msg = 'Unknown time step for time series given (%s). Only ''30min'' and ''1min'' are permitted. Skipping the calculation.'%(deltat,)
        print(msg)
        datout = pd.Series(index=dat.index)
    else:
        tstart = dt.datetime(time.dt.year.iloc[0],1,1).date()
        t = pd.date_range(start=tstart,end=time.dt.date.iloc[-1]+dt.timedelta(days=1), freq=deltat)
        t = t[0:-1]
        dat = dat.reindex(t)
        t = pd.Series(t,index=t)

        date = (t-t.iloc[0])/dt.timedelta(days=1)
        time2 = ((date-np.floor(date))*24)
        time2 = time2.round(decimals=5)
        time2 = time2.drop_duplicates().to_numpy()
        date = np.floor(date).unique()

        z3 = pd.DataFrame(data=dat.to_numpy().reshape(len(date),len(time2)))
        
        if method=='mean':
            z4 = pd.Series(data=z3.rolling(deltadays,center=True,min_periods=Nmin).mean().to_numpy().reshape(len(t),).T,index=t)
        elif method=='median':
            z4 = pd.Series(data=z3.rolling(deltadays,center=True,min_periods=Nmin).median().to_numpy().reshape(len(t),).T,index=t)
        elif method=='std':
            z4 = pd.Series(data=z3.rolling(deltadays,center=True,min_periods=Nmin).std().to_numpy().reshape(len(t),).T,index=t)


        datout = z4.reindex(index=time.index)

    return datout


def gapfill(df):
    

    dftmp = pd.DataFrame(index=df.index,columns=['time','FC','FC_GF_OMEGA','FC_GF_USTAR','NEE','NEE_GF_OMEGA','NEE_GF_USTAR','TA','SW_IN','VPD'])
    dftmp['time'] = dftmp.index
    dftmp['FC'] = df['FC'].copy()
    dftmp['FC_GF_OMEGA'] = df['FC_GF_OMEGA'].copy()
    dftmp['FC_GF_USTAR'] = df['FC_GF_USTAR'].copy()
    dftmp['NEE'] = df['NEE'].copy()
    dftmp['NEE_GF_OMEGA'] = df['NEE_GF_OMEGA'].copy()
    dftmp['NEE_GF_USTAR'] = df['NEE_GF_USTAR'].copy()

    dftmp['TA'] = df[[i for i in list(df.columns) if 'TA_' in i[0:3]]].mean(axis=1)
    if 'VPD' not in df.columns:
        if 'RH' not in df.columns:
            dftmp['RH'] = df[[i for i in list(df.columns) if 'RH_' in i[0:3]]].mean(axis=1)
        else:
            dftmp['RH'] = df['RH'].copy()
        dftmp['VPD'] = calculate_VPD(dftmp['RH'],dftmp['TA'])
    else:
        dftmp['VPD'] = df['VPD'].copy()
    if 'PPFD_IN' in df.columns:
        dftmp['SW_IN'] = df['PPFD_IN'].copy()*0.5
        if 'SW_IN' in df.columns:
            dftmp.loc[dftmp['SW_IN'].isna(),'SW_IN'] = df.loc[dftmp['SW_IN'].isna(),'SW_IN']
    elif 'SW_IN' in df.columns:
        dftmp['SW_IN'] = df['SW_IN'].copy()

    dftmp = InterpolateShortGaps(dftmp,'VPD',10)
    dftmp = InterpolateShortGaps(dftmp,'TA',10)
    dftmp = InterpolateShortGaps(dftmp,'SW_IN',10)
    dftmp['VPD_orig'] = dftmp['VPD'].copy()
    dftmp['TA_orig'] = dftmp['TA'].copy()
    dftmp['SW_IN_orig'] = dftmp['SW_IN'].copy()
    for deltaday in [5,10,30]:
        dftmp = MDC(dftmp,'VPD_orig','VPD',deltadays=deltaday)
        dftmp = MDC(dftmp,'TA_orig','TA',deltadays=deltaday)
        dftmp = MDC(dftmp,'SW_IN_orig','SW_IN',deltadays=deltaday)

    tol = dict()
    tol['SW_IN'] = [[[-np.inf,50],20],[[50,np.inf],50]]
    tol['TA'] = [[[-np.inf,np.inf],2.5]]
    tol['VPD'] = [[[-np.inf,np.inf],0.5]]
    if (pd.isnull(dftmp['FC']).sum()<len(dftmp)) & (pd.isnull(dftmp['SW_IN']).sum()<len(dftmp)) & (pd.isnull(dftmp['TA']).sum()<len(dftmp)) & (pd.isnull(dftmp['VPD']).sum()<len(dftmp)):
        dftmp = MDS(dftmp,'FC','FC_GF_OMEGA',['SW_IN','TA','VPD'],tol)
        dftmp = MDS(dftmp,'FC','FC_GF_USTAR',['SW_IN','TA','VPD'],tol)
    else:
        dftmp['FC_GF_OMEGA'] = np.nan
        dftmp['FC_GF_USTAR'] = np.nan
        print(['All FC data missing. No gapfilling'])
        
    if (pd.isnull(dftmp['NEE']).sum()<len(dftmp)) & (pd.isnull(dftmp['SW_IN']).sum()<len(dftmp)) & (pd.isnull(dftmp['TA']).sum()<len(dftmp)) & (pd.isnull(dftmp['VPD']).sum()<len(dftmp)):
        dftmp = MDS(dftmp,'NEE','NEE_GF_OMEGA',['SW_IN','TA','VPD'],tol)
        dftmp = MDS(dftmp,'NEE','NEE_GF_USTAR',['SW_IN','TA','VPD'],tol)
    else:
        dftmp['NEE_GF_OMEGA'] = np.nan
        dftmp['NEE_GF_USTAR'] = np.nan
        print(['All NEE data missing. No gapfilling'])

    df['FC_GF_OMEGA'] = dftmp['FC_GF_OMEGA']
    df['FC_GF_USTAR'] = dftmp['FC_GF_USTAR']
    df['NEE_GF_OMEGA'] = dftmp['NEE_GF_OMEGA']
    df['NEE_GF_USTAR'] = dftmp['NEE_GF_USTAR']

    del dftmp

    return df
    

def InterpolateShortGaps(df,vr,N,method=None):
    """Filling short gaps within the time series with interpolation.

    Args:
        df (DataFrame): data
        vr (str): column in df to be filled
        N (int): max length for the gaps to be filled
        method (str, optional): method used for interpolation. Defaults to None.

    Returns:
        DataFrame: dataframe where the short gaps in column are filled with selected interpolation method
    """


    if method is None:
        method = 'linear'

    if (method!='linear') & (method!='nearest') & (method!='quadratic') & (method!='cubic') & (method!='spline'):
        msg = 'Unknown method given (%s) for interpolating short gaps. Using ''linear'' instead.'%(method,)
        print(msg)
        method = 'linear'
    
    dat = df[vr].copy()

    dat_interpolated = dat.interpolate(method=method)

    mask = dat.isna()
    x = (mask.groupby((mask != mask.shift()).cumsum()).transform(lambda x: len(x) > N)*mask)
    dat_interpolated = dat_interpolated.loc[~x]

    df[vr] = dat_interpolated

    return df


def reorganise_for_fingerprint(dat,tstart,deltat=None):
    
    if deltat is None:
        deltat = dt.timedelta(minutes=30)
    
    time = dat.index

    t = pd.date_range(start=tstart,end=time.date[-1]+dt.timedelta(days=1), freq=deltat)
    t = t[0:-1]
    z2 = dat.reindex(t)
    t = pd.Series(t,index=t)

    date = (t-t.iloc[0])/dt.timedelta(days=1)
    time2 = ((date-np.floor(date))*24)
    time2 = time2.round(decimals=5)
    time2 = time2.drop_duplicates().to_numpy()
    date = np.floor(date).unique()

    xlm1 = (dt.datetime.combine(time.date[0],dt.datetime.min.time())-tstart)/dt.timedelta(days=1)
    xlm1 = np.floor(xlm1)
    xlm2 = (dt.datetime.combine(time.date[-1],dt.datetime.min.time())-tstart)/dt.timedelta(days=1)
    xlm2 = np.floor(xlm2)
    xlms = [xlm1,xlm2]

    zmat = z2.to_numpy().reshape(len(date),len(time2)).T

    return zmat,date,time2,xlms



def calculate_clear_sky_radiation(lat,lon,ts,te,deltat,tz=None,C=None):
    """Calculation of potential (clear-sky) incoming short wave radiation. Calculations are based on
    Masters, G. M.: The Solar Resource, in: Renewable and Efficient Electric Power Systems, 
    385–443, https://doi.org/10.1002/0471668826.ch7, 2004.

    Args:
        lat (float64): latitude
        lon (float64): longitude
        ts (datetime): start of the period for which radiation is calculated
        te (datetime): end of the period for which radiation is calculated
        deltat (timedelta): time step
        tz (datetime, optional): timezone. Defaults to None.
        C (Series, optional): Time series for sky diffuse factor. Defaults to None.

    Returns:
        Series: Potential (clear-sky) incoming short wave radiation on flat surface including direct beam and estimate for diffuse radiation
    """    

    if tz is None:
        tz = dt.timezone.utc

    t = pd.date_range(ts,te,freq=deltat,tz=tz)
    alt = pysolar.solar.get_altitude_fast(lat,lon,t+deltat/2)
    alt[alt<0] = np.nan
    # direct beam radiation on a plane normal to the beam
    DBnorm = pd.Series(data=pysolar.radiation.get_radiation_direct(t+deltat/2,alt),index=t)
    DBnorm[DBnorm.isna()] = 0
    # direct beam radiation on earth surface, Masters, p. 414
    DBsurf = DBnorm*np.sin(alt*np.pi/180)
    DBsurf[DBsurf.isna()] = 0
    # sky diffuse factor
    if C is None:
        day = pysolar.numeric.tm_yday(t)
        # Masters, p. 416
        C = 0.095+0.04*np.sin(2*np.pi/365*(day-100))
    # estimate for the diffuse radiation on earth surface, Masters p. 416
    diffsurf = DBnorm*C

    SWINPOT = DBsurf+diffsurf
    SWINPOT.name = 'SWIN_POT'
    df = SWINPOT.reset_index()
    df = df.rename(mapper={'index':'time'},axis=1)
    return df
