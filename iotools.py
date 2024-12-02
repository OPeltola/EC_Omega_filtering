import os
import yaml
import pandas as pd
import gzip
import h5py
import numpy as np
import datetime as dt
import glob

import utils
import constants as c

def read_processed_data(filein):
    if os.path.isfile(filein):
        datain = pd.read_csv(filein,header=0)
        datain[datain==c.missing_value] = np.nan
        datain['TIME'] = pd.to_datetime(datain['TIME'],format='%Y-%m-%d %H:%M:%S')
        datain = datain.set_index('TIME')
    else:
        print('File ' + filein + ' not found.')
        datain = pd.DataFrame()
    return datain

def output_processed_data(df,file):
    dfout = df.copy()
    if 'TIME' in dfout.columns:
        dfout.drop('TIME', axis=1, inplace=True)
    dfout[dfout.isna()] = c.missing_value

    dfout.to_csv(file)




def load_ICOS_LAI(DATADIR,siteID):
    """Loading LAI from ICOS sites

    Args:
        DATADIR (str): directory for the LAI files
        siteID (str): Site ID

    Returns:
        DataFrame: Leaf area index and metadata
    """    
    files = next(os.walk(os.path.join(DATADIR,siteID)))[2]
    files = [i for i in files if siteID in i]
    files = [i for i in files if '_ANCILLARY.csv' in i]

    LAI = pd.Series(dtype=float)
    if len(files)>0:
        fin = os.path.join(DATADIR,siteID,files[0])
        try:
            df = pd.read_csv(fin, sep=',', header=0,na_values=-9999)
        except:
            print('Problem loading ' + files[0] + '. Trying with different encoding')
            df = pd.read_csv(fin, sep=',', header=0,na_values=-9999, encoding='latin1')
        df = df[df['VARIABLE_GROUP']=='GRP_LAI']

        for gid in df['GROUP_ID'].unique():
            dftmp = df[df['GROUP_ID']==gid]
            if (dftmp.loc[dftmp['VARIABLE']=='LAI_STATISTIC','DATAVALUE'].iloc[0]=='Mean'):
                if (dftmp['VARIABLE']=='LAI_DATE').any():
                    datestr = df.loc[(df['GROUP_ID']==gid) & (df['VARIABLE']=='LAI_DATE'),'DATAVALUE'].astype('str').iloc[0]
                elif (dftmp['VARIABLE']=='LAI_DATE_START').any():
                    datestr = df.loc[(df['GROUP_ID']==gid) & (df['VARIABLE']=='LAI_DATE_START'),'DATAVALUE'].astype('str').iloc[0]

                if len(datestr)==8:
                    time = pd.to_datetime(datestr,format='%Y%m%d')
                else:
                    time = pd.to_datetime(datestr,format='%Y%m%d%H%M')

                val = df.loc[(df['GROUP_ID']==gid) & (df['VARIABLE']=='LAI'),'DATAVALUE'].astype('float').iloc[0]

                LAI = pd.concat([LAI,pd.Series(index=[time],data=val)])

    LAI = LAI.sort_index()
    return LAI


def load_NEON_LAI(DATADIR,siteID):
    """Loading LAI from NEON sites, source: https://gbov.acri.fr/home/

    Args:
        DATADIR (str): directory for the LAI files
        siteID (str): Site ID

    Returns:
        DataFrame: Leaf area index and metadata
    """    

    sname = utils.NEON_site_name(siteID)

    LAI = pd.DataFrame()
    if len(sname)>0:
        files = next(os.walk(DATADIR))[2]
        files = [i for i in files if sname in i]
        files = [i for i in files if '.csv' in i]

        for file in files:
            tmp = pd.read_csv(os.path.join(DATADIR,file),sep=';')
            LAI = pd.concat([LAI,tmp])

        if 'TIME_IS' in LAI.columns:
            LAI['TIME'] = pd.to_datetime(LAI['TIME_IS'].astype('str')).dt.tz_localize(None)
            LAI = LAI.set_index('TIME',drop=False)
            LAI = LAI.sort_index()

    LAI[LAI==-999] = np.nan

    return LAI




def load_hdf5_dataset(dataset,field):
    dat = pd.DataFrame(dataset[()])
    dat['TIME'] = pd.to_datetime(dat['timeBgn'].astype('str')).dt.tz_localize(None)
    dat = dat[['TIME',field]]
    dat = dat.set_index('TIME')
    datout = dat.squeeze(axis=1)
    return datout


def load_NEON_data(DATADIR,site):
    """Loading NEON daa

    Args:
        DATADIR (str): directory for the data
        site (str): site ID
        TMPDIR (str): directory for temporary data
    """    
    dfout = pd.DataFrame()
    dfzout = pd.Series()
    site_info = dict()


    subdirs = next(os.walk(DATADIR))[1]
    subdirs = [i for i in subdirs if site in i]
    for subdir in subdirs:
        curddir = DATADIR+subdir
        gzfile = [f for f in os.listdir(curddir) if (os.path.isfile(os.path.join(curddir, f)) & ('.gz' in f))]
        if len(gzfile)>0:
            print('Loading ' + gzfile[0])
            with gzip.open(os.path.join(curddir, gzfile[0])) as fgz:
                with h5py.File(fgz, 'r') as fh5:
                    dftmp = pd.DataFrame()

                    attributes = list(fh5[site].attrs.keys())
                    # canopy height
                    site_info['canz'] = fh5[site].attrs['DistZaxsCnpy'].astype(float)[0]
                    # displacement height
                    site_info['d'] = fh5[site].attrs['DistZaxsDisp'].astype(float)[0]
                    # latitude
                    site_info['lat'] = fh5[site].attrs['LatTow'].astype(float)[0]
                    # longitude
                    site_info['lon'] = fh5[site].attrs['LonTow'].astype(float)[0]
                    # EC height
                    site_info['z'] = fh5[site].attrs['DistZaxsTow'].astype(float)[0]
                    # T heights
                    site_info['zT'] = fh5[site].attrs['DistZaxsLvlMeasTow'].astype(float)

                    # info = pd.DataFrame(fh5['objDesc'][()])
                    # info['Object'] = info['Object'].astype(str)
                    # info['Description'] = info['Description'].astype(str)
                    
                    # fluxes
                    dftmp['FC'] = load_hdf5_dataset(fh5[site]['dp04']['data']['fluxCo2']['turb'],'flux')
                    dftmp['NEE'] = load_hdf5_dataset(fh5[site]['dp04']['data']['fluxCo2']['nsae'],'flux')
                    dftmp['LE'] = load_hdf5_dataset(fh5[site]['dp04']['data']['fluxH2o']['nsae'],'flux')
                    dftmp['USTAR'] = load_hdf5_dataset(fh5[site]['dp04']['data']['fluxMome']['turb'],'veloFric')
                    dftmp['H'] = load_hdf5_dataset(fh5[site]['dp04']['data']['fluxTemp']['turb'],'flux')

                    # sonic data
                    groups = list(fh5[site]['dp01']['data']['soni'])
                    sonic_height = groups[1]
                    dftmp['WD'] = load_hdf5_dataset(fh5[site]['dp01']['data']['soni'][groups[1]]['angZaxsErth'],'mean')
                    dftmp['WS'] = load_hdf5_dataset(fh5[site]['dp01']['data']['soni'][groups[1]]['veloXaxsYaxsErth'],'mean')
                    dftmp['W_SIGMA'] = load_hdf5_dataset(fh5[site]['dp01']['data']['soni'][groups[1]]['veloZaxsErth'],'vari')**0.5
                    dftmp['U_SIGMA'] = load_hdf5_dataset(fh5[site]['dp01']['data']['soni'][groups[1]]['veloXaxsErth'],'vari')**0.5
                    dftmp['V_SIGMA'] = load_hdf5_dataset(fh5[site]['dp01']['data']['soni'][groups[1]]['veloYaxsErth'],'vari')**0.5
                    dftmp['T_SIGMA'] = load_hdf5_dataset(fh5[site]['dp01']['data']['soni'][groups[1]]['tempSoni'],'vari')**0.5

                    dfzout['FC'] = site_info['z']
                    dfzout['NEE'] = site_info['z']
                    dfzout['LE'] = site_info['z']
                    dfzout['USTAR'] = site_info['z']
                    dfzout['H'] = site_info['z']
                    dfzout['WD'] = site_info['z']
                    dfzout['WS'] = site_info['z']
                    dfzout['W_SIGMA'] = site_info['z']
                    dfzout['U_SIGMA'] = site_info['z']
                    dfzout['V_SIGMA'] = site_info['z']
                    dfzout['T_SIGMA'] = site_info['z']

                    # air pressure
                    groups = list(fh5[site]['dp01']['data']['presBaro'])
                    dftmp['PA'] = load_hdf5_dataset(fh5[site]['dp01']['data']['presBaro'][groups[1]]['presAtm'],'mean')
                    dfzout['PA'] = site_info['z']

                    # temperature profile
                    groups = list(fh5[site]['dp01']['data']['tempAirLvl'])
                    groups = [i for i in groups if '_30m' in i]
                    for i in range(len(groups)):
                        dftmp['TA_'+str(i)] = load_hdf5_dataset(fh5[site]['dp01']['data']['tempAirLvl'][groups[i]]['temp'],'mean')
                        dfzout['TA_'+str(i)] = site_info['zT'][i]
                    groups = list(fh5[site]['dp01']['data']['tempAirTop'])
                    dftmp['TA_'+str(i+1)] = load_hdf5_dataset(fh5[site]['dp01']['data']['tempAirTop'][groups[1]]['temp'],'mean')
                    dfzout['TA_'+str(i+1)] = site_info['zT'][i+1]

                    # humidity profile
                    groups = list(fh5[site]['dp01']['data']['h2oStor'])
                    groups = [i for i in groups if '_30m' in i]
                    groups = [i for i in groups if '000' in i]
                    for i in range(len(groups)):
                        dftmp['H2O_'+str(i)] = load_hdf5_dataset(fh5[site]['dp01']['data']['h2oStor'][groups[i]]['rtioMoleDryH2o'],'mean')
                        dfzout['H2O_'+str(i)] = site_info['zT'][i]

                    # CO2 profile
                    groups = list(fh5[site]['dp01']['data']['co2Stor'])
                    groups = [i for i in groups if '_30m' in i]
                    groups = [i for i in groups if '000' in i]
                    for i in range(len(groups)):
                        dftmp['CO2_'+str(i)] = load_hdf5_dataset(fh5[site]['dp01']['data']['co2Stor'][groups[i]]['rtioMoleDryCo2'],'mean')
                        dfzout['CO2_'+str(i)] = site_info['zT'][i]
                    
                    # radiation components
                    groups = list(fh5[site]['dp01']['data']['radiNet'])
                    dftmp['LW_IN'] = load_hdf5_dataset(fh5[site]['dp01']['data']['radiNet'][groups[1]]['radiLwIn'],'mean')
                    dftmp['LW_OUT'] = load_hdf5_dataset(fh5[site]['dp01']['data']['radiNet'][groups[1]]['radiLwOut'],'mean')
                    dftmp['SW_IN'] = load_hdf5_dataset(fh5[site]['dp01']['data']['radiNet'][groups[1]]['radiSwIn'],'mean')
                    dftmp['SW_OUT'] = load_hdf5_dataset(fh5[site]['dp01']['data']['radiNet'][groups[1]]['radiSwOut'],'mean')
                    dfzout['LW_IN'] = site_info['z']
                    dfzout['LW_OUT'] = site_info['z']
                    dfzout['SW_IN'] = site_info['z']
                    dfzout['SW_OUT'] = site_info['z']


                    dfout = pd.concat([dfout,dftmp])

                    del dftmp
                    # co2flx2 = pd.DataFrame(fh5[site]['dp04']['data']['fluxCo2']['nsae'][()])

                    # co2flx = pd.DataFrame(fh5[site]['dp04']['data']['fluxCo2']['turb'][()])
                    # co2flx2 = pd.DataFrame(fh5[site]['dp04']['data']['fluxCo2']['nsae'][()])

                    # = fh5[site]['dp04']['data']['radiNet']['000_040_30m']['radiLwIn'][()]

                    # data = list(fh5[site]['dp01']['data']['radiNet']['000_040_30m']['radiLwIn'])
                    # data = fh5[site]['dp01']['data']['radiNet']['000_040_30m']['radiLwIn'][()]


            # print('Unzipping ' + os.path.join(curddir, gzfile[0]))
            # try:
            #     with ZipFile(os.path.join(curddir, gzfile[0]), 'r') as zipfl:
            #         uzfile = gzfile[0].replace('.gz','')
            #         zipfl.extract(uzfile,TMPDIR)
            #         # zipfl.extractall(TMPDIR)        
            # except Exception as e:
            #     print(e)

    if dfout.empty:
        print('Data not found for ' + site)
    return dfout,dfzout,site_info

def load_FI_Jok_data(DATADIR):

    dfout = pd.DataFrame()
    dfzout = pd.Series()

    filein = os.path.join(DATADIR,'FI_Jok.txt')

    if os.path.isfile(filein):
        dfin = pd.read_csv(filein, sep=';', header=None,na_values=-9999)

        tvec = dfin.apply(lambda row: dt.datetime(int(row[0]), int(row[1]), int(row[2]),int(row[3]),int(row[4]),0), axis=1)
        dfout['TIME'] = tvec
        cols = {'H':8,'TAU':9,'LE':12,'FC':11,'USTAR':10,'W_SIGMA':13,'ZL':7,'WS':5,'WD':6,'TA_1_1_1':19,'TA_2_1_1':20,'SW_IN':16,'R_NET':17,'PPFD_IN':18,'RH_1_1_1':25,'RH_2_1_1':26}
        for col in cols.keys():
            dfout[col] = dfin[cols[col]].copy()
        dfout['SC'] = np.nan
        dfout['NEE'] = np.nan

        dfzout['H'] = 3
        dfzout['TAU'] = 3
        dfzout['LE'] = 3
        dfzout['FC'] = 3
        dfzout['USTAR'] = 3
        dfzout['W_SIGMA'] = 3
        dfzout['ZL'] = 3
        dfzout['WS'] = 3
        dfzout['WD'] = 3
        dfzout['TA_1_1_1'] = 2.5
        dfzout['TA_2_1_1'] = 0.2
        dfzout['RH_1_1_1'] = 25
        dfzout['RH_2_1_1'] = 0.2
        dfzout['SW_IN'] = 3
        dfzout['R_NET'] = 3
        dfzout['PPFD_IN'] = 3


        dfout['MO_LENGTH'] = dfzout['H']/dfout['ZL']


        dfout = dfout.set_index('TIME',drop=False)

        # tvec_new = pd.date_range(start=dfout.index[0],end=dfout.index[-1],freq='30min')
    return dfout,dfzout

def load_FI_Ken_data(DATADIR):

    dfout = pd.DataFrame()
    dfzout = pd.Series()

    filein1 = os.path.join(DATADIR,'Kenttärova ustar test 2018-2020.xlsx')
    filein2 = os.path.join(DATADIR,'Kenttärova ustar test B 2018-2020.xlsx')

    if os.path.isfile(filein1):
        dfin = pd.read_excel(filein1, index_col=0)
        dfin.index.names = ['TIME']

        
        cols = {'H':'Sensible heat flux (W m-2)','LE':'Latent heat flux (W m-2)','FC':'CO2 flux (mgCO2 m-2 s-1) (kevyt karsinta)',\
                'USTAR':'friction velocity (m/s)','W_SIGMA':"rot<w'w'> ((m s-1)^2)",'SC':'Storage flux (mgCO2/m2s)'}

        for col in cols.keys():
            dfout[col] = dfin[cols[col]].copy()
            
        dfout['FC'] = dfout['FC']/c.Mco2
        dfout['SC'] = dfout['SC']/c.Mco2
        dfout['NEE'] = dfout['FC'] + dfout['SC']
        dfout['W_SIGMA'] = dfout['W_SIGMA']**0.5
        dfout['MO_LENGTH'] = 1/dfin['inverse Obukhov length: L-1 (m-1)']

      
    if os.path.isfile(filein2):
        dfin = pd.read_excel(filein2, index_col=0)

        
        cols = {'TA_1_1_1':'t_2m','TA_2_1_1':'t_10m','TA_3_1_1':'t_20m',\
                'RH_1_1_1':'rh_2m','WS':'wspeed','WD':'wdir','PA':'press','PPFD_IN':'PPFD','SW_IN':'avg(Global rad)'}

        for col in cols.keys():
            dfout[col] = dfin[cols[col]].copy()
            

    dfzout['H'] = 23
    dfzout['LE'] = 23
    dfzout['FC'] = 23
    dfzout['USTAR'] = 23
    dfzout['W_SIGMA'] = 23
    dfzout['ZL'] = 23
    dfzout['WS'] = 23
    dfzout['WD'] = 23
    dfzout['TA_1_1_1'] = 2
    dfzout['TA_2_1_1'] = 10
    dfzout['TA_3_1_1'] = 20
    dfzout['RH_1_1_1'] = 2
    dfzout['SW_IN'] = 24
    dfzout['PPFD_IN'] = 24
    dfzout['PA'] = 24


    dfout = dfout.asfreq(freq='30min', fill_value=np.nan)
    dfout['TIME'] = dfout.index


        # tvec_new = pd.date_range(start=dfout.index[0],end=dfout.index[-1],freq='30min')
    return dfout,dfzout



def load_FI_Sod_data(DATADIR):

    dfout = pd.DataFrame()
    dfzout = pd.Series()

    filein1 = os.path.join(DATADIR,'Sodankylä ustar test 2018-2020.xlsx')
    filein2 = os.path.join(DATADIR,'MET0002_2006-01-01_2012-12-31_100521095442.txt')

    if os.path.isfile(filein1):
        dfin = pd.read_excel(filein1, index_col=0)
        dfin.index = dfin.index-dt.timedelta(minutes=30)
        dfin = dfin[dfin.index>=dt.datetime(2006,1,1)]
        dfin.index.names = ['TIME']

        
        cols = {'H':'Sensible heat flux (W m-2)','LE':'Latent heat flux (W m-2)','FC':'CO2 flux (mgCO2 m-2 s-1) (kevyt karsinta)',\
                'USTAR':'friction velocity (m/s)','W_SIGMA':"<w'w'> (unrotated)",'SC':'CO2 storage flux (mgm-2s-1)','TA_1_1_1':'T 3m (degC)',\
                'TA_2_1_1':'T 8m (degC)','TA_3_1_1':'T 18m (degC)','TA_4_1_1':'T 32m (degC)','TA_5_1_1':'T 48m (degC)','SW_IN':'Global radiation (W/m2) (Luotaamo)',\
                'PPFD_IN':'PPFD (umol/m2s)','WS':'wind speed (m/s) EC','WD':'wind dir EC'}

        for col in cols.keys():
            dfout[col] = dfin[cols[col]].copy()
            
        dfout['FC'] = dfout['FC']/c.Mco2
        dfout['SC'] = dfout['SC']/c.Mco2
        dfout['NEE'] = dfout['FC'] + dfout['SC']
        dfout['W_SIGMA'] = dfout['W_SIGMA']**0.5
        dfout['MO_LENGTH'] = 1/dfin['inverse Obukhov length: L-1 (m-1)']
        
        dfout['PA'] = 101.325*np.exp(-179/8431)

            

    dfzout['H'] = 24.5
    dfzout['LE'] = 24.5
    dfzout['FC'] = 24.5
    dfzout['USTAR'] = 24.5
    dfzout['W_SIGMA'] = 24.5
    dfzout['ZL'] = 24.5
    dfzout['WS'] = 24.5
    dfzout['WD'] = 24.5
    dfzout['TA_1_1_1'] = 3
    dfzout['TA_2_1_1'] = 8
    dfzout['TA_3_1_1'] = 18
    dfzout['TA_4_1_1'] = 32
    dfzout['TA_5_1_1'] = 48
    dfzout['SW_IN'] = 24
    dfzout['PPFD_IN'] = 24
    dfzout['PA'] = 24


    dfout = dfout.asfreq(freq='30min', fill_value=np.nan)
    dfout['TIME'] = dfout.index


        # tvec_new = pd.date_range(start=dfout.index[0],end=dfout.index[-1],freq='30min')
    return dfout,dfzout

def load_US_MRf_data(DATADIR):

    dfout = pd.DataFrame()
    dfzout = pd.Series()
    dfflxout = pd.DataFrame()

    filein = os.path.join(DATADIR,'US-Fir_2006-2011_fluxmet_gapfilled_ct_07Feb13.csv')

    if os.path.isfile(filein):
        dfflxin = pd.read_csv(filein, sep=',',skiprows=[7],header=[6],na_values='NaN')
        dfflxin[dfflxin==-999] = np.nan

        tvec = (pd.to_datetime(dfflxin['YEAR'] * 1000 + dfflxin['DOY'], format='%Y%j') +
           pd.to_timedelta(np.floor(dfflxin['HRMIN']*1e-2), unit='h') +
           pd.to_timedelta(dfflxin['HRMIN']-np.floor(dfflxin['HRMIN']*1e-2)*1e2-15, unit='min'))
        
        dfflxout['TIME'] = tvec

        cols = {'H':'H_38_CSAT3','LE':'LE_38_CSAT3&LI7000&LI7500','FC':'FC_38_CSAT&LI7000&LI7500',\
                'USTAR':'UST_38_CSAT','W_SIGMA':'SIGMAW_38_CSAT','WS':'WS_38_CSAT&CUP','WD':'WD_38_CSAT','TA_1_1_1':'TA_38_PRT1000&HMP',\
                'SW_IN':'Rg_37_CNR1','R_NET':'RNET_37_Q7&NRLITE&CNR1','PPFD_IN':'PAR_37_LI190&PARlite','RH_1_1_1':'RH_38_HMP','RH_2_1_1':'RH_4_HMP',\
                'SC':'SFC_38_Li840&Li7500','FCqf':'FCflag_38',\
                'NEE':'NEE_38_CSAT3&LI7000&LI7500','NEEqf':'NEEflag_38','VPD_1_1_1':'VPD_37_HMP&LI7000&LI7500',\
                'PA':'PA_37','H2O':'H2O_37_LI7000&Li7500'}
        for col in cols.keys():
            dfflxout[col] = dfflxin[cols[col]].copy()

        H2Omr = dfflxout['H2O'].copy()
        H2Omr = H2Omr*(c.R*(dfflxout['TA_1_1_1']+c.NT)/(dfflxout['PA']*1e3))
        H2Omr[pd.isnull(H2Omr)] = H2Omr.median()
        dfflxout['MO_LENGTH'] = utils.calculate_Obukhov_length(dfflxout['H'],dfflxout['USTAR'],dfflxout['PA'],dfflxout['TA_1_1_1'],H2O=H2Omr)
        
        # removing gapfilled data            
        dfflxout.loc[dfflxout['FCqf']==1,'FC'] = np.nan
        dfflxout.loc[dfflxout['NEEqf']==1,'NEE'] = np.nan

        dfflxout = dfflxout.set_index('TIME')


    filein = os.path.join(DATADIR,'MF_temperature_profile_corrected_2007_DOY001to074_30minavg.csv')
    if os.path.isfile(filein):
        dfTin1 = pd.read_csv(filein, sep=',',header=[2],na_values='NaN')
        dfTin1[dfTin1==-999] = np.nan
        tvec = (pd.to_datetime(dfTin1['yyyy'] * 1000 + dfTin1['DOY_center'].astype(int), format='%Y%j') +
           pd.to_timedelta(np.floor(dfTin1['HH']), unit='h') +
           pd.to_timedelta(dfTin1['MM']-15, unit='min'))
        dfTin1['TIME'] = tvec
        dfTin1 = dfTin1.set_index('TIME',drop=False)
    else:
        dfTin1 = pd.DataFrame()

    filein = os.path.join(DATADIR,'MF_temperature_profile_corrected_2007_DOY074to298_30minavg.csv')
    if os.path.isfile(filein):
        dfTin2 = pd.read_csv(filein, sep=',',header=[2],na_values='NaN')
        dfTin2[dfTin2==-999] = np.nan
        tvec = (pd.to_datetime(dfTin2['Year'] * 1000 + dfTin2['DOY_center'].astype(int), format='%Y%j') +
           pd.to_timedelta(np.floor(dfTin2['hh']), unit='h') +
           pd.to_timedelta(dfTin2['mm']-15, unit='min'))
        dfTin2['TIME'] = tvec
        dfTin2 = dfTin2.set_index('TIME',drop=False)
    else:
        dfTin2 = pd.DataFrame()

    filein = os.path.join(DATADIR,'MF_temperature_profile_corrected_2006_DOY306to365_30minavg.csv')
    if os.path.isfile(filein):
        dfTin3 = pd.read_csv(filein, sep=',',header=[2],na_values='NaN')
        dfTin3[dfTin3==-999] = np.nan
        tvec = (pd.to_datetime(dfTin3['yyyy'] * 1000 + dfTin3['DOY_center'].astype(int), format='%Y%j') +
           pd.to_timedelta(np.floor(dfTin3['HH']), unit='h') +
           pd.to_timedelta(dfTin3['MM']-15, unit='min'))
        dfTin3['TIME'] = tvec
        dfTin3 = dfTin3.set_index('TIME',drop=False)
    else:
        dfTin3 = pd.DataFrame()

    cols = {'MF0_hobo112':'TA_2_1_1','MF1_hobo62':'TA_3_1_1','MF2_hobo65':'TA_4_1_1',\
            'MF3_hobo73':'TA_5_1_1','MF4_hobo75':'TA_6_1_1','MF5_hobo90':'TA_7_1_1',\
            'MF6_hobo95':'TA_8_1_1','MF7_097':'TA_9_1_1'}
    dfTcomb = pd.concat([dfTin1,dfTin2,dfTin3])
    dfTcomb = dfTcomb.rename(columns=cols)
    dfTcomb = dfTcomb.sort_index()
    dfTcomb = dfTcomb[~dfTcomb.index.duplicated(keep='first')]

    dfout = pd.DataFrame(index=dfTcomb.index)
    dfout[['TA_2_1_1','TA_3_1_1','TA_4_1_1','TA_5_1_1','TA_6_1_1','TA_7_1_1','TA_8_1_1','TA_9_1_1']] = \
    dfTcomb[['TA_2_1_1','TA_3_1_1','TA_4_1_1','TA_5_1_1','TA_6_1_1','TA_7_1_1','TA_8_1_1','TA_9_1_1']]

    dfout = pd.merge(dfout,dfflxout,left_index=True,right_index=True,how='inner')

    dfzout['H'] = 38
    dfzout['LE'] = 38
    dfzout['FC'] = 38
    dfzout['USTAR'] = 38
    dfzout['W_SIGMA'] = 38
    # dfzout['ZL'] = 38
    dfzout['MO_LENGTH'] = 38
    dfzout['WS'] = 38
    dfzout['WD'] = 38
    dfzout['TA_1_1_1'] = 38
    dfzout['TA_2_1_1'] = 0.5
    dfzout['TA_3_1_1'] = 2
    dfzout['TA_4_1_1'] = 4.5
    dfzout['TA_5_1_1'] = 9.3
    dfzout['TA_6_1_1'] = 13
    dfzout['TA_7_1_1'] = 19
    dfzout['TA_8_1_1'] = 26.8
    dfzout['TA_9_1_1'] = 37.5
    dfzout['RH_1_1_1'] = 38
    dfzout['RH_2_1_1'] = 4
    dfzout['SW_IN'] = 37
    dfzout['R_NET'] = 37
    dfzout['PPFD_IN'] = 37

    dfout = dfout.asfreq(freq='30min', fill_value=np.nan)
    dfout['TIME'] = dfout.index


    return dfout,dfzout


def load_ICOS_data(DATADIR,site):
    """Loading ICOS data

    Args:
        DATADIR (str): directory for the data
        site (str): site ID
    """    
    dfout = pd.DataFrame()
    dfzout = pd.Series()

    # file containing fluxes
    flxfile = DATADIR + site + '/ICOSETC_' + site + '_FLUXES_01.csv'
    # file containing meteo
    metfile = DATADIR + site + '/ICOSETC_' + site + '_METEO_01.csv'
    # file containing flux info
    flxinfofile = DATADIR + site + '/ICOSETC_' + site + '_VARINFO_FLUXES_01.csv'
    # file containing meteo info
    metinfofile = DATADIR + site + '/ICOSETC_' + site + '_VARINFO_METEO_01.csv'

    if os.path.isfile(flxfile) & os.path.isfile(metfile):
        print('Loading ' + flxfile)
        dfflx = pd.read_csv(flxfile, sep=',', header=0,na_values=-9999)
        print('Loading ' + metfile)
        dfmet = pd.read_csv(metfile, sep=',', header=0,na_values=-9999)
        
        print('Loading ' + flxinfofile)
        dfinfo = pd.read_csv(flxinfofile, sep=',', header=0,na_values=-9999)
        dfflxz = pd.Series()
        for gid in dfinfo['GROUP_ID'].unique():
            if (dfinfo.loc[(dfinfo['GROUP_ID']==gid)]['VARIABLE']=='VAR_INFO_HEIGHT').any():
                dfflxz[dfinfo.loc[(dfinfo['GROUP_ID']==gid) & (dfinfo['VARIABLE']=='VAR_INFO_VARNAME')]['DATAVALUE'].iloc[0]] = dfinfo.loc[(dfinfo['GROUP_ID']==gid) & (dfinfo['VARIABLE']=='VAR_INFO_HEIGHT')]['DATAVALUE'].astype(float).iloc[0]
        
        print('Loading ' + metinfofile)
        dfinfo = pd.read_csv(metinfofile, sep=',', header=0,na_values=-9999)
        dfmetz = pd.Series()
        for gid in dfinfo['GROUP_ID'].unique():
            if (dfinfo.loc[(dfinfo['GROUP_ID']==gid)]['VARIABLE']=='VAR_INFO_HEIGHT').any():
                dfmetz[dfinfo.loc[(dfinfo['GROUP_ID']==gid) & (dfinfo['VARIABLE']=='VAR_INFO_VARNAME')]['DATAVALUE'].iloc[0]] = dfinfo.loc[(dfinfo['GROUP_ID']==gid) & (dfinfo['VARIABLE']=='VAR_INFO_HEIGHT')]['DATAVALUE'].astype(float).iloc[0]


        dfflx['TIME'] = pd.to_datetime(dfflx['TIMESTAMP_START'].astype('str'),format='%Y%m%d%H%M')
        dfflx = dfflx.set_index('TIME',drop=False)
        dfflx = dfflx.sort_index()
        dfflx = dfflx.drop(['TIMESTAMP_START','TIMESTAMP_END'],axis=1)
        
        dfmet['TIME'] = pd.to_datetime(dfmet['TIMESTAMP_START'].astype('str'),format='%Y%m%d%H%M')
        dfmet = dfmet.set_index('TIME',drop=False)
        dfmet = dfmet.sort_index()
        dfmet = dfmet.drop(['TIMESTAMP_START','TIMESTAMP_END'],axis=1)
        dfmet = dfmet.drop(dfmet.columns[dfmet.columns.str.endswith('_N')],axis=1)
        dfmet = dfmet.drop(dfmet.columns[dfmet.columns.str.endswith('_SD')],axis=1)

        dfzout = pd.concat([dfmetz,dfflxz])
        
        dfout = pd.merge(dfflx.reset_index(drop=True), dfmet.reset_index(drop=True), on='TIME',how='outer',suffixes=[None,'_met'])
        dfout = dfout.set_index('TIME',drop=False)
        dfout = dfout.sort_index()
    else:
        print('Data files missing for '+site)


    if ('WS' not in dfout.columns) & (not dfout.empty):
        # if no WS data, then taking the closest height to EC
        WScols = [i for i in dfout.columns if 'WS_' in i]
        boolean_for_col = (dfzout['FC']-dfzout[WScols]).abs()==(dfzout['FC']-dfzout[WScols]).abs().min()
        WScol = boolean_for_col.index[boolean_for_col][0]
        dfzout['WS'] = dfzout[WScol]
        dfout['WS'] = dfout[WScol]
        dfzout['WD'] = dfzout[WScol.replace('WS','WD')]
        dfout['WD'] = dfout[WScol.replace('WS','WD')]
    if ('FC' not in dfzout.index) & (not dfout.empty):
        dfzout['FC'] = dfzout['USTAR']
        dfzout['NEE'] = dfzout['USTAR']

    if dfout.empty:
        print('Data not found for ' + site)

    if 'NEE' not in dfout.columns:
        dfout['NEE'] = np.nan

    return dfout,dfzout


def load_config(fin):
    """Loads config-file.

    Args:
        fin (str): file name

    Returns:
        dict: calculation configuration
    """    

    config = dict()
    if os.path.isfile(fin):
        with open(fin, 'r') as stream:
            try:
                config = yaml.safe_load(stream)
            except Exception as exc:
                print('Error when reading file ' + fin)
                print(exc)
    else:
        print('File ' + fin + ' does not exist.')

    return config

def load_ONEflux_threshold(DATADIR,site,flag):

    USTARth = pd.DataFrame()

    ONEFLUXDIR = os.path.join(DATADIR,'ONEflux',flag+'_USTARth')
    fnameformat = site.replace('-','_')+'_'+flag+'_%Y_ut.txt'
    fnameformat = site+'_%Yrid_ut.txt'

    fname, fformat_ext = os.path.splitext(fnameformat)
    fname = os.path.join(ONEFLUXDIR,fnameformat)

    now = dt.datetime.now()
    namelen = len(now.strftime(fnameformat))-len(fformat_ext)
    msk=os.path.join(ONEFLUXDIR,'?'*namelen + '%s'%(fformat_ext))
    msk=os.path.join(ONEFLUXDIR,site+'_????rid_ut.txt')
    fls=glob.glob(msk)
    tstamps = pd.to_datetime(fls,format=fname,errors='coerce')
    yrs = tstamps.year.to_list()

    for indx in range(len(fls)):
        fl = fls[indx]
        with open(fl) as fid:
            lines=fid.readlines()
        start_row = [i for i in range(len(lines)) if lines[i]=='-- percentiles forward mode 2\n']
        if len(start_row)==1:
            if (len(lines)>start_row[0]+10):
                vals = []
                for row in np.linspace(start_row[0]+1,start_row[0]+9,9):
                    val = lines[int(row)]
                    val = val.replace('\n','')
                    val = val.replace(' 50%','')
                    vals.append(float(val))
            else:
                vals = (np.ones((9,))*np.nan).tolist()
            # USTARth[tstamps[indx]] = vals
        else:
            vals = (np.ones((9,))*np.nan).tolist()
        
        USTARth[yrs[indx]] = vals



    return USTARth