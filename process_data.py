import numpy as np
import pandas as pd
import datetime as dt
import os
import matplotlib.pyplot as plt
from dotenv import load_dotenv;load_dotenv()

import iotools
import utils
import constants as c

config_file = 'config.yml'
config = iotools.load_config(config_file)    
    
ICOSDATADIR = os.getenv('ICOSDATADIR')
FIJOKDATADIR = os.getenv('FIJOKDATADIR')
USMRFDATADIR = os.getenv('USMRFDATADIR')
FIKENDATADIR = os.getenv('FIKENDATADIR')
FISODDATADIR = os.getenv('FISODDATADIR')
NEONDATADIR = os.getenv('NEONDATADIR')
ONEFLUXDIR = os.getenv('ONEFLUXDIR')
NEONLAIDIR = os.getenv('NEONLAIDIR')


OMEGAth = 0.66

def main():

    for site in config['sites']:

        # loading
        if config['source'][site]=='ICOS':
            data,z = iotools.load_ICOS_data(ICOSDATADIR,site)

            LAI = iotools.load_ICOS_LAI(ICOSDATADIR,site)

            if site=='FI-Hyy':
                # Hyytiälä LAI based on Aslan et al. 
                LAI = pd.Series(index=data.index)
                LAI[LAI.index<dt.datetime(2019,4,1)] = 3.9
                LAI[(LAI.index>=dt.datetime(2019,4,1)) & (LAI.index<dt.datetime(2020,1,1))] = 3.6
                LAI[LAI.index>=dt.datetime(2020,1,1)] = 2.1

            if len(LAI)>0:
                data['LAI'] = utils.interpolate_LAI(LAI,data.index)
            else:
                data['LAI'] = config['LAI'][site]


            if site in config['canz'].keys():
                data['canz'] = config['canz'][site]
                data['d'] = config['canz'][site]*0.66
            else:
                print('Canopy height missing for ' + site)
                data['canz'] = 0
                data['d'] = 0

        elif config['source'][site]=='NEON':

            data,z,site_info = iotools.load_NEON_data(NEONDATADIR,site)
            if not data.empty:
                H2Okeys = [i for i in list(z.keys()) if 'H2O_' in i]
                boolean_for_key = np.abs((z['FC']-z[H2Okeys]))==np.min(np.abs((z['FC']-z[H2Okeys])))
                H2Okey = boolean_for_key.index[boolean_for_key][0]
                del H2Okeys,boolean_for_key

                # calculating relative humidity
                data['RH'] = utils.calculate_RH(data[H2Okey],data['TA'+H2Okey[3:]],data['PA'])
                data.loc[data['RH']>100,'RH'] = 100
                z['RH'] = z['FC']
                # calculating Obukhov length and MO stability parameter
                data['MO_LENGTH'] = utils.calculate_Obukhov_length(data['H'],data['USTAR'],data['PA'],data['TA'+H2Okey[3:]],data[H2Okey])
                data['ZL'] = (site_info['z']-site_info['d'])/data['MO_LENGTH']
                z['MO_LENGTH'] = z['FC']
                z['ZL'] = z['FC']

                LAI = iotools.load_NEON_LAI(NEONLAIDIR,site)
                if not LAI.empty:
                    data['LAI'] = utils.interpolate_LAI(LAI[config['LAI_estimate']],data.index)
                    if pd.isnull(data['LAI']).sum()==len(data):
                        print('No good LAI values for ' + site + '. Assuming LAI=0.')
                        data['LAI'] = 0

                else:
                    data['LAI'] = 0
                data['canz'] = site_info['canz']
                data['d'] = site_info['d']

                # plt.plot(LAI[config['LAI_estimate']],'k.')
                # plt.plot(LAI['LAI_Miller_up'],'ro')
                # plt.show()
        elif config['source'][site]=='PI':

            if site=='FI-Jok':
                data,z = iotools.load_FI_Jok_data(FIJOKDATADIR)
                data['LAI'] = 0
                data['d'] = 0
                data['canz'] = 0
            elif site=='US-MRf':
                data,z = iotools.load_US_MRf_data(USMRFDATADIR)
                data['canz'] = config['canz'][site]
                data['d'] = 0.66*config['canz'][site]
                data['ZL'] = (z['FC']-data['d'])/data['MO_LENGTH']
                data['LAI'] = config['LAI'][site]
            elif site=='FI-Ken':
                data,z = iotools.load_FI_Ken_data(FIKENDATADIR)
                data['canz'] = config['canz'][site]
                data['d'] = 0.66*config['canz'][site]
                data['ZL'] = (z['FC']-data['d'])/data['MO_LENGTH']
                data['LAI'] = config['LAI'][site]
            elif site=='FI-Sod':
                data,z = iotools.load_FI_Sod_data(FISODDATADIR)
                data['canz'] = config['canz'][site]
                data['d'] = 0.66*config['canz'][site]
                data['ZL'] = (z['FC']-data['d'])/data['MO_LENGTH']
                data['LAI'] = config['LAI'][site]


            else:
                data = pd.DataFrame()
                print('Unidentified site. Not loading data.')


        if not data.empty:

            data = utils.quality_filter(data,site)

            data['Uh'] = utils.calculate_Uh(data['WS'],data['ZL'],data['USTAR'],z['WS'],data['canz'],data['d'])
            data['Uh'] = data['Uh'].interpolate(method='linear',limit=5,limit_direction='both')

            TAcols = list(data.columns[data.columns.str.startswith('TA_')])
            TPOT = utils.convert_to_potential_temperature(data.loc[:,TAcols],z[TAcols])
            for col in TPOT.columns:
                data[col.replace('TA','TPOT')] = TPOT[col]
                z[col.replace('TA','TPOT')] = z[col]

            TPOTcols = list(data.columns[data.columns.str.startswith('TPOT_')])
            data.loc[:,TPOTcols],TPOTbias = utils.fix_temperature_bias(data.loc[:,TPOTcols],data['H'],data['USTAR'],z[TPOTcols])
            data[TPOTcols] = data[TPOTcols].interpolate(method='linear',limit=5,limit_direction='both')

            omega,wec,reftheta = utils.calculate_omega(data['W_SIGMA'],data.loc[:,TPOTcols]+c.NT,data['Uh'],data['LAI'],data['canz'].iloc[0],z['FC'],z[TPOTcols])
            data['OMEGA'] = omega
            data['WEC'] = wec

            data['FC_REF'] = utils.calculate_reference_flux(data['FC'],data['SW_IN'],data['USTAR'])

            data['z'] = z['FC']
            TPOTz = ['z_' + s for s in TPOTcols]
            for indx in range(len(TPOTcols)):data[TPOTz[indx]] = z[TPOTcols[indx]]



            data['FC_GF_OMEGA'] = data['FC'].copy()
            data['FC_GF_USTAR'] = data['FC'].copy()
            data['NEE_GF_OMEGA'] = data['NEE'].copy()
            data['NEE_GF_USTAR'] = data['NEE'].copy()
            data['OMEGA_FLAG'] = 0
            data['FC_USTAR_FLAG'] = 0
            data['NEE_USTAR_FLAG'] = 0
            data.loc[(pd.isnull(data['USTAR'])) | (pd.isnull(data['OMEGA'])),'FC_GF_USTAR'] = np.nan
            data.loc[(pd.isnull(data['USTAR'])) | (pd.isnull(data['OMEGA'])),'NEE_GF_USTAR'] = np.nan
            data.loc[(pd.isnull(data['USTAR'])) | (pd.isnull(data['OMEGA'])),'FC_GF_OMEGA'] = np.nan
            data.loc[(pd.isnull(data['USTAR'])) | (pd.isnull(data['OMEGA'])),'NEE_GF_OMEGA'] = np.nan
            data.loc[(pd.isnull(data['USTAR'])) | (pd.isnull(data['OMEGA'])),'OMEGA_FLAG'] = np.nan
            data.loc[(pd.isnull(data['USTAR'])) | (pd.isnull(data['OMEGA'])),'FC_USTAR_FLAG'] = np.nan
            data.loc[(pd.isnull(data['USTAR'])) | (pd.isnull(data['OMEGA'])),'NEE_USTAR_FLAG'] = np.nan

            for flx in ['FC','NEE']:
                
                USTARth = iotools.load_ONEflux_threshold(os.path.join(os.getcwd(),'data'),site,flx)
                if not USTARth.empty:
                    USTth_average = USTARth.iloc[4].mean()
                else:
                    USTth_average = np.nan
                if pd.isnull(USTth_average):
                    print(['No valid USTAR threshold estimated by ONEflux for '+site+' from '+flx])
                    data[flx+'_GF_USTAR'] = np.nan
                for yr in data.index.year.unique():
                    if yr in USTARth.columns:
                        USTth = USTARth[yr].iloc[4]
                    else:
                        USTth = np.nan

                    if pd.isnull(USTth):
                        USTth = USTth_average

                    data.loc[(data.index.year==yr) & (data['USTAR']<USTth),flx+'_GF_USTAR'] = np.nan
                    data.loc[(data.index.year==yr) & (data['USTAR']<USTth),flx+'_USTAR_FLAG'] = 1


            data.loc[data['OMEGA']<OMEGAth,'FC_GF_OMEGA'] = np.nan
            data.loc[data['OMEGA']<OMEGAth,'NEE_GF_OMEGA'] = np.nan
            data.loc[data['OMEGA']<OMEGAth,'OMEGA_FLAG'] = 1

            outdir = os.path.join(os.getcwd(),'data')
            if not os.path.isdir(outdir):
                os.mkdir(outdir)
            fileout = os.path.join(outdir,site+'.csv')
            iotools.output_processed_data(data,fileout)

            # TPOTbias.to_csv(os.path.join(outdir,'TPOT_bias_'+site+'.csv'),header=False)



if __name__ == "__main__":
    main()