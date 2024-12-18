import numpy as np
import pandas as pd
import datetime as dt
import os
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess
from cmcrameri import cm as cmp
import matplotlib as mpl
from scipy import optimize
from scipy.stats import gaussian_kde
import scipy
from dotenv import load_dotenv;load_dotenv()

import iotools
import utils
import constants as c

sns.set_theme(style='ticks',context='paper')
cm = 1/2.54

var='a'
al = []
al = [(chr(ord(var)+i)) for i in range(26)]

# OMEGAth = 0.66
# OMEGAth_std = 0.06
# OMEGAth_skewness = 0.52
# OMEGAth_quantiles = [0.59,0.61,0.63,0.64,0.65,0.67,0.69,0.71,0.73]

OMEGAth = 0.59
OMEGAth_iqr = [0.50,0.63]
OMEGAth_std = 0.10
OMEGAth_skewness = 0.35
OMEGAth_quantiles = [0.44,0.49,0.51,0.56,0.59,0.60,0.62,0.64,0.68]
config_file = 'config.yml'
config = iotools.load_config(config_file)   

figfold = os.path.join(os.getcwd(),'figs')
if not os.path.isdir(figfold):
    try:
        os.mkdir(figfold)
    except:
        print('Cannot create directory '+figfold)


# figsize = [18,12]
figsize = [18,8]
figsize2 = [12,8]
figsize3 = [18,8]


def fig3_4_5(data,clrs,sites=None):


    p0 = [0.6,0.02,1]
    bounds = ([0.02,-0.1,0.5],[1.5,1,np.Inf])

    if sites is None:
        sites = list(data.keys())
    bins = dict()
    bins['USTAR'] = np.linspace(0,1,20)
    bins['W_SIGMA'] = np.linspace(0,1,20)
    bins['OMEGA'] = np.linspace(0,2,20)
    xlms = {'USTAR':[0,1],'W_SIGMA':[0,1],'OMEGA':[0,1.7]}
    xlabels = {'USTAR':'u$_{*}$ (m/s)','W_SIGMA':'$\sigma_w$ (m/s)','OMEGA':'$\Omega$ (-)'}
    xticks = {'USTAR':[0,0.15,0.3,0.45,0.6,0.75,0.9],'W_SIGMA':[0,0.15,0.3,0.45,0.6,0.75,0.9],'OMEGA':[0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6]}
    Nmin = 30

    # fig, axs = plt.subplot_mosaic([['USTAR FC','W_SIGMA FC','OMEGA FC']],gridspec_kw={'bottom':0.17,'top':0.97,'right':0.97,'wspace':0.05,'hspace':0.05})
    fig, axs = plt.subplot_mosaic([['USTAR FC','W_SIGMA FC','OMEGA FC']],gridspec_kw={'bottom':0.25,'right':0.97,'wspace':0.05})
    # fig, ax = plt.subplot_mosaic([['A','B'],['cb','cb']],gridspec_kw={'wspace':0.05,'hspace':0.35,"height_ratios":[1, 0.05]})
    fig.set_size_inches(figsize[0]*cm, figsize[1]*cm,forward=True)
    for ax in axs:
        axs[ax].fill_between([0,10],[0.9,0.9],[1.1,1.1], linewidth=0,color=[0.8,0.8,0.8])
        axs[ax].plot([0,10],[1,1],'k--')
    
    # fig2, axs2 = plt.subplot_mosaic([['USTAR NEE','W_SIGMA NEE','OMEGA NEE']],gridspec_kw={'bottom':0.17,'top':0.97,'right':0.97,'wspace':0.05,'hspace':0.05})
    fig2, axs2 = plt.subplot_mosaic([['USTAR NEE','W_SIGMA NEE','OMEGA NEE']],gridspec_kw={'bottom':0.25,'right':0.97,'wspace':0.05})
    # fig, ax = plt.subplot_mosaic([['A','B'],['cb','cb']],gridspec_kw={'wspace':0.05,'hspace':0.35,"height_ratios":[1, 0.05]})
    fig2.set_size_inches(figsize[0]*cm, figsize[1]*cm,forward=True)
    for ax in axs2:
        axs2[ax].fill_between([0,10],[0.9,0.9],[1.1,1.1], linewidth=0,color=[0.8,0.8,0.8])
        axs2[ax].plot([0,10],[1,1],'k--')
        # if 'OMEGA' in ax:
        #     axs[ax].plot([OMEGAth,OMEGAth],[-4,5],'k--')



    # FCbinned = pd.DataFrame(index=np.mean(bins['OMEGA'].reshape(-1, 2), axis=1),columns=sites)
            
    FCbinned = dict()
    NEEbinned = dict()
    for key in bins.keys():FCbinned[key] = pd.DataFrame(index=np.linspace(0,len(bins['OMEGA'])-2,len(bins['OMEGA'])-1),columns=sites)
    for key in bins.keys():NEEbinned[key] = pd.DataFrame(index=np.linspace(0,len(bins['OMEGA'])-2,len(bins['OMEGA'])-1),columns=sites)
    for site in sites:
    # for site in ['FR-Aur']:
        # data[site]['FC_REF'] = utils.calculate_reference_flux(data[site]['FC'],data[site]['SW_IN'],data[site]['USTAR'],SWINlm=10,window=dt.timedelta(days=20),min_periods=10)
        # fig, axs = plt.subplot_mosaic([['USTAR FC','W_SIGMA FC','OMEGA FC'],['USTAR NEE','W_SIGMA NEE','OMEGA NEE']],gridspec_kw={'wspace':0.05,'hspace':0.05})
        # # fig, ax = plt.subplot_mosaic([['A','B'],['cb','cb']],gridspec_kw={'wspace':0.05,'hspace':0.35,"height_ratios":[1, 0.05]})
        # fig.set_size_inches(21*cm, 12*cm,forward=True)
        # for ax in axs:
        #     axs[ax].fill_between([0,10],[0.9,0.9],[1.1,1.1], linewidth=0,color=[0.8,0.8,0.8])
        #     axs[ax].plot([0,10],[1,1],'k--')
        #     if 'OMEGA' in ax:
        #         axs[ax].plot([OMEGAth,OMEGAth],[-4,5],'k--')

        TA = data[site].loc[:,data[site].columns.str.startswith('TA_')].median(axis=1)
        TA = TA.rolling(window = dt.timedelta(days=2), center = True).quantile(.5)
        FCbins = dict()
        NEEbins = dict()
        for key in bins.keys():
            FCbins[key] = pd.DataFrame(index=range(len(bins[key])-1),columns=[key,'FC'])
            for i in range(len(bins[key])-1):

                indx = (data[site][key]>bins[key][i]) & (data[site][key]<=bins[key][i+1]) & (data[site]['WEC']!=0) & (data[site]['SW_IN']<10) & (TA>5) & (~data[site]['FC'].isna()) & (~data[site]['FC_REF'].isna())
                if indx.sum()>Nmin:
                    FCbins[key].loc[i,key] = data[site].loc[indx,key].median()
                    FCbins[key].loc[i,'FC'] = (data[site].loc[indx,'FC']/data[site].loc[indx,'FC_REF']).median()
                    # FCbins[key].loc[i,'FC'] = data[site].loc[indx,'FC'].median()/data[site].loc[indx,'FC_REF'].median()
                    
                    FCbinned[key].loc[i,site] = FCbins[key].loc[i,'FC']



            axlbl = key+' FC'
            if site=='WOOD':
                axs[axlbl].plot(FCbins[key].loc[:,key],FCbins[key].loc[:,'FC'],color=clrs[site],label='z < 10 m')
            elif site=='FI-Hyy':
                axs[axlbl].plot(FCbins[key].loc[:,key],FCbins[key].loc[:,'FC'],color=clrs[site],label='10 m < z < 30 m')
            elif site=='DE-Hai':
                axs[axlbl].plot(FCbins[key].loc[:,key],FCbins[key].loc[:,'FC'],color=clrs[site],label='z > 30 m')
            else:
                axs[axlbl].plot(FCbins[key].loc[:,key],FCbins[key].loc[:,'FC'],color=clrs[site])

            axs[axlbl].set_ylim([-0.05,1.5])
            if key!='USTAR':
                axs[axlbl].set_yticklabels([])
            else:
                axs[axlbl].set_ylabel('Normalised F$_c$ (-)')
            axs[axlbl].set_xticks(xticks[key])
            axs[axlbl].set_xlim(xlms[key])
            axs[axlbl].set_xlabel(xlabels[key])
            # axs[axlbl].set_title(site)

            
            NEEbins[key] = pd.DataFrame(index=range(len(bins[key])-1),columns=[key,'NEE'])
            for i in range(len(bins[key])-1):

                indx = (data[site][key]>bins[key][i]) & (data[site][key]<bins[key][i+1]) & (data[site]['WEC']!=0) & (data[site]['SW_IN']<10) & (TA>5) & (~data[site]['NEE'].isna()) & (~data[site]['FC_REF'].isna())
                if indx.sum()>Nmin:
                    NEEbins[key].loc[i,key] = data[site].loc[indx,key].median()
                    NEEbins[key].loc[i,'NEE'] = (data[site].loc[indx,'NEE']/data[site].loc[indx,'FC_REF']).median()
                    # FCbins[key].loc[i,'FC'] = data[site].loc[indx,'FC'].median()/data[site].loc[indx,'FC_REF'].median()
                    
                    NEEbinned[key].loc[i,site] = NEEbins[key].loc[i,'NEE']
            axlbl = key+' NEE'
            if site=='WOOD':
                axs2[axlbl].plot(NEEbins[key].loc[:,key],NEEbins[key].loc[:,'NEE'],color=clrs[site],label='z < 10 m')
            elif site=='FI-Hyy':
                axs2[axlbl].plot(NEEbins[key].loc[:,key],NEEbins[key].loc[:,'NEE'],color=clrs[site],label='10 m < z < 30 m')
            elif site=='DE-Hai':
                axs2[axlbl].plot(NEEbins[key].loc[:,key],NEEbins[key].loc[:,'NEE'],color=clrs[site],label='z > 30 m')
            else:
                axs2[axlbl].plot(NEEbins[key].loc[:,key],NEEbins[key].loc[:,'NEE'],color=clrs[site])
            # axs2[axlbl].plot(NEEbins[key].loc[:,key],NEEbins[key].loc[:,'NEE'],color=clrs[site])
            axs2[axlbl].set_ylim([-0.05,1.5])
            if key!='USTAR':
                axs2[axlbl].set_yticklabels([])
            else:
                axs2[axlbl].set_ylabel('Normalised F$_c$+S$_c$ (-)')
            axs2[axlbl].set_xticks(xticks[key])
            axs2[axlbl].set_xlim(xlms[key])
            axs2[axlbl].set_xlabel(xlabels[key])
            

    xvar = 'OMEGA'
    omth1 = pd.DataFrame(index=FCbinned.keys(),columns=['val','unc'])
    coefs = pd.DataFrame(index=FCbinned.keys(),columns=[0,1,2])
    coefs2 = pd.DataFrame(index=FCbinned.keys(),columns=[0,1,2])
    for key in FCbinned[xvar].keys():
        # print(key)
        indx = FCbinned[xvar][key].index[~FCbinned[xvar][key].isnull()].astype(int).to_numpy()
        x = bins[xvar][0:-1]+(bins[xvar][1]-bins[xvar][0])/2
        if len(indx)>4:
            # coef,e = optimize.curve_fit(utils.piecewise_linear, x[indx],FCbinned[xvar][key].to_numpy()[indx],p0=[0.6,1,1,0],bounds=([0.2,0.7,0,0],[0.8,1.3,30,10]))
            coef,e = optimize.curve_fit(utils.piecewise_linear, x[indx],FCbinned[xvar][key].to_numpy()[indx],p0=p0,bounds=bounds,nan_policy='omit')
            # coef,e = optimize.curve_fit(utils.piecewise_linear, x[indx],FCbinned[xvar][key].to_numpy()[indx],p0=[0.1,1,1,0],bounds=([0.03,0.7,0,0],[0.8,1.3,30,10]))
            omth1.loc[key,'val'] = coef[0]
            omth1.loc[key,'unc'] = e[0,0]
            coefs.loc[key,0] = coef[0]
            coefs.loc[key,1] = coef[1]
            coefs.loc[key,2] = coef[2]
        else:
            omth1.loc[key,'val'] = np.nan
            omth1.loc[key,'unc'] = np.nan
    
    omth2 = pd.DataFrame(index=NEEbinned.keys(),columns=['val','unc'])
    for key in NEEbinned[xvar].keys():
        # print(key)
        indx = NEEbinned[xvar][key].index[~NEEbinned[xvar][key].isnull()].astype(int).to_numpy()
        x = bins[xvar][0:-1]+(bins[xvar][1]-bins[xvar][0])/2
        if len(indx)>4:
            # coef,e = optimize.curve_fit(utils.piecewise_linear, x[indx],FCbinned[xvar][key].to_numpy()[indx],p0=[0.6,1,1,0],bounds=([0.2,0.7,0,0],[0.8,1.3,30,10]))
            coef,e = optimize.curve_fit(utils.piecewise_linear, x[indx],NEEbinned[xvar][key].to_numpy()[indx],p0=p0,bounds=bounds,nan_policy='omit')
            # coef,e = optimize.curve_fit(utils.piecewise_linear, x[indx],FCbinned[xvar][key].to_numpy()[indx],p0=[0.1,1,1,0],bounds=([0.03,0.7,0,0],[0.8,1.3,30,10]))
            omth2.loc[key,'val'] = coef[0]
            omth2.loc[key,'unc'] = e[0,0]
            coefs2.loc[key,0] = coef[0]
            coefs2.loc[key,1] = coef[1]
            coefs2.loc[key,2] = coef[2]
        else:
            omth2.loc[key,'val'] = np.nan
            omth2.loc[key,'unc'] = np.nan
        # fig3, axs3 = plt.subplot_mosaic([[key]])
        # axs3[key].plot(x,NEEbinned[xvar][key].to_numpy(),'b.')
        # axs3[key].plot(x,FCbinned[xvar][key].to_numpy(),'r.')
        # if key in coefs2.index:
        #     axs3[key].plot(np.linspace(0,2,100),utils.piecewise_linear(np.linspace(0,2,100), coefs2.loc[key,0], coefs2.loc[key,1], coefs2.loc[key,2]),'b-')
        #     axs3[key].plot([coefs2.loc[key,0],coefs2.loc[key,0]],[0,1.4],'b:')
        # if key in coefs.index:
        #     axs3[key].plot(np.linspace(0,2,100),utils.piecewise_linear(np.linspace(0,2,100), coefs.loc[key,0], coefs.loc[key,1], coefs.loc[key,2]),'r-')
        #     axs3[key].plot([coefs.loc[key,0],coefs.loc[key,0]],[0,1.4],'r:')
        # axs3[key].plot([OMEGAth,OMEGAth],[0,1.4],'k:')
        # axs3[key].set_title(key+' '+str(omth1.loc[key,'val'])+' '+str(omth2.loc[key,'val']))
        # fig3.savefig(figfold+'/' + 'fig3_' + key + '.png', dpi=200,bbox_inches='tight')

    # omth_for_plot = omth1.loc[~omth1['val'].isnull(),'val']
    # axs['OMEGA FC'].violinplot([omth_for_plot.to_list()],showmeans=False,showmedians=True,vert=False,positions=[0.2],widths=[0.2],alpha=1)

    # xvar = 'OMEGA'
    # x = bins[xvar][0:-1]+(bins[xvar][1]-bins[xvar][0])/2
    # y = FCbinned[xvar].mean(axis=1).astype(float).to_numpy()
    # yerr = FCbinned[xvar].std(axis=1).astype(float).to_numpy()
    # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,sigma=yerr,p0=[0.6,1,1],bounds=([0.2,0,0,0],[0.8,2,30]))
    # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,sigma=yerr,p0=[0.6,1,1],bounds=([0.02,0,0.9],[1,10,1.1]),nan_policy='omit')
    # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y-yerr,sigma=yerr)
    # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y+yerr,sigma=yerr)
    # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,p0=[0.6,1,1],bounds=([0.2,0.7,0],[0.8,1.3,30]))
    # omth.loc[key,'val'] = coef[0]
    # omth.loc[key,'unc'] = e[0,0]**0.5

    
    xvar = 'OMEGA'
    n = 1000
    omth = pd.DataFrame(index=np.linspace(0,n-1,n),columns=['val','unc'])
    x = bins[xvar][0:-1]+(bins[xvar][1]-bins[xvar][0])/2
    coefs3 = pd.DataFrame(index=np.linspace(0,n-1,n),columns=[0,1,2])
    for i in range(n):
        # print(i)
        indx_for_y = np.random.randint(0, high=len(FCbinned[xvar].columns)-1, size=(len(FCbinned[xvar].index),1), dtype=int)
        y = pd.Series(index=FCbinned[xvar].index)
        for j in range(len(y)):
            y[j] = FCbinned[xvar].iloc[j,indx_for_y[j]].iloc[0]
            # y[j] = cols[indx_for_y[j]]
        if (~y.isnull()).sum()>5:
            y = y.astype(float).to_numpy()
            # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,p0=[0.6,1,1,0],bounds=([0.1,0.5,0,0],[0.9,1.5,30,10]),nan_policy='omit')
            # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,p0=[0.6,1,1,0],bounds=([0.02,0,0,0],[1,10,30,10]),nan_policy='omit')
            coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,p0=p0,bounds=bounds,nan_policy='omit')
            coefs3.loc[i,0] = coef[0]
            coefs3.loc[i,1] = coef[1]
            coefs3.loc[i,2] = coef[2]
            # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,nan_policy='omit')
            omth.loc[i,'val'] = coef[0]
            omth.loc[i,'unc'] = e[0,0]**0.5
            # plt.plot(x,y,'o-')
            # plt.plot(x,utils.piecewise_linear(x, coef[0], coef[1], coef[2]),'-')
    parts = axs['OMEGA FC'].violinplot([omth['val'].to_list()],showmeans=False,showmedians=False,vert=False,showextrema=False,positions=[0.15],widths=[0.2])
    quantiles = omth['val'].quantile([0.25,0.75])
    for pc in parts['bodies']:
        # pc.set_facecolor('#D43F3A')
        pc.set_facecolor((0.9,0.9,0.9,0))
        pc.set_edgecolor('black')
        pc.set_alpha(1)
        pc.set_zorder(2)

        
    for ax in axs:
        if 'OMEGA' in ax:
            axs[ax].plot([omth['val'].median(),omth['val'].median()],[-4,5],'k--',zorder=1)
        
    axs['OMEGA FC'].scatter(omth['val'].median(),0.15, marker='o', color='k', s=10, zorder=3)
    # axs['OMEGA FC'].hlines(0.2, quantiles[0.25], quantiles[0.75], color='k', linestyle='-', lw=5)
    axs['OMEGA FC'].hlines(0.15, omth['val'].min(), omth['val'].max(), color='k', linestyle='-', lw=1, zorder=3)
    
    # xvar = 'OMEGA'
    # n = 500
    # omth = pd.DataFrame(index=np.linspace(0,n-1,n),columns=['val','unc'])
    # x = bins[xvar][0:-1]+(bins[xvar][1]-bins[xvar][0])/2
    # for i in range(n):
    #     # print(i)
    #     indx_for_y = np.random.randint(0, high=len(FCbinned[xvar].columns)-1, size=(len(FCbinned[xvar].index),1), dtype=int)
    #     y = pd.Series(index=FCbinned[xvar].index)
    #     for j in range(len(y)):
    #         y[j] = FCbinned[xvar].iloc[j,indx_for_y[j]].iloc[0]
    #         # y[j] = cols[indx_for_y[j]]
    #     if (~y.isnull()).sum()>5:
    #         y = y.astype(float).to_numpy()
    #         # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,p0=[0.6,1,1,0],bounds=([0.1,0.5,0,0],[0.9,1.5,30,10]),nan_policy='omit')
    #         # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,p0=[0.6,1,1,0],bounds=([0.02,0,0,0],[1,10,30,10]),nan_policy='omit')
    #         coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,p0=[0.6,1,1,0],bounds=([0.02,0,0.8,0],[1,10,1.2,1e-8]),nan_policy='omit')
    #         # coef,e = optimize.curve_fit(utils.piecewise_linear, x,y,nan_policy='omit')
    #         omth.loc[i,'val'] = coef[0]
    #         omth.loc[i,'unc'] = e[0,0]**0.5
    #         # plt.plot(x,y,'o-')
    #         plt.plot(x,utils.piecewise_linear(x, coef[0], coef[1], coef[2], coef[3]),'-')

    # axs[key+' FC'].plot(bins[key][0:-1]+(bins[key][1]-bins[key][0])/2,utils.piecewise_linear(bins[key][0:-1]+(bins[key][1]-bins[key][0])/2, coef[0], coef[1], coef[2], coef[3]),'k-')
    # plt.plot(x,utils.piecewise_linear(x, coef[0], coef[1], coef[2], coef[3]),'k-')

    fig.suptitle('Turbulence-physics point of view')
    fig2.suptitle('Biological point of view')

    
    fig.text(.95, .05, 'a)', ha='right', va='bottom', transform=axs['USTAR FC'].transAxes)
    fig.text(.95, .05, 'b)', ha='right', va='bottom', transform=axs['W_SIGMA FC'].transAxes)
    fig.text(.95, .05, 'c)', ha='right', va='bottom', transform=axs['OMEGA FC'].transAxes)
    
    fig2.text(.95, .05, 'a)', ha='right', va='bottom', transform=axs2['USTAR NEE'].transAxes)
    fig2.text(.95, .05, 'b)', ha='right', va='bottom', transform=axs2['W_SIGMA NEE'].transAxes)
    fig2.text(.95, .05, 'c)', ha='right', va='bottom', transform=axs2['OMEGA NEE'].transAxes)

    axlbl = 'USTAR FC'
    handles, labels = axs[axlbl].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center',ncols=3)
    axlbl = 'USTAR NEE'
    handles, labels = axs2[axlbl].get_legend_handles_labels()
    fig2.legend(handles, labels, loc='lower center',ncols=3)

    for key in bins.keys():
        axs[key+' FC'].errorbar(bins[key][0:-1]+(bins[key][1]-bins[key][0])/2,FCbinned[key].mean(axis=1),color='k',yerr=FCbinned[key].std(axis=1),marker='o',linestyle='none',markerfacecolor='w')
        
    # plt.show()

    omth1["z"] = np.nan
    for site in omth1.index:
        if site in data.keys():
            # print(site)
            omth1.loc[site,'z'] = data[site]["z"].mean()
            omth1.loc[site,'LAI'] = data[site]["LAI"].max()
    omth2["z"] = np.nan
    for site in omth2.index:
        if site in data.keys():
            # print(site)
            omth2.loc[site,'z'] = data[site]["z"].mean()
            omth2.loc[site,'LAI'] = data[site]["LAI"].max()

    # omth1.to_csv(os.getcwd()+os.path.sep+'omegath.csv')
    # omth2.to_csv(os.getcwd()+os.path.sep+'omegath2.csv')
    
            
    # fig3, ax3 = plt.subplot_mosaic([["A"]])
    
    fig3, ax3 = plt.subplot_mosaic([['A','cb','cb2']],gridspec_kw={'wspace':0.05,"width_ratios":[1, 0.02, 0.02]})
    pos = ax3["cb"].get_position()
    pos2 = ax3["cb2"].get_position()
    pos2.x0 = pos.x1
    pos2.x1 = pos2.x0+pos.x1-pos.x0
    ax3["cb2"].set_position(pos2)
    fig3.set_size_inches(figsize2[0]*cm, figsize2[1]*cm,forward=True)
    # ax3["A"].fill_between([-10,omth1["z"].max()+10],[OMEGAth-OMEGAth_std,OMEGAth-OMEGAth_std],[OMEGAth+OMEGAth_std,OMEGAth+OMEGAth_std], linewidth=0,color=[0.8,0.8,0.8])
    ax3["A"].fill_between([-10,omth1["z"].max()+10],[OMEGAth_iqr[0],OMEGAth_iqr[0]],[OMEGAth_iqr[1],OMEGAth_iqr[1]], linewidth=0,color=[0.8,0.8,0.8])
    ax3["A"].plot([0,omth1["z"].max()+10],[OMEGAth,OMEGAth],'k--',label="$\Omega_{th}$="+f"{OMEGAth:.2f}")
    # ax3["A"].errorbar(omth2["z"],omth2["val"],color="r",yerr=omth2["unc"],marker="^",linestyle="none",markerfacecolor="w",label="$F_c+S_c$")
    # ax3["A"].errorbar(omth1["z"],omth1["val"],color="k",yerr=omth1["unc"],marker="o",linestyle="none",markerfacecolor="w",label="$F_c$")
    cmap = mpl.colormaps['Grays']
    edgecolor = cmap(160)
    sc1 = ax3["A"].scatter(omth2["z"],omth2["val"],s=12,c=omth2["LAI"],marker="^",edgecolor=edgecolor,label="$F_c+S_c$",cmap="Grays",vmin=0,vmax=10,alpha=0.7)
    cmap = mpl.colormaps['Blues']
    edgecolor = cmap(250)
    sc2 = ax3["A"].scatter(omth1["z"],omth1["val"],c=omth1["LAI"],marker="o",edgecolor=edgecolor,label="$F_c$",cmap="Blues",vmin=0,vmax=10)
    cb = fig3.colorbar(sc1, cax=ax3["cb"])
    cb2 = fig3.colorbar(sc2, cax=ax3["cb2"])
    cb.set_ticks([])
    cb2.set_label('PAI maximum (m$^2$ m$^{-2}$)')
    ax3["A"].set_xlim([0,omth1["z"].max()+2])
    ax3["A"].set_ylim([0,np.max([omth1["val"].max(),omth2["val"].max()])+0.05])
    ax3["A"].set_ylabel("$\Omega$ threshold (-)")
    ax3["A"].set_xlabel("EC measurement height (m)")
    ax3["A"].legend(loc="upper right")
    # plt.show()
    
    fig.savefig(figfold+'/fig3.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig3.pdf', dpi=200,bbox_inches='tight')
    fig2.savefig(figfold+'/fig4.png', dpi=200,bbox_inches='tight')
    fig2.savefig(figfold+'/fig4.pdf', dpi=200,bbox_inches='tight')
    
    fig3.savefig(figfold+'/fig5.png', dpi=200,bbox_inches='tight')
    fig3.savefig(figfold+'/fig5.pdf', dpi=200,bbox_inches='tight')



def fig7(data,USTARth,clrs):

    fig, axs = plt.subplot_mosaic([['MCC','COMP']],gridspec_kw={'wspace':0.2,'hspace':0.05,'right':0.95,'left':0.1,'top':0.95})
    fig.set_size_inches(figsize3[0]*cm, figsize3[1]*cm,forward=True)
    

    USTAR_limits = pd.DataFrame(index=data.keys(),columns=['MCC','ONEflux','MCC2'])
    for site in data.keys():
        MCC = pd.Series(index=np.linspace(0.01,1,100))
        flags = pd.DataFrame(index=data[site].index,columns=['USTAR_flag','OMEGA_flag'])
        flags.loc[(data[site]['OMEGA']<OMEGAth),'OMEGA_flag'] = 0
        flags.loc[(data[site]['OMEGA']>=OMEGAth),'OMEGA_flag'] = 1
        flags.loc[data[site]['OMEGA'].isna(),'OMEGA_flag'] = np.nan
        flags.loc[data[site]['USTAR'].isna(),'OMEGA_flag'] = np.nan
        flags.loc[data[site]['SW_IN']>10,'OMEGA_flag'] = np.nan

        for th in MCC.index:
            flags['USTAR_flag'] = np.nan
            flags.loc[(data[site]['USTAR']<th),'USTAR_flag'] = 0
            flags.loc[(data[site]['USTAR']>=th),'USTAR_flag'] = 1
            flags.loc[data[site]['OMEGA'].isna(),'USTAR_flag'] = np.nan
            flags.loc[data[site]['USTAR'].isna(),'USTAR_flag'] = np.nan
            flags.loc[data[site]['SW_IN']>10,'USTAR_flag'] = np.nan

            MCC[th] = utils.calculate_MCC(flags['OMEGA_flag'],flags['USTAR_flag'])
        if (~data[site]['OMEGA'].isna()).sum()==0:
            print('Error in omega calculation at ' + site + '. Not calculating MCC.')
        else:

            if site=='WOOD':
                axs['MCC'].plot(MCC,color=clrs[site],label='z < 10 m')
            elif site=='FI-Hyy':
                axs['MCC'].plot(MCC,color=clrs[site],label='10 m < z < 30 m')
            elif site=='DE-Hai':
                axs['MCC'].plot(MCC,color=clrs[site],label='z > 30 m')
            else:
                axs['MCC'].plot(MCC,color=clrs[site])
            axs['MCC'].plot(MCC.index[MCC==MCC.max()][0],MCC.max(),'ko',markerfacecolor=clrs[site])

            USTAR_limits.loc[site,'MCC'] = MCC.index[MCC==MCC.max()][0]
            USTAR_limits.loc[site,'MCC_max'] = MCC.max()
            USTAR_limits.loc[site,'ONEflux'] = USTARth['FC'][site].iloc[4].mean()
            # USTAR_limits.loc[site,'MCC2'] = data[site].loc[(data[site]['OMEGA']>OMEGAth-0.02) & (data[site]['OMEGA']<OMEGAth+0.02) & (data[site]['SW_IN']<10),'USTAR'].mean()
            USTAR_limits.loc[site,'ONEflux'] = USTARth['FC'][site].iloc[4].mean()
            lerr = USTARth['FC'][site].iloc[4].mean()-USTARth['FC'][site].iloc[1].mean()
            uerr = USTARth['FC'][site].iloc[8].mean()-USTARth['FC'][site].iloc[4].mean()
            errlm = np.ones((2,1))*np.nan
            errlm[0,0] = lerr
            errlm[1,0] = uerr
            axs['COMP'].errorbar(USTAR_limits.loc[site,'MCC'],USTAR_limits.loc[site,'ONEflux'],yerr=errlm,fmt='o',markerfacecolor=clrs[site],markeredgecolor='k',color=clrs[site])

            
            if (not USTARth['NEE'][site].empty):
                USTAR_limits.loc[site,'ONEflux_NEE'] = USTARth['NEE'][site].iloc[4].mean()
                # lerr = USTARth['NEE'][site].iloc[4].mean()-USTARth['NEE'][site].iloc[1].mean()
                # uerr = USTARth['NEE'][site].iloc[8].mean()-USTARth['NEE'][site].iloc[4].mean()
                # errlm = np.ones((2,1))*np.nan
                # errlm[0,0] = lerr
                # errlm[1,0] = uerr
                # axs['COMP'].errorbar(USTAR_limits.loc[site,'MCC'],USTAR_limits.loc[site,'ONEflux_NEE'],yerr=errlm,fmt='^',markerfacecolor=clrs[site],markeredgecolor='k',color=clrs[site])

    corr = USTAR_limits.corr()
    
    plt.text(.05, .95, 'a)', ha='left', va='top', transform=axs['MCC'].transAxes)
    plt.text(.05, .95, f"b) r = {corr.loc['MCC','ONEflux']:.2f}", ha='left', va='top', transform=axs['COMP'].transAxes)
    # axs['COMP'].text(0.05,0.65,f"r = {corr.loc['MCC','ONEflux']:.2f}")
    
    axs['COMP'].plot([0,1],[0,1],'k--')
    axs['COMP'].set_ylim([0,0.9])
    axs['COMP'].set_xlim([0,0.9])
    axs['COMP'].set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    axs['COMP'].set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    axs['COMP'].set_xlabel('u$_{*th}$, maximum\nagreement with $\Omega$ (m/s)')
    axs['COMP'].set_ylabel('u$_{*th}$, ONEflux (m/s)')

    axs['MCC'].set_ylim([0,1])
    axs['MCC'].set_xlim([0.05,0.9])
    axs['MCC'].set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    axs['MCC'].set_xlabel('u$_*$ threshold (m/s)')
    axs['MCC'].set_ylabel('MCC (-)')
    axs['MCC'].legend(loc="upper right")


    # plt.show()
    fig.savefig(figfold+'/fig7.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig7.pdf', dpi=200,bbox_inches='tight')

def fig2(data,sID,USTARth):

    xvar = 'USTAR'
    fig, ax = plt.subplot_mosaic([['A'],['cb']],gridspec_kw={'wspace':0.05,'hspace':0.35,"height_ratios":[1, 0.05]})
    fig.set_size_inches(12*cm, 12*cm,forward=True)

    indx = ~(np.isinf(data[sID]['OMEGA']))

    sc = ax['A'].scatter(data[sID].loc[indx,xvar],data[sID].loc[indx,'OMEGA'],c=data[sID].loc[indx,'WEC'],s=20,cmap=cmp.batlow,alpha=0.1,linewidth=0,vmax=0)
    
    ax['A'].plot([-4,5],[OMEGAth,OMEGAth],'k--')
    if xvar=='USTAR':
        ax['A'].plot([USTARth,USTARth],[-4,5],'k--')
    else:
        ax['A'].plot([data[sID].loc[indx,xvar].median(),data[sID].loc[indx,xvar].median()],[-4,5],'k--')
    ax['A'].set_xlim([0,data[sID].loc[indx,xvar].quantile(q=.99)])
    ax['A'].set_ylim([0,data[sID].loc[indx,'OMEGA'].quantile(q=.99)])
    if xvar=='USTAR':
        ax['A'].set_xlabel('u$_\star$ (m/s)')
    else:
        ax['A'].set_xlabel('$\sigma_w$ (m/s)')
    ax['A'].set_ylabel('$\Omega$ (-)')
    ax['A'].set_title(sID)
    
    plt.text(.01, .95, 'Coupled-Decoupled', ha='left', va='top', transform=ax['A'].transAxes)
    plt.text(.99, .95, 'Coupled-Coupled', ha='right', va='top', transform=ax['A'].transAxes)
    plt.text(.99, .05, 'Decoupled-Coupled', ha='right', va='bottom', transform=ax['A'].transAxes)
    plt.text(.01, .05, 'Decoupled-Decoupled', ha='left', va='bottom', transform=ax['A'].transAxes)
    # plt.text(.01, .95, snames[ID], ha='left', va='top', transform=ax['A'].transAxes)
    # plt.text(.01, .95, snames[ID], ha='left', va='top', transform=ax['A'].transAxes)
    # plt.text(.01, .95, snames[ID], ha='left', va='top', transform=ax['A'].transAxes)
    
    cb = fig.colorbar(sc,cax=ax['cb'], orientation="horizontal")
    cb.set_label('w$_{e,crit}$ (m/s)')
    cb.solids.set(alpha=1)
    # cb.ax.xaxis.set_major_locator(mdates.DayLocator(interval=20))
    # cb.ax.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m'))
    # plt.show()
    fig.savefig(figfold+'/fig2.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig2.pdf', dpi=200,bbox_inches='tight')

def fig9(data,clrs):

    df = pd.DataFrame(index=data.keys(),columns=['OMEGA','WEC','W_SIGMA','LAI','canz','z'])
    for site in data.keys():
        dat = data[site].copy()
        if (dat['canz'].iloc[0]<1):
            # print(site + ' ' + str(dat['LAI'].iloc[0]))
            dat['LAI'] = 0
        if ((dat['canz'].iloc[0]>1) & (dat['LAI'].iloc[0]>0)) | ((dat['canz'].iloc[0]<1)):
        # if ((dat['canz'].iloc[0]>1) & (dat['LAI'].iloc[0]>0)):
            indx = (dat['SW_IN']<10) & (dat['WEC']<0) & (dat['ZL'].abs()<0.05)
            
            # indx = (dat['SW_IN']<10) & (dat['WEC']<0)
            LAI_quantiles = dat.loc[indx,'LAI'].quantile([.25,.75])
            if LAI_quantiles[0.25]!=LAI_quantiles[0.75]:
                # LAI changes in time
                indx = (dat['SW_IN']<10) & (dat['WEC']<0) & (dat['ZL'].abs()<0.05) & (dat['LAI']<=dat['LAI'].quantile(0.25))
                df.loc[site+' low LAI','OMEGA'] = dat.loc[indx,'OMEGA'].mean()
                df.loc[site+' low LAI','WEC'] = dat.loc[indx,'WEC'].mean()
                # df.loc[site,'W_SIGMA'] = (dat.loc[indx,'W_SIGMA']/dat.loc[indx,'Uh']).mean()
                df.loc[site+' low LAI','W_SIGMA'] = (dat.loc[indx,'W_SIGMA']).mean()
                df.loc[site+' low LAI','LAI'] = dat.loc[indx,'LAI'].mean()
                df.loc[site+' low LAI','z'] = dat.loc[indx,'z'].mean()
                df.loc[site+' low LAI','canz'] = dat.loc[indx,'canz'].mean()
                
                indx = (dat['SW_IN']<10) & (dat['WEC']<0) & (dat['ZL'].abs()<0.05) & (dat['LAI']<dat['LAI'].quantile(0.75)) & (dat['LAI']>dat['LAI'].quantile(0.25))
                df.loc[site+' medium LAI','OMEGA'] = dat.loc[indx,'OMEGA'].mean()
                df.loc[site+' medium LAI','WEC'] = dat.loc[indx,'WEC'].mean()
                # df.loc[site,'W_SIGMA'] = (dat.loc[indx,'W_SIGMA']/dat.loc[indx,'Uh']).mean()
                df.loc[site+' medium LAI','W_SIGMA'] = (dat.loc[indx,'W_SIGMA']).mean()
                df.loc[site+' medium LAI','LAI'] = dat.loc[indx,'LAI'].mean()
                df.loc[site+' medium LAI','z'] = dat.loc[indx,'z'].mean()
                df.loc[site+' medium LAI','canz'] = dat.loc[indx,'canz'].mean()
                
                indx = (dat['SW_IN']<10) & (dat['WEC']<0) & (dat['ZL'].abs()<0.05) & (dat['LAI']>=dat['LAI'].quantile(0.75))
                df.loc[site+' high LAI','OMEGA'] = dat.loc[indx,'OMEGA'].mean()
                df.loc[site+' high LAI','WEC'] = dat.loc[indx,'WEC'].mean()
                # df.loc[site,'W_SIGMA'] = (dat.loc[indx,'W_SIGMA']/dat.loc[indx,'Uh']).mean()
                df.loc[site+' high LAI','W_SIGMA'] = (dat.loc[indx,'W_SIGMA']).mean()
                df.loc[site+' high LAI','LAI'] = dat.loc[indx,'LAI'].mean()
                df.loc[site+' high LAI','z'] = dat.loc[indx,'z'].mean()
                df.loc[site+' high LAI','canz'] = dat.loc[indx,'canz'].mean()
            else:
                df.loc[site,'OMEGA'] = dat.loc[indx,'OMEGA'].mean()
                df.loc[site,'WEC'] = dat.loc[indx,'WEC'].mean()
                # df.loc[site,'W_SIGMA'] = (dat.loc[indx,'W_SIGMA']/dat.loc[indx,'Uh']).mean()
                df.loc[site,'W_SIGMA'] = (dat.loc[indx,'W_SIGMA']).mean()
                df.loc[site,'LAI'] = dat.loc[indx,'LAI'].mean()
                df.loc[site,'z'] = dat.loc[indx,'z'].mean()
                df.loc[site,'canz'] = dat.loc[indx,'canz'].mean()
        del dat
    

    fig, ax = plt.subplot_mosaic([['LAI']],gridspec_kw={'wspace':0.2,'hspace':0.05,'right':0.95,'left':0.1,'top':0.95})
    fig.set_size_inches(figsize2[0]*cm, figsize2[1]*cm,forward=True)

    filtered = lowess(df['WEC'], df['LAI'],frac=0.7)
    filtered2 = lowess(df['W_SIGMA'], df['LAI'],frac=0.7)

    # for site in data.keys():
    #     if site+ ' medium LAI' in df.index:
            
    #         ax['LAI'].plot(df.loc[site+' low LAI','LAI'],-df.loc[site+' low LAI','WEC'],'k^',markerfacecolor=clrs[site],label='|w$_{e,crit}$|')
    #         ax['LAI'].plot(df.loc[site+' low LAI','LAI'],df.loc[site+' low LAI','W_SIGMA'],'ko',markerfacecolor=clrs[site],label='$\sigma_w$')
    #         ax['LAI'].plot(df.loc[site+' medium LAI','LAI'],-df.loc[site+' medium LAI','WEC'],'k^',markerfacecolor=clrs[site],label='|w$_{e,crit}$|')
    #         ax['LAI'].plot(df.loc[site+' medium LAI','LAI'],df.loc[site+' medium LAI','W_SIGMA'],'ko',markerfacecolor=clrs[site],label='$\sigma_w$')
    #         ax['LAI'].plot(df.loc[site+' high LAI','LAI'],-df.loc[site+' high LAI','WEC'],'k^',markerfacecolor=clrs[site],label='|w$_{e,crit}$|')
    #         ax['LAI'].plot(df.loc[site+' high LAI','LAI'],df.loc[site+' high LAI','W_SIGMA'],'ko',markerfacecolor=clrs[site],label='$\sigma_w$')
    #     else:
    #         ax['LAI'].plot(df.loc[site,'LAI'],-df.loc[site,'WEC'],'k^',markerfacecolor=clrs[site],label='|w$_{e,crit}$|')
    #         ax['LAI'].plot(df.loc[site,'LAI'],df.loc[site,'W_SIGMA'],'ko',markerfacecolor=clrs[site],label='$\sigma_w$')

    # # omega, assuming constant turbulence intensity (sigma_w/U=0.35)
    # ome = 1/(2*0.277)*0.3*1/(0.2*df['LAI'])

    ax['LAI'].plot(df['LAI'],-df['WEC'],'k^',markerfacecolor='w',label='|w$_{e,crit}$|')
    ax['LAI'].plot(filtered[:,0],-filtered[:,1],'k-')
    ax['LAI'].plot(df['LAI'],df['W_SIGMA'],'ro',markerfacecolor='w',label='$\sigma_w$')
    ax['LAI'].plot(filtered2[:,0],filtered2[:,1],'r-')
    # ax['LAI'].plot(df['LAI'],df['W_SIGMA']/ome,'bo',markerfacecolor='w')
    ax['LAI'].set_xlabel('PAI (m$^2$ m$^{-2}$)')
    ax['LAI'].set_ylabel('|w$_{e,crit}$| or $\sigma_w$ (m s$^{-1}$)')
    ax['LAI'].set_title('Night time near neutral conditions')
    ax['LAI'].legend(loc="upper left")
    ax['LAI'].set_xlim([-0.05,10])
    # plt.show()
    fig.savefig(figfold+'/fig9.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig9.pdf', dpi=200,bbox_inches='tight')

def fig10(data,clrs):

    df = pd.DataFrame(index=data.keys(),columns=['OMEGA','WEC','W_SIGMA','LAI','canz','z'])
    for site in data.keys():
        indx = (data[site]['SW_IN']<10) & (data[site]['WEC']<0)
        # indx = (data[site]['SW_IN']<10) & (data[site]['WEC']<0)
        df.loc[site,'OMEGA'] = data[site].loc[indx,'OMEGA'].mean()
        df.loc[site,'WEC'] = data[site].loc[indx,'WEC'].mean()
        # df.loc[site,'W_SIGMA'] = (data[site].loc[indx,'W_SIGMA']/data[site].loc[indx,'Uh']).mean()
        df.loc[site,'W_SIGMA'] = (data[site].loc[indx,'W_SIGMA']).mean()
        df.loc[site,'LAI'] = data[site].loc[indx,'LAI'].mean()
        df.loc[site,'z'] = data[site].loc[indx,'z'].mean()
        df.loc[site,'canz'] = data[site].loc[indx,'canz'].mean()
    
    filtered = lowess(df['WEC'], df['z'],frac=0.5)
    filtered2 = lowess(df['W_SIGMA'], df['z'],frac=0.5)

    fig, ax = plt.subplot_mosaic([['z']],gridspec_kw={'wspace':0.2,'hspace':0.05,'right':0.95,'left':0.1,'top':0.95})
    fig.set_size_inches(figsize2[0]*cm, figsize2[1]*cm,forward=True)
    # ax['A'].plot(df['LAI'],df['WEC'].abs(),'k.')
    ax['z'].plot(df['z'],-df['WEC'],'k^',markerfacecolor='w',label='|w$_{e,crit}$|')
    ax['z'].plot(filtered[:,0],-filtered[:,1],'k-')
    # ax['z'].plot(filtered[:,0],-filtered[:,1]*OMEGAth,'k--')
    ax['z'].plot(df['z'],df['W_SIGMA'],'ro',markerfacecolor='w',label='$\sigma_w$')
    ax['z'].plot(filtered2[:,0],filtered2[:,1],'r-')
    ax['z'].set_xlabel('EC height (m)')
    ax['z'].set_ylabel('|w$_{e,crit}$| or $\sigma_w$ (m s$^{-1}$)')
    ax['z'].set_title('Night time')
    ax['z'].legend(loc="upper left")

    fig.savefig(figfold+'/fig10.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig10.pdf', dpi=200,bbox_inches='tight')


def fig8(data,clrs):


    bins = np.logspace(-2,0.4,30)
    xlms = [1e-2,2]
    Nmin = 50

    fig, axs = plt.subplot_mosaic([['A']],gridspec_kw={'wspace':0.2,'hspace':0.05,'right':0.95,'left':0.1,'top':0.95})
    fig.set_size_inches(figsize2[0]*cm, figsize2[1]*cm,forward=True)




    for site in data.keys():

        OMEGAbins = pd.DataFrame(columns=['OMEGA','ZL'])
        for i in range(len(bins)-1):
            indx = (data[site]['ZL']>bins[i]) & (data[site]['ZL']<=bins[i+1]) & (data[site]['WEC']!=0) & (data[site]['SW_IN']<10)
            if indx.sum()>Nmin:
                OMEGAbins.loc[i,'OMEGA'] = data[site].loc[indx,'OMEGA'].median()
                OMEGAbins.loc[i,'ZL'] = data[site].loc[indx,'ZL'].median()
                    
        if site=='WOOD':
            axs['A'].plot(OMEGAbins.loc[:,'ZL'],OMEGAbins.loc[:,'OMEGA'],color=clrs[site],label='z < 10 m')
        elif site=='FI-Hyy':
            axs['A'].plot(OMEGAbins.loc[:,'ZL'],OMEGAbins.loc[:,'OMEGA'],color=clrs[site],label='10 m < z < 30 m')
        elif site=='DE-Hai':
            axs['A'].plot(OMEGAbins.loc[:,'ZL'],OMEGAbins.loc[:,'OMEGA'],color=clrs[site],label='z > 30 m')
        else:
            axs['A'].plot(OMEGAbins.loc[:,'ZL'],OMEGAbins.loc[:,'OMEGA'],color=clrs[site])

    axs['A'].plot([1e-5,5],[OMEGAth,OMEGAth],'k--')
    axs['A'].set_ylim([0,3])
    axs['A'].set_xscale('log')
    axs['A'].set_xlim(xlms)
    axs['A'].set_ylabel('$\Omega$ (-)')
    axs['A'].set_xlabel('Stability parameter, $\zeta$ (-)')

    
    # handles, labels = axs['A'].get_legend_handles_labels()
    axs['A'].legend(loc='upper right',ncols=1)

    # plt.show()
    fig.savefig(figfold+'/fig8.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig8.pdf', dpi=200,bbox_inches='tight')



def fig12(data,clrs):

    data_cov = dict()
    data_cov['day'] = pd.DataFrame(index=data.keys(),columns=['OMEGA_filt_FC','USTAR_filt_FC','OMEGA_filt_NEE','USTAR_filt_NEE','LAI','canz','z','WEC'])
    data_cov['night'] = pd.DataFrame(index=data.keys(),columns=['OMEGA_filt_FC','USTAR_filt_FC','OMEGA_filt_NEE','USTAR_filt_NEE','LAI','canz','z','WEC'])
    for site in data.keys():
        # omega_FC_filt = ((data[site]['SW_IN']<10) & (data[site]['OMEGA_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['FC_USTAR_FLAG'])))).sum()
        # ustar_FC_filt = ((data[site]['SW_IN']<10) & (data[site]['FC_USTAR_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['FC_USTAR_FLAG'])))).sum()
        # omega_NEE_filt = ((data[site]['SW_IN']<10) & (data[site]['OMEGA_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['NEE_USTAR_FLAG'])))).sum()
        # ustar_NEE_filt = ((data[site]['SW_IN']<10) & (data[site]['NEE_USTAR_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['NEE_USTAR_FLAG'])))).sum()
        # FC_dc = ((data[site]['SW_IN']<10) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['FC_USTAR_FLAG'])))).sum()
        # NEE_dc = ((data[site]['SW_IN']<10) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['NEE_USTAR_FLAG'])))).sum()

        omega_FC_filt = ((data[site]['SW_IN']<10) & (data[site]['OMEGA_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['FC_USTAR_FLAG']))) & (~pd.isnull((data[site]['FC'])))).sum()
        ustar_FC_filt = ((data[site]['SW_IN']<10) & (data[site]['FC_USTAR_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['FC_USTAR_FLAG']))) & (~pd.isnull((data[site]['FC'])))).sum()
        omega_NEE_filt = ((data[site]['SW_IN']<10) & (data[site]['OMEGA_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['NEE_USTAR_FLAG']))) & (~pd.isnull((data[site]['NEE'])))).sum()
        ustar_NEE_filt = ((data[site]['SW_IN']<10) & (data[site]['NEE_USTAR_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['NEE_USTAR_FLAG']))) & (~pd.isnull((data[site]['NEE'])))).sum()
        FC_dc = ((data[site]['SW_IN']<10) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['FC_USTAR_FLAG']))) & (~pd.isnull((data[site]['FC'])))).sum()
        NEE_dc = ((data[site]['SW_IN']<10) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['NEE_USTAR_FLAG']))) & (~pd.isnull((data[site]['NEE'])))).sum()
        # FC_dc = len(data[site])
        # NEE_dc = len(data[site])

        if FC_dc>200:
            data_cov['night'].loc[site,'OMEGA_filt_FC'] = omega_FC_filt/FC_dc
            data_cov['night'].loc[site,'USTAR_filt_FC'] = ustar_FC_filt/FC_dc
        if NEE_dc>200:
            data_cov['night'].loc[site,'OMEGA_filt_NEE'] = omega_NEE_filt/NEE_dc
            data_cov['night'].loc[site,'USTAR_filt_NEE'] = ustar_NEE_filt/NEE_dc
            
        data_cov['night'].loc[site,'LAI'] = data[site]['LAI'].mean()
        data_cov['night'].loc[site,'z'] = data[site]['z'].mean()
        data_cov['night'].loc[site,'canz'] = data[site]['canz'].mean()
        data_cov['night'].loc[site,'WEC'] = -data[site].loc[(data[site]['SW_IN']<10),'WEC'].mean()

        

        omega_FC_filt = ((data[site]['SW_IN']>10) & (data[site]['OMEGA_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['FC_USTAR_FLAG']))) & (~pd.isnull((data[site]['FC'])))).sum()
        ustar_FC_filt = ((data[site]['SW_IN']>10) & (data[site]['FC_USTAR_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['FC_USTAR_FLAG']))) & (~pd.isnull((data[site]['FC'])))).sum()
        omega_NEE_filt = ((data[site]['SW_IN']>10) & (data[site]['OMEGA_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['NEE_USTAR_FLAG']))) & (~pd.isnull((data[site]['NEE'])))).sum()
        ustar_NEE_filt = ((data[site]['SW_IN']>10) & (data[site]['NEE_USTAR_FLAG']==1) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['NEE_USTAR_FLAG']))) & (~pd.isnull((data[site]['NEE'])))).sum()
        FC_dc = ((data[site]['SW_IN']>10) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['FC_USTAR_FLAG']))) & (~pd.isnull((data[site]['FC'])))).sum()
        NEE_dc = ((data[site]['SW_IN']>10) & (~pd.isnull((data[site]['OMEGA_FLAG']))) & (~pd.isnull((data[site]['NEE_USTAR_FLAG']))) & (~pd.isnull((data[site]['NEE'])))).sum()
        # FC_dc = len(data[site])
        # NEE_dc = len(data[site])

        if FC_dc>200:
            data_cov['day'].loc[site,'OMEGA_filt_FC'] = omega_FC_filt/FC_dc
            data_cov['day'].loc[site,'USTAR_filt_FC'] = ustar_FC_filt/FC_dc
        if NEE_dc>200:
            data_cov['day'].loc[site,'OMEGA_filt_NEE'] = omega_NEE_filt/NEE_dc
            data_cov['day'].loc[site,'USTAR_filt_NEE'] = ustar_NEE_filt/NEE_dc
            
        data_cov['day'].loc[site,'LAI'] = data[site]['LAI'].mean()
        data_cov['day'].loc[site,'z'] = data[site]['z'].mean()
        data_cov['day'].loc[site,'canz'] = data[site]['canz'].mean()
        data_cov['day'].loc[site,'WEC'] = -data[site].loc[(data[site]['SW_IN']<10),'WEC'].mean()
    
    
    xvars = ['WEC','z']
    xvars = ['z']
    xvars = ['day','night']
    fig, axs = plt.subplot_mosaic([xvars],gridspec_kw={'wspace':0.05,'hspace':0.05,'right':0.95,'left':0.1,'top':0.95})
    if len(xvars)==1:
        fig.set_size_inches(figsize2[0]*cm, figsize2[1]*cm,forward=True)
    else:
        fig.set_size_inches(figsize3[0]*cm, figsize3[1]*cm,forward=True)
    for xvar in xvars:
        filtered = lowess(data_cov[xvar]['OMEGA_filt_NEE'], data_cov[xvar]['z'],frac=0.7)
        filtered2 = lowess(data_cov[xvar]['USTAR_filt_NEE'], data_cov[xvar]['z'],frac=0.7)
        axs[xvar].plot(data_cov[xvar]['z'],data_cov[xvar]['USTAR_filt_NEE']*1e2,'ro',markerfacecolor='w',label='u$_{*}$ filtering')
        axs[xvar].plot(filtered2[:,0],filtered2[:,1]*1e2,'r-')
        
        axs[xvar].plot(data_cov[xvar]['z'],data_cov[xvar]['OMEGA_filt_NEE']*1e2,'k^',markerfacecolor='w',label='$\Omega$ filtering')
        axs[xvar].plot(filtered[:,0],filtered[:,1]*1e2,'k-')
        axs[xvar].set_xlabel('EC height (m)')
        # ax['z'].set_title('Night time')
            
        if xvar==xvars[0]:
            axs[xvar].legend(loc="upper left")
            axs[xvar].set_ylabel('Data filtered out (%)')
        else:
            axs[xvar].set_yticklabels([])
        axs[xvar].set_ylim([0,100])

    axs['day'].set_title('Daytime')
    axs['night'].set_title('Night time')
    
    fig.text(.95, .05, al[0]+')', ha='right', va='bottom', transform=axs['day'].transAxes,backgroundcolor='w')    
    fig.text(.95, .05, al[1]+')', ha='right', va='bottom', transform=axs['night'].transAxes,backgroundcolor='w')    

    # plt.show()
    fig.savefig(figfold+'/fig12.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig12.pdf', dpi=200,bbox_inches='tight')



def figC15(data,clrs,forest_site,grass_site):


    def sample_heights(zvec,n,ztop):
        zvecout = pd.Series()
        n = int(n)
        if n<2:
            print('Give higher n')
        else:
            z_sampled = np.linspace(0,ztop,n)
            # z_sampled = np.logspace(np.log10(5e-1),np.log10(ztop),n)
            for z1 in z_sampled:

                zdiff_sorted = (zvec-z1).abs().sort_values()
                zvalout = None
                indx = 0
                while (zvalout is None) & (indx<len(zdiff_sorted)-1):
                    if len(zvecout)>0:
                        if zvecout.isin([zvec[zdiff_sorted.index[indx]]]).any():
                        # if zvec.iloc[indx] not in zvecout:
                            indx = indx + 1
                        else:
                            zvalindex = zdiff_sorted.index[zdiff_sorted.iloc[indx]==zdiff_sorted]
                            zvalout = zvec[zvalindex]
                    else:
                        zvalindex = zdiff_sorted.index[zdiff_sorted.iloc[indx]==zdiff_sorted]
                        zvalout = zvec[zvalindex]
                # print(str(z1) + ' ' + str(zvalout.iloc[0]) + ' ' + str(len(zvalout)))
                if len(zvalout)>1:
                    zvalout = pd.Series(index=[zvalout.index[0]],data=zvalout.iloc[0])
                if z1==z_sampled[0]:
                    zvecout = zvalout
                else:
                    zvecout = pd.concat([zvecout,zvalout])
        zvecout = zvecout.sort_values()
        return zvecout

    clrs = dict()
    clrs['TP'] = cmp.vik(0.2)
    clrs['TN'] = cmp.vik(0.4)
    clrs['FP'] = cmp.vik(0.6)
    clrs['FN'] = cmp.vik(0.8)

    labels = dict()
    labels['TP'] = 'Coupled-Coupled'
    labels['TN'] = 'Decoupled-Decoupled'
    labels['FP'] = 'Decoupled-Coupled'
    labels['FN'] = 'Coupled-Decoupled'

    ticklabels = dict()
    ticklabels['constant_LAI'] = 'constant PAI'
    ticklabels['05_LAI'] = '0.5*PAI'
    ticklabels['02_LAI'] = '0.8*PAI'
    ticklabels['12_LAI'] = '1.2*PAI'
    ticklabels['15_LAI'] = '1.5*PAI'
    for i in np.linspace(2,20,19):
        ticklabels['T_'+str(int(i))] = str(int(i)) + ' T heights'


    # forest site
    TPOTcols = list(data[forest_site].columns[data[forest_site].columns.str.startswith('TPOT_')])
    TPOTzcols = ['z_'+x for x in TPOTcols]
    zTPOT = data[forest_site][TPOTzcols].iloc[0]
    zTPOT.index = zTPOT.index.str.replace('z_','')
    # zTPOT = dict(zTPOT)
    
    omega_values = dict()
    omega_values[forest_site] = dict()
    omega_values[forest_site]['ref'] = data[forest_site]['OMEGA']
    # constant LAI
    LAI = data[forest_site]['LAI'].mean()
    # omega2,wec2,reftheta2 = utils.calculate_omega_ver2(data[forest_site]['W_SIGMA'],data[forest_site].loc[:,TPOTcols]+c.NT,data[forest_site]['Uh'],LAI,data[forest_site]['canz'].iloc[0],data[forest_site]['z'].iloc[0],dict(zTPOT))
    omega,wec,reftheta = utils.calculate_omega(data[forest_site]['W_SIGMA'],data[forest_site].loc[:,TPOTcols]+c.NT,data[forest_site]['Uh'],LAI,data[forest_site]['canz'].iloc[0],data[forest_site]['z'].iloc[0],dict(zTPOT))
    omega_values[forest_site]['constant_LAI'] = omega
    # 50 % underestimated LAI
    LAI = data[forest_site]['LAI']*0.5
    omega,wec,reftheta = utils.calculate_omega(data[forest_site]['W_SIGMA'],data[forest_site].loc[:,TPOTcols]+c.NT,data[forest_site]['Uh'],LAI,data[forest_site]['canz'].iloc[0],data[forest_site]['z'].iloc[0],dict(zTPOT))
    omega_values[forest_site]['05_LAI'] = omega
    # 20 % underestimated LAI
    LAI = data[forest_site]['LAI']*0.8
    omega,wec,reftheta = utils.calculate_omega(data[forest_site]['W_SIGMA'],data[forest_site].loc[:,TPOTcols]+c.NT,data[forest_site]['Uh'],LAI,data[forest_site]['canz'].iloc[0],data[forest_site]['z'].iloc[0],dict(zTPOT))
    omega_values[forest_site]['02_LAI'] = omega
    # 20 % overestimated LAI
    LAI = data[forest_site]['LAI']*1.2
    omega,wec,reftheta = utils.calculate_omega(data[forest_site]['W_SIGMA'],data[forest_site].loc[:,TPOTcols]+c.NT,data[forest_site]['Uh'],LAI,data[forest_site]['canz'].iloc[0],data[forest_site]['z'].iloc[0],dict(zTPOT))
    omega_values[forest_site]['12_LAI'] = omega
    # 50 % overestimated LAI
    LAI = data[forest_site]['LAI']*1.5
    omega,wec,reftheta = utils.calculate_omega(data[forest_site]['W_SIGMA'],data[forest_site].loc[:,TPOTcols]+c.NT,data[forest_site]['Uh'],LAI,data[forest_site]['canz'].iloc[0],data[forest_site]['z'].iloc[0],dict(zTPOT))
    omega_values[forest_site]['15_LAI'] = omega
    # T measurement heights
    # for n in np.linspace(2,len(TPOTzcols),5):
    for n in np.linspace(2,len(TPOTzcols)-3,len(TPOTzcols)-5):
        if int(n)<len(TPOTzcols):
            if ('T_'+str(int(n))) not in omega_values[forest_site].keys():
                zTPOT2 = sample_heights(zTPOT,int(n),data[forest_site]['z'].iloc[0])
                TPOTcols2 = list(zTPOT2.index)
                TPOTcols2 = [sub.replace('z_', '') for sub in TPOTcols2]
                omega,wec,reftheta = utils.calculate_omega(data[forest_site]['W_SIGMA'],data[forest_site].loc[:,TPOTcols2]+c.NT,data[forest_site]['Uh'],data[forest_site]['LAI'],data[forest_site]['canz'].iloc[0],data[forest_site]['z'].iloc[0],dict(zTPOT2))
                omega_values[forest_site]['T_'+str(int(n))] = omega
    
    # grass site
    TPOTcols = list(data[grass_site].columns[data[grass_site].columns.str.startswith('TPOT_')])
    TPOTzcols = ['z_'+x for x in TPOTcols]
    zTPOT = data[grass_site][TPOTzcols].iloc[0]
    zTPOT.index = zTPOT.index.str.replace('z_','')

    omega_values[grass_site] = dict()
    omega_values[grass_site]['ref'] = data[grass_site]['OMEGA']
    # T measurement heights
    for n in np.linspace(2,len(TPOTzcols),len(TPOTzcols)-1):
        if int(n)<len(TPOTzcols):
            if ('T_'+str(int(n))) not in omega_values[grass_site].keys():
                zTPOT2 = sample_heights(zTPOT,int(n),data[grass_site]['z'].iloc[0])
                TPOTcols2 = list(zTPOT2.index)
                TPOTcols2 = [sub.replace('z_', '') for sub in TPOTcols2]
                omega,wec,reftheta = utils.calculate_omega(data[grass_site]['W_SIGMA'],data[grass_site].loc[:,TPOTcols2]+c.NT,data[grass_site]['Uh'],data[grass_site]['LAI'],data[grass_site]['canz'].iloc[0],data[grass_site]['z'].iloc[0],dict(zTPOT2))
                omega_values[grass_site]['T_'+str(int(n))] = omega
            


    MCC_values = dict()
    accuracy = dict()
    for site in omega_values.keys():
        MCC_values[site] = dict()
        accuracy[site] = dict()
        for key in omega_values[site].keys():
            ref_flag = pd.Series(index=data[site].index,data=np.nan)
            ref_flag[omega_values[site]['ref']<OMEGAth] = 1
            ref_flag[omega_values[site]['ref']>=OMEGAth] = 0
            # ref_flag[data[site]['SW_IN']>=10] = np.nan
            if key!='ref':
                MCC_values[site][key] = pd.Series()
                flag = pd.Series(index=data[site].index,data=np.nan)
                flag[omega_values[site][key]<OMEGAth] = 1
                flag[omega_values[site][key]>=OMEGAth] = 0
                # flag[data[site]['SW_IN']>=10] = np.nan

                indx = (~ref_flag.isna()) & (~flag.isna())
                ref_flag = ref_flag[indx]
                flag = flag[indx]

                MCC_values[site][key]['MCC'] = utils.calculate_MCC(ref_flag,flag)
                MCC_values[site][key]['TP'] = ((ref_flag==1) & (flag==1)).sum()/len(ref_flag)
                MCC_values[site][key]['TN'] = ((ref_flag==0) & (flag==0)).sum()/len(ref_flag)
                MCC_values[site][key]['FP'] = ((ref_flag==0) & (flag==1)).sum()/len(ref_flag)
                MCC_values[site][key]['FN'] = ((ref_flag==1) & (flag==0)).sum()/len(ref_flag)
                accuracy[site][key] = MCC_values[site][key]['TP']+MCC_values[site][key]['TN']


    fig, ax = plt.subplot_mosaic([[grass_site],[forest_site],[forest_site],[forest_site]],gridspec_kw={'wspace':0.05,'hspace':0.5,'right':0.95,'left':0.2,'top':0.9,'bottom':0.15})
    fig.set_size_inches(figsize2[0]*cm, figsize2[1]*cm,forward=True)
    
    for site in omega_values.keys():
        ticklbl = []
        for indx in range(len(MCC_values[site])):
            key = list(MCC_values[site].keys())[indx]
            ticklbl.append(ticklabels[key])
            bottom = 0
            for var in MCC_values[site][key].index:
                if var!='MCC':
                    if indx==0:
                        ax[site].barh(indx,100*MCC_values[site][key][var],left=bottom,color=clrs[var],edgecolor='k',label=labels[var])
                    else:
                        ax[site].barh(indx,100*MCC_values[site][key][var],left=bottom,color=clrs[var],edgecolor='k')
                    bottom = bottom+100*MCC_values[site][key][var]
        ax[site].set_xlim([0,100])
        
        ax[site].set_yticks(np.linspace(0,len(MCC_values[site])-1,len(MCC_values[site])))
        # ax[site].set_xticklabels(ticklbl,rotation=45)
        ax[site].set_yticklabels(ticklbl)

        if site!=grass_site:
            ax[site].set_xlabel('Fraction of data (%)')
            ax[site].legend(loc='upper left')
        else:
            ax[site].set_xticklabels([])


    if config['source'][forest_site]=='NEON':
        ax[forest_site].set_title('Forest site (z=' + f"{data[forest_site]['z'].iloc[0]: .1f}" + ' m), ' + utils.NEON_to_ameriflux(forest_site))
    else:
        ax[forest_site].set_title('Forest site (z=' + f"{data[forest_site]['z'].iloc[0]: .1f}" + ' m), ' + forest_site)
    # ax[forest_site].set_title('Forest site (z=' + f"{data[forest_site]['z'].iloc[0]: .1f}" + ' m), ' + forest_site)
    # ax[grass_site].set_title('Grassland site (z=' + f"{data[grass_site]['z'].iloc[0]: .1f}" + ' m), ' + grass_site)
    if config['source'][grass_site]=='NEON':
        ax[grass_site].set_title('Grassland site (z=' + f"{data[grass_site]['z'].iloc[0]: .1f}" + ' m), ' + utils.NEON_to_ameriflux(grass_site))
    else:
        ax[grass_site].set_title('Grassland site (z=' + f"{data[grass_site]['z'].iloc[0]: .1f}" + ' m), ' + grass_site)


    # plt.show()


    fig.savefig(figfold+'/figC15.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/figC15.pdf', dpi=200,bbox_inches='tight')




    fig, ax = plt.subplot_mosaic([[grass_site],[forest_site],[forest_site],[forest_site]],gridspec_kw={'wspace':0.05,'hspace':0.7,'right':0.95,'left':0.2,'top':0.9,'bottom':0.15})
    fig.set_size_inches(figsize2[0]*cm, figsize2[1]*cm,forward=True)
    
    for site in omega_values.keys():
        ticklbl = []
        for indx in range(len(MCC_values[site])):
            key = list(MCC_values[site].keys())[indx]
            ticklbl.append(ticklabels[key])
            bottom = 0
            for var in MCC_values[site][key].index:
                if var=='MCC':
                    if indx==0:
                        ax[site].plot(MCC_values[site][key][var],indx,'o',markeredgecolor='k',markerfacecolor='w')
                    else:
                        ax[site].plot(MCC_values[site][key][var],indx,'o',markeredgecolor='k',markerfacecolor='w')
        
        ax[site].set_xlim([0.6,1])
        ax[site].set_yticks(np.linspace(0,len(MCC_values[site])-1,len(MCC_values[site])))
        # ax[site].set_xticklabels(ticklbl,rotation=45)
        ax[site].set_yticklabels(ticklbl)

        if site!=grass_site:
            ax[site].set_xlabel('MCC (-)')
        else:
            ax[site].set_xticklabels([])


    ax[forest_site].set_title('Forest site (z=' + f"{data[forest_site]['z'].iloc[0]: .1f}" + ' m), ' + forest_site)
    ax[grass_site].set_title('Grassland site (z=' + f"{data[grass_site]['z'].iloc[0]: .1f}" + ' m), ' + grass_site)

    fig.savefig(figfold+'/fig9_ver2.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig9_ver2.pdf', dpi=200,bbox_inches='tight')

    # plt.show()

def fig11(data_all,clrs,sitesplot=None):

    if sitesplot is None:
        sitesplot = list(data_all.keys())

    wec = dict()
    flag = dict()
    flag2 = dict()
    for site in sitesplot:
        wec[site] = pd.DataFrame(index=data_all[site].index,columns=['buoyancy','drag'])
        flag[site] = pd.Series(index=data_all[site].index)
        flag2[site] = pd.Series(index=data_all[site].index)

        TPOTcols = list(data_all[site].columns[data_all[site].columns.str.startswith('TPOT_')])
        TPOTzcols = ['z_'+x for x in TPOTcols]
        zTPOT = data_all[site][TPOTzcols].iloc[0]
        zTPOT.index = zTPOT.index.str.replace('z_','')

        if data_all[site]['canz'].iloc[0]>1:
                #removing canopy drag
                omega,wec1,reftheta = utils.calculate_omega(data_all[site]['W_SIGMA'],data_all[site].loc[:,TPOTcols]+c.NT,data_all[site]['Uh'],0,data_all[site]['canz'].iloc[0],data_all[site]['z'].iloc[0],dict(zTPOT))
                wec[site]['drag'] = wec1

                #removing buoyancy
                TPOTmat = pd.DataFrame(index=data_all[site].index)
                for col in TPOTcols:
                    TPOTmat[col] = data_all[site].loc[:,TPOTcols].mean(axis=1)+c.NT
                    
                omega,wec1,reftheta = utils.calculate_omega(data_all[site]['W_SIGMA'],TPOTmat,data_all[site]['Uh'],data_all[site]['LAI'],data_all[site]['canz'].iloc[0],data_all[site]['z'].iloc[0],dict(zTPOT))
                wec[site]['buoyancy'] = wec1

                flag[site][wec[site]['drag']<wec[site]['buoyancy']] = 0
                flag[site][wec[site]['drag']>wec[site]['buoyancy']] = 1
                flag[site][data_all[site]['OMEGA']>OMEGAth] = 2

                
                flag2[site][data_all[site]['OMEGA']<OMEGAth] = 3
                flag2[site][wec[site]['buoyancy']/data_all[site]['WEC']>0.8] = 1
                flag2[site][wec[site]['drag']/data_all[site]['WEC']>0.8] = 2
                flag2[site][data_all[site]['OMEGA']>OMEGAth] = 4

    xlm = [dt.datetime(2019,1,1),dt.datetime(2020,12,31,23,30,0)]
    fig, ax = plt.subplot_mosaic([sitesplot,['cb' for i in range(len(sitesplot))]],gridspec_kw={'bottom':0.17,'top':0.97,'right':0.97,'wspace':0.04,'hspace':0.24,"height_ratios":[1, 0.05]})
    fig.set_size_inches(figsize3[0]*cm, figsize3[1]*cm,forward=True)
    locator = mpl.dates.MonthLocator((1,5,9))
    fmt = mpl.dates.DateFormatter('%m/%Y')
    for i in range(len(sitesplot)):
        site = sitesplot[i]
        if data_all[site]['canz'].iloc[0]>1:
            dat = data_all[site].copy()
            wec_tmp = wec[site].copy()
            flag2_tmp = flag2[site].copy()
            if site=='SERC':
                dat.index = dat.index-dt.timedelta(hours=5)
                wec_tmp.index = wec_tmp.index-dt.timedelta(hours=5)
                flag2_tmp.index = flag2_tmp.index-dt.timedelta(hours=5)



            tstart = pd.to_datetime(flag2_tmp.index[0]).floor('1D')
            SWINmat,datevec,time2,xlms = utils.reorganise_for_fingerprint(dat['SW_IN'],tstart)


            ratio = wec_tmp['drag']/dat['WEC']
            ratio[dat['WEC']==0] = 0
            zmat,datevec,time2,xlms = utils.reorganise_for_fingerprint(flag2_tmp,tstart)
            levels = mpl.ticker.MaxNLocator(nbins=4).tick_values(1, 4)
            cmap = cmp.batlowW
            norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            datevec = flag2_tmp.index[0]+pd.to_timedelta(datevec,unit='day')
            
            pcm = ax[site].pcolormesh(datevec,time2,zmat,cmap=cmap, norm=norm,rasterized=True)
            
            ax[site].contour(datevec,time2,SWINmat, levels=[10],linestyles='dashed',colors='k',linewidths=0.7)
            
            ax[site].patch.set_facecolor('grey')
            # ax[site].set_xlim([np.min(datevec[~np.all(np.isnan(zmat),0)]),np.max(datevec[~np.all(np.isnan(zmat),0)])])
            ax[site].set_xlim(xlm)
            if site==sitesplot[0]:
                ax[site].set_ylabel('Local time')
            else:
                ax[site].set_yticklabels([])
            # ax[site].tick_params(axis='both', labelsize=6)

            if config['source'][site]=='NEON':
                ax[site].set_title(utils.NEON_to_ameriflux(site))
            else:
                ax[site].set_title(site)
            
            ax[site].xaxis.set_major_locator(locator)
            ax[site].xaxis.set_major_formatter(fmt)
            if site==sitesplot[0]:
                cb = fig.colorbar(pcm,cax=ax['cb'], orientation="horizontal")
                # cb.set_ticks(ticks=[1.2,2,2.8,3.6], labels=['Decoupled,\ndrag','Decoupled,\nbuoyancy','Decoupled,\ndrag+buoyancy','Coupled'],fontsize=6)
                cb.set_ticks(ticks=[1.2,2,2.8,3.6], labels=['Decoupled,\ndrag','Decoupled,\nbuoyancy','Decoupled,\ndrag+buoyancy','Coupled'])
                # labels = cb.ax.get_xticklabels()
                # cb.set_ticks(ticks=[0,0.2,0.4,0.6,0.8,1],labels=labels,fontsize=6)
                # cb.set_label('w$_{e,crit,buoyancy}/$w$_{e,crit}$ (-)')


            fig.text(.95, .05, al[i]+')', ha='right', va='bottom', transform=ax[site].transAxes,backgroundcolor='w')    
            # plt.show()

            del flag2_tmp,wec_tmp,dat

        fig.savefig(figfold+'/fig11.png', dpi=200,bbox_inches='tight')
        fig.savefig(figfold+'/fig11.pdf', dpi=200,bbox_inches='tight')


def fig13_14(data_all,clrs,sitesplot=None):

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    if sitesplot is None:
        # sitesplot = list(data_all.keys())

        z = pd.Series(index=data_all.keys())
        LAIpeak = pd.Series(index=data_all.keys())
        LAIstd = pd.Series(index=data_all.keys())
        for site in data_all.keys():
            z[site] = data_all[site]["z"].mean()
            LAIpeak[site] = data_all[site]["LAI"].max()
            LAIstd[site] = data_all[site]["LAI"].std()
        sitesplot = list(z.sort_values().index)
        # sitesplot = list(LAIpeak.sort_values().index)
        # sitesplot = list(LAIstd.sort_values().index)
        del z, LAIpeak, LAIstd

    sitesplot1 = []
    for site in sitesplot:
        if ((~data_all[site]['NEE'].isnull()) & (~data_all[site]['NEE_USTAR_FLAG'].isnull()) & (~data_all[site]['OMEGA_FLAG'].isnull())).sum()>100:
            sitesplot1.append(site)

    night_df = pd.DataFrame(index = sitesplot1,columns=['pvalue','mean_diff','median_diff','std_diff','iqr_diff'])
    day_df = pd.DataFrame(index = sitesplot1,columns=['pvalue','mean_diff','median_diff','std_diff','iqr_diff'])
    for i in range(len(sitesplot1)):
        site = sitesplot1[i]
        indx = (data_all[site]['SW_IN']<10) & (data_all[site]['NEE_USTAR_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        indx_omega = (data_all[site]['SW_IN']<10) & (data_all[site]['OMEGA_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        # indx = ((data_all[site].index.month<5) | (data_all[site].index.month>9)) & (data_all[site]['NEE_USTAR_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        # indx_omega = ((data_all[site].index.month<5) | (data_all[site].index.month>9)) & (data_all[site]['OMEGA_FLAG']==0) & (~data_all[site]['NEE'].isnull())

        x_ust = data_all[site].loc[indx,'NEE']
        x_omega = data_all[site].loc[indx_omega,'NEE']

        kstest_results = scipy.stats.kstest(x_ust,x_omega)
        night_df.loc[site,'pvalue'] = kstest_results.pvalue
        night_df.loc[site,'mean_diff'] = (x_ust.mean()-x_omega.mean())/x_ust.mean()
        night_df.loc[site,'median_diff'] = x_ust.median()-x_omega.median()
        night_df.loc[site,'std_diff'] = (x_ust.std()-x_omega.std())/x_ust.std()
        night_df.loc[site,'iqr_diff'] = (x_ust.quantile(q=0.75)-x_ust.quantile(q=0.25))-(x_omega.quantile(q=0.75)-x_omega.quantile(q=0.25))
        
        indx = (data_all[site]['SW_IN']>10) & (data_all[site]['NEE_USTAR_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        indx_omega = (data_all[site]['SW_IN']>10) & (data_all[site]['OMEGA_FLAG']==0) & (~data_all[site]['NEE'].isnull())

        x_ust = data_all[site].loc[indx,'NEE']
        x_omega = data_all[site].loc[indx_omega,'NEE']

        kstest_results = scipy.stats.kstest(x_ust,x_omega)
        day_df.loc[site,'pvalue'] = kstest_results.pvalue
        day_df.loc[site,'mean_diff'] = (x_ust.mean()-x_omega.mean())/x_ust.mean()
        day_df.loc[site,'median_diff'] = x_ust.median()-x_omega.median()
        day_df.loc[site,'std_diff'] = (x_ust.std()-x_omega.std())/x_ust.std()
        day_df.loc[site,'iqr_diff'] = (x_ust.quantile(q=0.75)-x_ust.quantile(q=0.25))-(x_omega.quantile(q=0.75)-x_omega.quantile(q=0.25))

    asd = pd.Series(index=data_all.keys())
    for site in sitesplot:
        asd[site] = data_all[site]['LAI'].max()-data_all[site]['LAI'].min()

    mosaic = []
    indx = 0
    for i in range(7):
        row = []
        for j in range(6):
            if indx<=len(sitesplot1):
                row.append(sitesplot1[indx])
            else:
                row.append('X')
            indx = indx + 1
        mosaic.append(row)


    # fig, ax = plt.subplot_mosaic(mosaic,empty_sentinel='X',gridspec_kw={'bottom':0.17,'top':0.97,'right':0.97,'wspace':0,'hspace':0})
    fig, ax = plt.subplot_mosaic(mosaic,empty_sentinel='X',gridspec_kw={'wspace':0,'hspace':0})
    fig.set_size_inches(21*cm, 29.7*0.7*cm,forward=True)
    for i in range(len(sitesplot1)):
        site = sitesplot1[i]
        # indx = (data_all[site]['SW_IN']<10) & (data_all[site]['NEE_USTAR_FLAG']==0) & (~data_all[site]['NEE'].isnull()) & (data_all[site].index.month>4) & (data_all[site].index.month<10)
        # indx_omega = (data_all[site]['SW_IN']<10) & (data_all[site]['OMEGA_FLAG']==0) & (~data_all[site]['NEE'].isnull()) & (data_all[site].index.month>4) & (data_all[site].index.month<10)
        indx = (data_all[site]['SW_IN']<10) & (data_all[site]['NEE_USTAR_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        indx_omega = (data_all[site]['SW_IN']<10) & (data_all[site]['OMEGA_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        # indx = (data_all[site].index.month>4) & (data_all[site].index.month<10) & (data_all[site]['NEE_USTAR_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        # indx_omega = (data_all[site].index.month>4) & (data_all[site].index.month<10) & (data_all[site]['OMEGA_FLAG']==0) & (~data_all[site]['NEE'].isnull())


        x_ust = data_all[site].loc[indx,'NEE']
        x_omega = data_all[site].loc[indx_omega,'NEE']

        x_mean = x_ust.mean()
        x_std = x_ust.std()
        x_ust = (x_ust-x_mean)/x_std
        x_omega = (x_omega-x_mean)/x_std
        used_quantiles = np.linspace(.01,.99,99)
        q_ust = x_ust.quantile(q=used_quantiles)
        q_omega = x_omega.quantile(q=used_quantiles)

        density_omega = gaussian_kde(x_omega)
        density_ust = gaussian_kde(x_ust)
        xvals = np.linspace(-4,4,100)
        ylms = [0,np.max([np.max(density_omega(xvals)),np.max(density_ust(xvals))])]
        ax[site].plot([0,0],[-4,4],'--',color=(0.5,0.5,0.5))
        # ax[site].boxplot(x_omega,showfliers=False,vert=False,positions=[ylms[1]*0.1],widths=[ylms[1]*0.05],boxprops={'color':colors[0]},medianprops={'color':colors[0]})
        # ax[site].boxplot(x_ust,showfliers=False,vert=False,positions=[ylms[1]*0.2],widths=[ylms[1]*0.05],boxprops={'color':colors[1]},medianprops={'color':colors[1]})
        bp_omega = ax[site].boxplot(x_omega, patch_artist=True,showfliers=False,vert=False,positions=[ylms[1]*0.1],widths=[ylms[1]*0.05],capprops={'linewidth':0.5},whiskerprops={'linewidth':0.5},boxprops={'color':'k','linewidth':0.5},medianprops={'color':'k','linewidth':0.5})
        bp_ust = ax[site].boxplot(x_ust, patch_artist=True,showfliers=False,vert=False,positions=[ylms[1]*0.2],widths=[ylms[1]*0.05],capprops={'linewidth':0.5},whiskerprops={'linewidth':0.5},boxprops={'color':'k','linewidth':0.5},medianprops={'color':'k','linewidth':0.5})
        for patch in bp_omega['boxes']:
            patch.set(facecolor=colors[0])
        for patch in bp_ust['boxes']:
            patch.set(facecolor=colors[1])
        ax[site].plot(xvals,density_omega(xvals),label="$\Omega$ filtering",color=colors[0])
        ax[site].plot(xvals,density_ust(xvals),label="$u_*$ filtering",color=colors[1])
        # ax[site].hist(x_ust,bins=np.linspace(-4,4,10),density=True,alpha)

        # ax[site].plot([-4,4],[-4,4],'k--')
        # ax[site].plot([-4,4],[0,0],'--',color=(0.5,0.5,0.5))
        # ax[site].plot(q_ust,q_omega,'k.')
        # ax[site].plot(x_ust.median(),x_omega.median(),'ko',markerfacecolor='r')
        # ax[site].plot(x_ust.mean(),x_omega.mean(),'k^',markerfacecolor='b')
        ax[site].set_xlim([-4,4])
        ax[site].set_ylim(ylms)
        ax[site].set_yticks([])
        # ax[site].set_ylim([-4,4])
        ax[site].set_xticks([-3,-2,-1,0,1,2,3])
        ax[site].tick_params(axis="x", labelsize=8)
        # ax[site].set_yticks([-4,-2,0,2,4])
        if i==0:
            ax[site].legend(loc='upper left',ncols=2,bbox_to_anchor=(0.0, 1.35))
        txt = ""
        if config['source'][site]=='NEON':
            txt = txt + utils.NEON_to_ameriflux(site) + "\n"
        else:
            txt = txt + site + "\n"
        txt = txt + "z=" + f"{data_all[site]['z'].mean():.1f}" + " m\n"
        if data_all[site]['LAI'].min()==data_all[site]['LAI'].max():
            PAItxt = "PAI=" + f"{data_all[site]['LAI'].mean():.1f}" + "\n"
        else:
            PAItxt = "PAI=" + f"{data_all[site]['LAI'].min():.1f}" + "..." + f"{data_all[site]['LAI'].max():.1f}" + "\n"
        txt = txt + PAItxt
        if night_df.loc[site,'pvalue']<0.05:
            fig.text(.95, .95, txt, ha='right', va='top', transform=ax[site].transAxes,fontsize=6,color='r')    
        else:            
            fig.text(.95, .95, txt, ha='right', va='top', transform=ax[site].transAxes,fontsize=6) 

        left_bool = False
        bottom_bool = False
        for i in range(len(mosaic)):
            
            row = mosaic[i]

            if (any(site in x for x in row)) & (i==len(mosaic)-1):
                bottom_bool = True 
            
            if (row[0]==site):
                left_bool = True

        if not bottom_bool:
            ax[site].set_xticklabels([])
        if not left_bool:
            ax[site].set_yticklabels([])

        # fig.suptitle('Night time')
    
    # fig.supylabel('Standardised night time F$_c$+S$_c$, $\Omega$ filtering (-)',x=0.06)
    fig.supylabel('Probability density function (-)',x=0.09)
    fig.supxlabel('Standardised night time F$_c$+S$_c$ (-)',y=0.06)
    fig.savefig(figfold+'/fig13.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig13.pdf', dpi=200,bbox_inches='tight')
    

    # fig, ax = plt.subplot_mosaic(mosaic,empty_sentinel='X',gridspec_kw={'bottom':0.17,'top':0.97,'right':0.97,'wspace':0,'hspace':0})
    fig, ax = plt.subplot_mosaic(mosaic,empty_sentinel='X',gridspec_kw={'wspace':0,'hspace':0})
    fig.set_size_inches(21*cm, 29.7*0.7*cm,forward=True)
    for i in range(len(sitesplot1)):
        site = sitesplot1[i]
        # indx = (data_all[site]['SW_IN']>10) & (data_all[site]['NEE_USTAR_FLAG']==0) & (~data_all[site]['NEE'].isnull()) & (data_all[site].index.month>4) & (data_all[site].index.month<10)
        # indx_omega = (data_all[site]['SW_IN']>10) & (data_all[site]['OMEGA_FLAG']==0) & (~data_all[site]['NEE'].isnull()) & (data_all[site].index.month>4) & (data_all[site].index.month<10)
        # indx = ((data_all[site].index.month<5) | (data_all[site].index.month>9)) & (data_all[site]['NEE_USTAR_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        # indx_omega = ((data_all[site].index.month<5) | (data_all[site].index.month>9)) & (data_all[site]['OMEGA_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        indx = (data_all[site]['SW_IN']>10) & (data_all[site]['NEE_USTAR_FLAG']==0) & (~data_all[site]['NEE'].isnull())
        indx_omega = (data_all[site]['SW_IN']>10) & (data_all[site]['OMEGA_FLAG']==0) & (~data_all[site]['NEE'].isnull())

        x_ust = data_all[site].loc[indx,'NEE']
        x_omega = data_all[site].loc[indx_omega,'NEE']

        x_mean = x_ust.mean()
        x_std = x_ust.std()
        x_ust = (x_ust-x_mean)/x_std
        x_omega = (x_omega-x_mean)/x_std
        used_quantiles = np.linspace(.01,.99,99)
        q_ust = x_ust.quantile(q=used_quantiles)
        q_omega = x_omega.quantile(q=used_quantiles)
        density_omega = gaussian_kde(x_omega)
        density_ust = gaussian_kde(x_ust)
        xvals = np.linspace(-4,4,100)
        ylms = [0,np.max([np.max(density_omega(xvals)),np.max(density_ust(xvals))])]
        ax[site].plot([0,0],[-4,4],'--',color=(0.5,0.5,0.5))
        # ax[site].boxplot(x_omega,showfliers=False,vert=False,positions=[ylms[1]*0.1],widths=[ylms[1]*0.05],boxprops={'color':colors[0]},medianprops={'color':colors[0]})
        # ax[site].boxplot(x_ust,showfliers=False,vert=False,positions=[ylms[1]*0.2],widths=[ylms[1]*0.05],boxprops={'color':colors[1]},medianprops={'color':colors[1]})
        bp_omega = ax[site].boxplot(x_omega, patch_artist=True,showfliers=False,vert=False,positions=[ylms[1]*0.1],widths=[ylms[1]*0.05],capprops={'linewidth':0.5},whiskerprops={'linewidth':0.5},boxprops={'color':'k','linewidth':0.5},medianprops={'color':'k','linewidth':0.5})
        bp_ust = ax[site].boxplot(x_ust, patch_artist=True,showfliers=False,vert=False,positions=[ylms[1]*0.2],widths=[ylms[1]*0.05],capprops={'linewidth':0.5},whiskerprops={'linewidth':0.5},boxprops={'color':'k','linewidth':0.5},medianprops={'color':'k','linewidth':0.5})
        for patch in bp_omega['boxes']:
            patch.set(facecolor=colors[0])
        for patch in bp_ust['boxes']:
            patch.set(facecolor=colors[1])
        ax[site].plot(xvals,density_omega(xvals),label="$\Omega$ filtering",color=colors[0])
        ax[site].plot(xvals,density_ust(xvals),label="$u_*$ filtering",color=colors[1])
        # ax[site].hist(x_ust,bins=np.linspace(-4,4,10),density=True,alpha)

        # ax[site].plot([-4,4],[-4,4],'k--')
        # ax[site].plot([-4,4],[0,0],'--',color=(0.5,0.5,0.5))
        # ax[site].plot(q_ust,q_omega,'k.')
        # ax[site].plot(x_ust.median(),x_omega.median(),'ko',markerfacecolor='r')
        # ax[site].plot(x_ust.mean(),x_omega.mean(),'k^',markerfacecolor='b')
        ax[site].set_xlim([-4,4])
        ax[site].set_ylim([0,np.max([np.max(density_omega(xvals)),np.max(density_ust(xvals))])])
        # ax[site].set_ylim([-4,4])
        ax[site].set_xticks([-3,-2,-1,0,1,2,3])
        ax[site].tick_params(axis="x", labelsize=8)
        ax[site].set_yticks([])
        # ax[site].set_yticks([-4,-2,0,2,4])
        if i==0:
            ax[site].legend(loc='upper left',ncols=2,bbox_to_anchor=(0.0, 1.35))
        txt = ""
        if config['source'][site]=='NEON':
            txt = txt + utils.NEON_to_ameriflux(site) + "\n"
        else:
            txt = txt + site + "\n"
        txt = txt + "z=" + f"{data_all[site]['z'].mean():.1f}" + " m\n"
        if data_all[site]['LAI'].min()==data_all[site]['LAI'].max():
            PAItxt = "PAI=" + f"{data_all[site]['LAI'].mean():.1f}" + "\n"
        else:
            PAItxt = "PAI=" + f"{data_all[site]['LAI'].min():.1f}" + "..." + f"{data_all[site]['LAI'].max():.1f}" + "\n"
        txt = txt + PAItxt
        if day_df.loc[site,'pvalue']<0.05:
            fig.text(.05, .95, txt, ha='left', va='top', transform=ax[site].transAxes,fontsize=6,color='r')    
        else:            
            fig.text(.05, .95, txt, ha='left', va='top', transform=ax[site].transAxes,fontsize=6) 
        left_bool = False
        bottom_bool = False
        for i in range(len(mosaic)):
            
            row = mosaic[i]

            if (any(site in x for x in row)) & (i==len(mosaic)-1):
                bottom_bool = True 
            
            if (row[0]==site):
                left_bool = True

        if not bottom_bool:
            ax[site].set_xticklabels([])
        if not left_bool:
            ax[site].set_yticklabels([])

        # fig.suptitle('Night time')
    
    fig.supylabel('Probability density function (-)',x=0.09)
    fig.supxlabel('Standardised daytime F$_c$+S$_c$ (-)',y=0.06)
    fig.savefig(figfold+'/fig14.png', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/fig14.pdf', dpi=200,bbox_inches='tight')
        



def plot_omega_schematic(data,sites,clrs,sitetext):
    sigma_w = np.linspace(0,2,1000)
    wecrit = np.linspace(-4,0,1000)
    text_loc = [1.4,-1.4/OMEGAth]


    X, Y = np.meshgrid(sigma_w, wecrit)
    Z = X/np.abs(Y)

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize2[0]*cm, figsize2[1]*cm,forward=True)
    p1 = ax.pcolormesh(X, Y, Z,vmin=0,vmax=OMEGAth*2,cmap=mpl.cm.BrBG,rasterized=True)

    contour_legend = []
    contour_text = []
    for site in sites:
        m1 = data[site]['W_SIGMA'].copy()
        m2 = data[site]['WEC'].copy()
        m2[m2==0] = np.nan
        indx = (~m1.isnull()) & (~m2.isnull())
        # indx = (~m1.isnull()) & (~m2.isnull()) & (data[site]['SW_IN']<10)
        m1 = m1[indx]
        m2 = m2[indx]
        xmin = 0
        xmax = m1.max()
        ymin = m2.min()
        ymax = 0
        X1, Y1 = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]
        positions = np.vstack([X1.ravel(), Y1.ravel()])
        values = np.vstack([m1, m2])
        kernel = scipy.stats.gaussian_kde(values)
        Z1 = np.reshape(kernel(positions).T, X1.shape)

    
        cntr1 = ax.contour(X1, Y1, Z1, levels=4,linestyles='-',colors=[clrs[site]], linewidths=1)
        h1,_ = cntr1.legend_elements()
        contour_legend.append(h1[0])
        if site in sitetext.keys():
            contour_text.append(sitetext[site])
        else:
            contour_text.append(site)

    ax.contour(X,Y,Z, levels=[OMEGAth],linestyles='dashed',colors='k')

    ax.text(text_loc[0],text_loc[1], '$\Omega=\Omega_{th}$',
            rotation=np.arctan(text_loc[1]/text_loc[0])*180/np.pi, rotation_mode='anchor',
            transform_rotates_text=True, ha='right', va='bottom')
    
    # ax.legend(contour_legend,sites, loc='lower right')
    ax.legend(contour_legend,contour_text, loc='lower right',fontsize=8)

    ax.set_xlim([0,1.5])
    ax.set_ylim([-3,0])
    ax.set_ylabel('$w_{e,crit}$ (m/s)')
    ax.set_xlabel('$\sigma_w$ (m/s)')
    cb = fig.colorbar(p1,extend='max')
    cax = cb.ax
    cax.hlines(OMEGAth, 0, 1, colors = 'k', linestyles = '--')
    cb.set_label('$\Omega$ (-)')
    # plt.show()
    fig.savefig(figfold+'/omega_schematic_ver2.pdf', dpi=200,bbox_inches='tight')
    fig.savefig(figfold+'/omega_schematic_ver2.png', dpi=200,bbox_inches='tight')

def define_colors(z):

    val = np.random.rand(1,1)[0][0]*0.05
    if z<10:
        clr = list(cmp.batlow(0.1))
    elif (z>=10) & (z<30):
        clr = list(cmp.batlow(0.5))
    elif (z>=30):
        clr = list(cmp.batlow(0.7))

    clr = list(clr[0:3]+val)
    clr.append(1.0)
    return tuple(clr)


def main():
    data_all = dict()
    USTARth = dict()
    USTARth['FC'] = dict()
    USTARth['NEE'] = dict()
    clrs = dict()
    for site in config['sites']:

        filein = os.path.join(os.getcwd(),'data',site+'.csv')
        datain = iotools.read_processed_data(filein)

        for var in ['FC','NEE']:
            USTARth[var][site] = iotools.load_ONEflux_threshold(os.path.join(os.getcwd(),'data'),site,var)

        datain['OMEGA_FLAG'] = 0
        datain.loc[(pd.isnull(datain['USTAR'])) | (pd.isnull(datain['OMEGA'])),'OMEGA_FLAG'] = np.nan
        datain.loc[datain['OMEGA']<OMEGAth,'OMEGA_FLAG'] = 1
        if not datain.empty:
            data_all[site] = datain
            clrs[site] = utils.define_colors(datain.loc[datain.index[0],'z'])
        del datain


    sID = 'DE-HoH'
    USTARth_example = USTARth['FC'][sID].iloc[4].mean()
    fig2(data_all,sID,USTARth_example)
    fig3_4_5(data_all,clrs)
    fig7(data_all,USTARth,clrs)
    fig8(data_all,clrs)
    fig9(data_all,clrs)
    fig10(data_all,clrs)
    site = ['FI-Hyy','SERC']
    fig11(data_all,clrs,site)
    fig12(data_all,clrs)
    fig13_14(data_all,clrs)
    forest_site = 'DE-HoH'
    grass_site = 'KONA'
    figC15(data_all,clrs,forest_site,grass_site)


if __name__ == "__main__":
    main()