import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

import iotools

sns.set_theme(style='ticks',context='paper')
# sns.set_palette(cmp.batlow(np.linspace(0,1,500)))
cm = 1/2.54

var='a'
al = []
al = [(chr(ord(var)+i)) for i in range(26)]

OMEGAth = 0.66
OMEGAth_std = 0.06
OMEGAth_skewness = 0.52
OMEGAth_quantiles = [0.59,0.61,0.63,0.64,0.65,0.67,0.69,0.71,0.73]
config_file = 'config.yml'
config = iotools.load_config(config_file)   

figfold = os.path.join(os.getcwd(),'figs')
# figsize = [18,12]
figsize = [18,8]
figsize2 = [12,8]
figsize3 = [18,8]
text_loc = [0.7,-0.7/OMEGAth]

# sigma_w = np.linspace(1e-2,1,100)
# wecrit = np.linspace(-4,-1e-1,100)
sigma_w = np.linspace(0,2,1000)
wecrit = np.linspace(-4,0,1000)


X, Y = np.meshgrid(sigma_w, wecrit)
Z = X/np.abs(Y)


fig, ax = plt.subplots()
fig.set_size_inches(figsize2[0]*cm, figsize2[1]*cm,forward=True)
# ax.pcolormesh(X, Y, Z,vmin=0,vmax=OMEGAth*2,cmap=mpl.cm.RdBu_r)
p1 = ax.pcolormesh(X, Y, Z,vmin=0,vmax=OMEGAth*2,cmap=mpl.cm.BrBG,rasterized=True)
# p1 = ax.contourf(X, Y, Z, levels=np.linspace(0,OMEGAth*2,10),cmap=mpl.cm.BrBG)
# ax.contour(X,Y,Z, levels=[0.1,0.2,0.30,0.4,0.5,OMEGAth,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2])
p2 = ax.contour(X,Y,Z, levels=[OMEGAth],linestyles='dashed',colors='k')

th2 = ax.text(text_loc[0],text_loc[1], '$\Omega=\Omega_{th}$',
              rotation=np.arctan(text_loc[1]/text_loc[0])*180/np.pi, rotation_mode='anchor',
              transform_rotates_text=True, ha='right', va='bottom')
# th2 = ax.text(text_loc[0],text_loc[1], 'text rotated correctly')
# fmt = {}
# strs = ['\Omega=\Omega$_{th}$']
# strs = ['asd']
# for l, s in zip(p2.levels, strs):
#     fmt[l] = s

# Label every other level using strings
# ax.clabel(p2, p2.levels, inline=True, fmt=fmt)

ax.set_xlim([0,1.5])
ax.set_ylim([-3,0])
ax.set_ylabel('$w_{e,crit}$ (m/s)')
ax.set_xlabel('$\sigma_w$ (m/s)')
cb = fig.colorbar(p1,extend='max')
cax = cb.ax
cax.hlines(OMEGAth, 0, 1, colors = 'k', linestyles = '--')
cb.set_label('$\Omega$ (-)')


# plt.show()
plt.text(.95, .95, 'Coupled', ha='right', va='top', transform=ax.transAxes,backgroundcolor='w')
plt.text(.05, .05, 'Decoupled', ha='left', va='bottom', transform=ax.transAxes,backgroundcolor='w')
fig.savefig(figfold+'/omega_schematic.pdf', dpi=200,bbox_inches='tight')