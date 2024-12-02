import numpy as np
import pandas as pd

import utils
import constants as c


# case #1: short tower above short vegetation (e.g. grassland)
# EC measurement height
z = 4
# temperature profile
Tprofile = pd.DataFrame(data={'TPOT_1':1,'TPOT_2':2,'TPOT_3':3},index=[0])
# temperature measurement heights
zT = {'TPOT_1':1,'TPOT_2':2,'TPOT_3':3}
# canopy height
canz = 0
# wind speed at canopy height
Uh = pd.Series(np.nan)
# plant area index
PAI = pd.Series(np.nan)
# standard deviation of vertical wind component
W_SIGMA = pd.Series(0.6)

omega,wec,reftheta = utils.calculate_omega(W_SIGMA,Tprofile+c.NT,Uh,PAI,canz,z,zT)
print('Omega for stably stratified flow above short canopy with short tower:'+str(omega.loc[0]))



# case #2: tall tower above short vegetation (e.g. grassland)
# EC measurement height
z = 40
# temperature profile
Tprofile = pd.DataFrame(data={'TPOT_1':1,'TPOT_2':2,'TPOT_3':3},index=[0])
# temperature measurement heights
zT = {'TPOT_1':10,'TPOT_2':20,'TPOT_3':30}
# canopy height
canz = 0
# wind speed at canopy height
Uh = pd.Series(np.nan)
# plant area index
PAI = pd.Series(np.nan)
# standard deviation of vertical wind component
W_SIGMA = pd.Series(0.6)

omega,wec,reftheta = utils.calculate_omega(W_SIGMA,Tprofile+c.NT,Uh,PAI,canz,z,zT)
print('Omega for stably stratified flow above short canopy with tall tower:'+str(omega.loc[0]))





# case #3: tall tower above sparse forest
# EC measurement height
z = 40
# temperature profile
Tprofile = pd.DataFrame(data={'TPOT_1':1,'TPOT_2':2,'TPOT_3':3},index=[0])
# temperature measurement heights
zT = {'TPOT_1':10,'TPOT_2':20,'TPOT_3':30}
# canopy height
canz = 20
# wind speed at canopy height
Uh = pd.Series(2)
# plant area index
PAI = pd.Series(1)
# standard deviation of vertical wind component
W_SIGMA = pd.Series(0.6)

omega,wec,reftheta = utils.calculate_omega(W_SIGMA,Tprofile+c.NT,Uh,PAI,canz,z,zT)
print('Omega for stably stratified flow above sparse forest:'+str(omega.loc[0]))




# case #4: tall tower above dense forest
# EC measurement height
z = 40
# temperature profile
Tprofile = pd.DataFrame(data={'TPOT_1':1,'TPOT_2':2,'TPOT_3':3},index=[0])
# temperature measurement heights
zT = {'TPOT_1':10,'TPOT_2':20,'TPOT_3':30}
# canopy height
canz = 20
# wind speed at canopy height
Uh = pd.Series(2)
# plant area index
PAI = pd.Series(5)
# standard deviation of vertical wind component
W_SIGMA = pd.Series(0.6)

omega,wec,reftheta = utils.calculate_omega(W_SIGMA,Tprofile+c.NT,Uh,PAI,canz,z,zT)
print('Omega for stably stratified flow above dense forest:'+str(omega.loc[0]))



# case #5: tall tower above dense forest
# EC measurement height
z = 40
# temperature profile
Tprofile = pd.DataFrame(data={'TPOT_1':15,'TPOT_2':15,'TPOT_3':15},index=[0])
# temperature measurement heights
zT = {'TPOT_1':10,'TPOT_2':20,'TPOT_3':30}
# canopy height
canz = 20
# wind speed at canopy height
Uh = pd.Series(2)
# plant area index
PAI = pd.Series(6)
# standard deviation of vertical wind component
W_SIGMA = pd.Series(0.6)

omega,wec,reftheta = utils.calculate_omega(W_SIGMA,Tprofile+c.NT,Uh,PAI,canz,z,zT)
print('Omega for neutrally stratified flow above dense forest:'+str(omega.loc[0]))