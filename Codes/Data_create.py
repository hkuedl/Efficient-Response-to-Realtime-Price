#%%
import numpy as np
import math
import pandas as pd
import os
import random
import copy
import scipy.stats as stats
def set_seed(seed):
    # seed init
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
set_seed(20)
import scipy.io

data_in = np.loadtxt("90zone_15min.csv",delimiter=",",skiprows=1,usecols=range(1,5*90+2))
T_rad0 = data_in[:,1+90*1:1+90*2]/1000
T_o0 = data_in[:,0].reshape(-1,1)
T_occ0 = data_in[:,1:1+90*1].copy()/1000

N_days = 10
delta = np.zeros((288,N_days))
for day_i in range(N_days):
    delta_0 = 0.3 + 0.3*T_occ0[day_i*96:(day_i+1)*96,0] + 0.3*T_rad0[day_i*96:(day_i+1)*96,0] + 0.05*T_o0[day_i*96:(day_i+1)*96,0]
    original_indices = np.arange(96)
    new_indices = np.linspace(0, 95, 288)
    delta_i = np.interp(new_indices, original_indices, delta_0[:]).reshape(-1,1)
    delta_i += np.random.uniform(low=-0.05, high=0.05, size=(12*24, 1))
    delta[:,day_i:day_i+1] = delta_i

#Jan: 1-31
i_day = ['01','02','03','04','05','06','07','08','09'] + [str(i) for i in range(10,32)]
Pri_rt_Jan = np.array(pd.read_csv('RT_2023/data_1.csv').iloc[:,10]).reshape(-1,1)
Pri_da_Jan_read = [pd.read_csv('DA_Jan_2023/202301'+i_day[i]+'damlbmp_zone.csv') for i in range(31)]
Pri_da_Jan_all = np.zeros((24*31,1))
for i in range(31):
    Pri_da_Jan_all[24*i:24*(i+1),0] = Pri_da_Jan_read[i][Pri_da_Jan_read[i]['Name'] == 'N.Y.C.'].iloc[:,3]
Pri_da_Jan = np.repeat(Pri_da_Jan_all, repeats=12, axis=0)

#Feb: 1-10
Pri_rt_Feb = np.array(pd.read_csv('RT_2023/data_2.csv').iloc[0*288:N_days*288,10]).reshape(-1,1)  #10 is N.Y.C zone
Pri_da_Feb_read = [pd.read_csv('DA_Feb_2023/202302'+i_day[i]+'damlbmp_zone.csv') for i in range(0,N_days)]
Pri_da_Feb_all = np.zeros((24*N_days,1))
for i in range(N_days):
    Pri_da_Feb_all[24*i:24*(i+1),0] = Pri_da_Feb_read[i][Pri_da_Feb_read[i]['Name'] == 'N.Y.C.'].iloc[:,3]
Pri_da_Feb = np.repeat(Pri_da_Feb_all, repeats=12, axis=0)

# scipy.io.savemat('Data.mat', {'delta': delta, 'Pri_rt_Jan': Pri_rt_Jan, 'Pri_da_Jan': Pri_da_Jan,\
#                               'Pri_rt_Feb': Pri_rt_Feb, 'Pri_da_Feb': Pri_da_Feb})

## generate the price distribution
def calculate_real_price_distribution(day_ahead_prices_30d, real_time_prices_30d,sam_i_num):
    hourly_diff = real_time_prices_30d - day_ahead_prices_30d
    percentiles = np.linspace(0, 1, sam_i_num+2)[1:-1]
    hourly_percentiles = np.zeros((288,sam_i_num))

    for hour in range(2):
        hour_diff = hourly_diff[:,hour]
        hourly_percentiles[hour,:] = np.percentile(hour_diff, percentiles * 100)
    return hourly_percentiles

sam_num = [10,20,30,40,50,60,70,80,90,100]
quantiles_all = {}
for sam_i in range(len(sam_num)):
    sam_i_num = sam_num[sam_i]
    hourly_diff_percentiles = calculate_real_price_distribution(Pri_da_Jan.reshape(31,288), Pri_rt_Jan.reshape(31,288),sam_i_num)
    quantiles_x = np.linspace(0, 1, sam_i_num+2)
    quantiles_all['x_'+str(sam_i_num)] = quantiles_x
    for dd in range(N_days):
        quantiles = np.zeros((288,sam_i_num))
        da_price = Pri_da_Feb[288*dd:288*(dd+1),:]
        quantiles = da_price + hourly_diff_percentiles
        quantiles_all['y_'+str(dd)+'_'+str(sam_i_num)] = quantiles
#scipy.io.savemat('Data_quan.mat', quantiles_all)




a1 = Pri_rt_Jan.reshape(31,288)-Pri_da_Jan.reshape(31,288)
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 5))
plt.plot(a1)
plt.ylim(-100,500)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(a1)
plt.ylim(-100,500)
plt.show()

#%%


#% Continue to formulate the price forecasting distribution
for N_scale in [20]:
    scale_random = N_scale*np.ones((12*24*N_days,1))
    center_random = np.random.uniform(-0.0, 0.0, (12*24*N_days,1))
    center_rt = Pri_rt_Feb + center_random
    
    sam_num = [10,20,30,40,50,60,70,80,90,100]
    quantiles_all = {}
    for sam_i in range(len(sam_num)):
        sam_i_num = sam_num[sam_i]
        quantiles_x = np.linspace(0, 1, sam_i_num+2)
        quantiles_all['x_'+str(sam_i_num)] = quantiles_x
        for dd in range(N_days):
            quantiles = np.zeros((288,sam_i_num))
            for ii in range(288):
                normal_dist = stats.norm(loc=center_rt[dd*288+ii,0], scale = scale_random[dd*288+ii,0])
                quantiles[ii,:] = normal_dist.ppf([i for i in quantiles_x[1:-1]])
            quantiles_all['y_'+str(dd)+'_'+str(sam_i_num)] = quantiles

    scipy.io.savemat('Data_quan_'+str(N_scale)+'.mat', quantiles_all)