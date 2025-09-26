#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as fm
from matplotlib import rcParams
import scipy
font_path = 'arial.ttf'
custom_font = fm.FontProperties(fname=font_path)
fm.fontManager.addfont(font_path)
rcParams['font.family'] = custom_font.get_name()
s_font = 18


#% SDDP&Proposed
from matplotlib.lines import Line2D

make = 'Sample' #'Time'  #'Sample'

if make == 'Time':
    SDDP = np.array([[150.0899349,305.3327509,455.601675,2425.750224,3811.04621],[51.95043137,106.3876838,158.936052,255.4874618,332.4472764]])
    Pro = np.array([[0.5955181,1.5590226,2.8049279,3.7569818,4.6367876],[51.95204986,106.387799,158.9371368,255.6275064,332.9608456]])
    labels = ['2h','4h','6h','8h','10h']
else:
    SDDP = np.array([[538.3051696,1525.248976,790.3121741,1382.501419,1390.807388],[79.68735508,79.68735508,79.68735508,79.68735508,79.68735508]])
    Pro = np.array([[1.0436567,0.9434197,0.956537,0.9519198,0.9695071],[79.6881922,79.6881922,79.6881922,79.6881922,79.6881922]])
    labels = ['20','40','60','80','100']

fig = plt.figure(figsize=(10, 6))
ax = plt.gca()
ax.set_xscale('log')
for ii in range(5):
    ax.scatter(Pro[0,ii], Pro[1,ii], color='black', marker=['o','^','v','*','s'][ii] , s=100)
    ax.scatter(SDDP[0,ii], SDDP[1,ii], color='black', marker=['o','^','v','*','s'][ii] ,s=100)    
ax.tick_params(labelsize=s_font)
ax.set_xlabel('Time (s)', fontsize=s_font+2)
ax.set_ylabel('Objective value ($)', fontsize=s_font+2)
if make == 'Time':
    ymin,ymax = 0,450
else:
    ymin,ymax = 75,85
ax.set_ylim(ymin,ymax)
ax.set_xlim(1e-1,1e4)
ax.grid(True,axis = 'y', which = "both", ls="--")
# main_legend = ax.legend(['Proposed','SDDP'], fontsize=s_font, loc='upper left')

custom_lines = [
    Line2D([0], [0], color='black', marker='o', linestyle='None', markersize=8, label=labels[0]),
    Line2D([0], [0], color='black', marker='^', linestyle='None', markersize=8, label=labels[1]),
    Line2D([0], [0], color='black', marker='v', linestyle='None', markersize=8, label=labels[2]),
    Line2D([0], [0], color='black', marker='*', linestyle='None', markersize=8, label=labels[3]),
    Line2D([0], [0], color='black', marker='s', linestyle='None', markersize=8, label=labels[4]),
]

# ax.add_artist(main_legend)
custom_legend = ax.legend(handles=custom_lines, loc='upper center', fontsize=s_font, ncol = 5)

ax.axvspan(0, 10, facecolor='lightpink', alpha=0.3)
ax.axvspan(10, 1e4, facecolor='lightblue', alpha=0.3)

ax.text(0.6, ymin+0.8*(ymax-ymin),'Proposed',fontsize=s_font,color='black')
ax.text(3e2, ymin+0.8*(ymax-ymin),'SDDP',fontsize=s_font,color='black')

plt.show()
fig.savefig('Figs/SDDP_'+make+'.pdf',format='pdf',dpi=600,bbox_inches='tight')


#%% price
data_in = scipy.io.loadmat('Data.mat')
Pri_rt_Jan,Pri_da_Jan,Pri_rt_Feb,Pri_da_Feb = data_in['Pri_rt_Jan'],data_in['Pri_da_Jan'],data_in['Pri_rt_Feb'],data_in['Pri_da_Feb']

X = [x for x in range(289) for _ in range(2)][1:-1]
dd = 6
Pri_da_Feb_list = [Pri_da_Feb[288*dd+x,0] for x in range(288) for _ in range(2)]
Pri_rt_Feb_list = [Pri_rt_Feb[288*dd+x,0] for x in range(288) for _ in range(2)]

diff = Pri_rt_Jan.reshape(31,288)-Pri_da_Jan.reshape(31,288)

fig = plt.figure(figsize=(12, 4))
plt.plot(X, Pri_da_Feb_list, label='DAP', color='steelblue', linestyle='--', linewidth=1.5)
plt.plot(X, Pri_rt_Feb_list, label='RTP', color='darkslateblue', linewidth=1.5)
for i in range(24):
    y = Pri_da_Feb[288*dd+i*12,0]+diff[:,i*12]
    plt.boxplot(y,positions=[i*12],showfliers=False,widths = 6)
plt.ylim(0,200)
plt.xlim(-5,288)
plt.xticks([12*0,12*6,12*12,12*18,12*24],['0:00', '6:00', '12:00', '18:00', '24:00'])
plt.tick_params(labelsize=s_font)
plt.xlabel('Hour', fontsize=s_font+2)
plt.legend(fontsize=s_font)
plt.ylabel('Price ($/MWh)', fontsize=s_font+2)
plt.show()
fig.savefig('Figs/price.pdf',format='pdf',dpi=600,bbox_inches='tight')


#%% decision results

data = np.random.rand(3, 4, 10)

data[0,:,:] = np.array([
[700.0235304, 543.180105, 1437.921634, 2242.953387, 1486.135144, 926.296417, 605.7651393, 521.7268733, 331.6392508, 274.0586095],
[713.5050048, 543.5366298, 1470.416783, 2249.778892, 1893.029537, 932.5753531, 619.5747301, 524.7032998, 332.82613, 280.9962702],
[713.5050048, 543.5366298, 1469.962214, 2246.275872, 1860.622457, 932.5753531, 619.5747301, 524.7032998, 332.82613, 280.9962702],
[713.5050048, 543.5366298, 1469.962214, 2246.275872, 1848.309085, 932.5753531, 619.5747301, 524.7032998, 332.82613, 280.9962702]])

data[1,:,:] = np.array([
[686.7562316, 531.3868579, 1436.048378, 2221.487717, 1456.732288, 910.0139076, 607.163557, 516.3225455, 326.6891666, 258.0343574],
[710.2524266, 546.8407745, 1573.583082, 2254.034549, 1715.705224, 955.4950919, 618.3962439, 526.4019736, 330.9717504, 266.3019134],
[713.2011023, 551.0781111, 1583.234409, 2280.582853, 1700.495177, 954.4068845, 620.2692766, 525.6203197, 330.195015, 265.5156845],
[713.9952204, 551.5201138, 1584.533236, 2284.619586, 1698.02321, 954.2140321, 620.6134205, 525.4199273, 330.0613392, 265.3498923]
])

data[2,:,:] = np.array([
[660.9896006, 527.6117533, 1351.247864, 2200.92024, 1334.746429, 876.7277514, 589.1147145, 507.836571, 325.7756139, 255.3701598],
[700.7949174, 573.4593395, 1535.425601, 2355.399595, 1864.108732, 956.2731058, 607.4013243, 523.5816665, 326.3157171, 268.7480734],
[712.4528031, 596.9323066, 1545.019495, 2381.355085, 1840.894705, 967.2528047, 609.9057681, 523.5816665, 326.3157171, 268.7480734],
[713.3022525, 598.7146798, 1542.979313, 2385.382584, 1840.894705, 967.2528047, 610.7782037, 523.5816665, 326.3157171, 268.7480734]
])

groups = [f'Day {i+1}' for i in range(10)]
bars = ['10 Samples', '50 Samples', '100 Samples']
fig, axes = plt.subplots(3, 1, figsize=(14.4, 6), constrained_layout=True)
for i, ax in enumerate(axes):
    subplot_data = data[i]
    x = np.arange(len(groups))
    width = 0.2
    positions = [x - width, x, x + width]
    for j in range(3):
        ax.bar(positions[j], subplot_data[1+j,:], width, label = bars[j],color=['#97B2DE',"#1A488E",'#092147'][j])
    ax.plot(positions[1], subplot_data[0,:], 'o--' ,color='black',label='Ground-truth')
    ax.set_title(['Linear','Quadratic','Piecewise'][i], fontsize=s_font+2, y=0.8)
    ax.set_ylabel('Objective ($)', fontsize=s_font+2)
    if i == 2:
        ax.set_xticks(x)
        ax.set_xticklabels(groups, fontsize=s_font+2)
    else:
        ax.set_xticklabels([])
    ax.tick_params(labelsize=s_font)
    ax.legend(loc = 'upper right', fontsize=s_font-2)
    # for pos_group in positions:
    #     for pos, val in zip(pos_group, subplot_data[:, j]):
    #         ax.text(pos, val + 2, f'{val:.1f}', ha='center')
plt.show()
fig.savefig('Figs/Decision.pdf',format='pdf',dpi=600,bbox_inches='tight')