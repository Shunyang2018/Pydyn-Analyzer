

import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pandas as pd
import time
import matplotlib.collections as collections
import file_proc_pkg.func as func
import file_proc_pkg.class_traj as t
import scipy.signal as signal

###### Read configuration ######
start_time = time.perf_counter()
        
#### Traj Set-up ####

direction = 'only'  # 'both' for traj of both directions
                    # 'only' for traj of only one direction
                    
#### what to plot ####
                    
plot_potential_energy_surface = 'on'
plot_timing_histogram = "off"

#### plot set-up ####

plot_size = (8,6)
diagram_size = (7.5,3)
timing_histogram = "off"
is_3d_plot = 'off'
color_bar_for_pes = 'on' # more customization can be manually done in the code below
traj_anal = True
pes_data_file = "dataplot.csv"

traj_plot_setup = 'all'
traj_plot_amount = 9

###### Plot PES #######
fig = plt.figure(figsize = plot_size, dpi=100)
if plot_potential_energy_surface == 'on':
    
    print('PES data Found in',pes_data_file)
    data = pd.read_csv(open(pes_data_file), delimiter=",")
    data.values
    data = np.array(data)
    # Make data.
    
    y_range = np.array([2.505, 1.905])
    x_range = np.array([2.5997290915,1.9])
    
    dim_row,dim_col = np.shape(data)
    Y = np.arange(x_range[0], x_range[1], (x_range[1]-x_range[0])/dim_row)
    X = np.arange(y_range[0],y_range[1], (y_range[1]-y_range[0])/dim_col)
    
    X, Y = np.meshgrid(X, Y)
    # Plot the surface.
    grad = np.gradient(data)
    surf = plt.contourf(X, Y, data, offset=-698.2, cmap=cm.Spectral_r)
    
    # Add a color bar which maps values to colors.
    if color_bar_for_pes == 'on': 
        fig.colorbar(surf, shrink=0.5, aspect=10)
    plt.title("Potential Energy Surface")
    plt.xlabel("Bond Length C1-C2")
    plt.ylabel("Bond Length C1-C3")
    
    
###### Main Loop ########
#---------------- read all trajs in folders ----------------#
    
all_traj = [] # store all the trajs
time_list = []
trajname_list = []
index = 0
if is_3d_plot == 'on':
    
    ax = fig.gca(projection='3d')

anal_rules = func.read_anal_file()

directory_info_list = list(os.walk(os.getcwd()))

for i in range(0,len(directory_info_list[0][1])+2):
    
    for filename in directory_info_list[i][2]:
        
        if filename[0:4] == 'traj' and direction == "both":
            
        
            dyn_coord = func.get_raw_dyn_downhill(directory_info_list[i][0]+'/'+filename)
            traj = t.Downhill_trajectory(dyn_coord[0],dyn_coord[1],filename,index,dyn_coord[2],dyn_coord[3])

            if traj.complete == True:
                
                all_traj.append(traj)
                time_list.append(traj.time)
                func.anal_traj_file(anal_rules[0],traj)
                func.assign_traj(traj,anal_rules[1],anal_rules[2])
                trajname_list.append(traj.end_point_name)
                
            
        if filename[0:4] == 'traj' and direction == "only":
            
            dyn_coord = func.get_raw_dyn_uphill(directory_info_list[i][0]+'/'+filename)
            traj = t.Uphill_trajectory(dyn_coord[0],dyn_coord[1],filename,index,dyn_coord[2],True)
            all_traj.append(traj)
            time_list.append(traj.time)            
            func.anal_traj_file(anal_rules[0],traj)
            func.assign_traj(traj,anal_rules[1],anal_rules[2])
            trajname_list.append(traj.end_point_name)

            index += 1    

#---------------- Filtering out traj that need to be plotted ----------------#

if traj_plot_setup == 'random':
    
    traj_index = np.random.randint(0,len(all_traj))
    all_traj[traj_index].plotting_status = True
    
else:
    
    traj_plot_amount = len(all_traj)
    
peak_list = []
   
for traj in all_traj:
    
    traj.traj_filter(traj_plot_setup,len(all_traj),traj_plot_amount)
    
    if traj.plotting_status == True :
        
#        plt.plot(traj.distance(3,2),traj.distance(2,16),linewidth=0.7)
#        plt.plot(traj.distance(7,5),traj.distance(3,5))
#        peak_list.append(len(signal.argrelextrema(np.array((traj.distance(7,5))),np.greater)[0]))
        pass
#print(peak_list)
#print(np.average(peak_list))


#---------------- Histogram for  ----------------#
###### Plot all Traj ######
            
#traj_to_plot = specify_traj_for_plot(traj_plot,len(all_traj))
me_list = []
et_list = []
traj_plotted = 0
time_list_me = []
time_list_h = []
theta_list_me,theta_list_h = [],[]
fig,ax = plt.subplots()
i = 0
atom = 16

'''
for traj in all_traj:
            
    if traj.structname in traj_to_plot[0] and traj.index in traj_to_plot[1] and traj_plotted < traj_to_plot_num :
        
#        me_list.append(traj.func_grp_ke([26,25,27,4])[0])
#        et_list.append(traj.func_grp_ke([3,31,32,33,34,35,36])[0])
#        plt.ylim((0,0.02))
#        plt.xlim((0,1))
#        plt.hist(me_list,1,density=True, histtype='bar', facecolor='dodgerblue', rwidth=1)
#        plt.hist(me_list,5,density=True, histtype='bar', facecolor='orange', rwidth=1)
#        plt.plot(range(traj.time-1),traj.func_grp_ke([26,25,27,4]))
#        plt.plot(traj.distance(3,2),traj.distance(16,2),linewidth=0.4)
        if traj_plotted in [4]:
            
            if traj.structname == 'Me_Product':
                
    #            atom = 16
                time_list_me.append(traj.time)
                ax.plot(range(traj.time-1),traj.force_along_direction(atom))
                collection = collections.BrokenBarHCollection.span_where(range(traj.time-1), ymin=0, ymax=0.1, where=np.array(traj.force_along_direction(atom)) > 0, facecolor='r', alpha=0.5)
                ax.add_collection(collection)
                collection = collections.BrokenBarHCollection.span_where(range(traj.time-1), ymin=-0.1, ymax=0, where=np.array(traj.force_along_direction(atom)) < 0, facecolor='g', alpha=0.5)
                ax.add_collection(collection)
                print('Me_Product')
                
            if traj.structname == 'H_Product':
                
                
                time_list_h.append(traj.time)
    #            plt.plot(range(traj.time-1),traj.theta(3))
                ax.plot(range(traj.time-1),traj.force_along_direction(atom))
                collection = collections.BrokenBarHCollection.span_where(range(traj.time-1), ymin=0, ymax=0.1, where=np.array(traj.force_along_direction(atom)) > 0, facecolor='r', alpha=0.5)
                ax.add_collection(collection)
                collection = collections.BrokenBarHCollection.span_where(range(traj.time-1), ymin=-0.1, ymax=0, where=np.array(traj.force_along_direction(atom)) < 0, facecolor='g', alpha=0.5)
                ax.add_collection(collection)
                print('H_Product')
                
        traj_plotted += 1 
'''       

###### Analyzing all Traj ######
trajtype_set = []
for each_set in trajname_list:
    
    if each_set not in trajtype_set:
        trajtype_set.append(each_set)
                   
print('Total {} Trajectories were found in folder {}'.format(len(all_traj),os.getcwd()))

for eachtype in trajtype_set:
    
    number_of_traj = trajname_list.count(eachtype)
    
    if len(eachtype) == 1:
    
        eachtype = list(eachtype)
        print( 'There are total',number_of_traj, eachtype[0],'to',eachtype[0],'Trajectories, ',\
              float(number_of_traj*100/len(all_traj)),'% of total')
        
    else:
        
        eachtype = list(eachtype)
        print( 'There are total',number_of_traj, eachtype[0],'to',eachtype[1],'Trajectories, ',\
              float(number_of_traj*100/len(all_traj)),'% of total')

          
########## Plot ###########

if plot_timing_histogram == "on":
    
    plt.figure(figsize=diagram_size, dpi=100)
    plt.ylabel('Probability density')
    plt.xlabel('Time/fs')
    plt.title('Probability Density Distribution of Ethyl Migration in 6')
    plt.hist(time_list_h,40,density=True, histtype='bar', facecolor='dodgerblue', rwidth=0.6)
    plt.hist(time_list_me,40,density=True, histtype='bar', facecolor='g', rwidth=0.5)
    plt.xlim((0, 300))
    plt.ylim((0, 0.1))
    print('mean',np.mean(time_list))
    print('stand deviation', np.std(time_list))

end_time = time.perf_counter()

###### Momentum Analysis ######



print("It takes ", end_time-start_time,'sec to process')
print("Normal Termination")

