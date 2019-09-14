#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 07:47:57 2019

@author: zhtfeng
"""
######## the class traj contains all the atrributes of a trajctory,
######## and with functions we can pull out "any" property needed 

import numpy as np
        
        
class Trajectory:
    
    def __init__(self,coord_list,energy_list,name,index,atom_list,completeness):
        
        self.coord = coord_list
        self.energy = energy_list
        self.index = index
        self.name = name
        self.time = len(self.energy)
        self.attribute = [1]
        self.end_point_name = set()
        self.atom = atom_list
        self.starting_attribute = [1]
        self.plotting_status = False
        self.complete = completeness
        self.recrossing = False
        self.end_point_type = set()
        
    def distance(self,a,b): 
        
        dist_list = []
        
        for point in range(self.time):
            
            dist = self.dist(a,b,point)    
            dist_list.append(dist)
        
        return dist_list
    
    def angle(self,a,b,c):
        
        angle_list = []
        
        for point in range(self.time):
            
            ang = self.ang(a,b,c,point)
            angle_list.append(ang)
            
        return angle_list
            
    def dihedral(self,a,b,c,d):
        
        dihedral_list = []
        
        for point in range(self.time):
            
            dih = self.dih(a,b,c,d,point)
            dihedral_list.append(dih)
            
        return dihedral_list
            
    
    def initialize_attribute(self,full):
        
        self.attribute = self.attribute*full
        self.starting_attribute = self.starting_attribute*full
        
        return self.attribute,self.starting_attribute
    
    def dist(self,a,b,point):
        
        a -= 1
        b -= 1
        cart = self.coord[int(point)]
        dist = np.sqrt((cart[a,0]-cart[b,0])**2 + (cart[a,1]-cart[b,1])**2 + \
                           (cart[a,2]-cart[b,2])**2)    
        
        return dist
    
    def ang(self,a,b,c,point):
        
        a -= 1
        b -= 1
        c -= 1
        cart = self.coord[int(point)]
        vec_a = np.array([(cart[a,0]-cart[b,0]),(cart[a,1]-cart[b,1]),(cart[a,2]-cart[b,2])])
        vec_b = np.array([(cart[c,0]-cart[b,0]),(cart[c,1]-cart[b,1]),(cart[c,2]-cart[b,2])])
        cos_theta = vec_a.dot(vec_b)/(np.sqrt(vec_a.dot(vec_a)*vec_b.dot(vec_b)))
        theta = np.arccos(cos_theta)
        theta = 180*theta/np.pi

        return theta
    
    def dih(self,a,b,c,d,point):
        
        a -= 1
        b -= 1
        c -= 1
        d -= 1
        cart = self.coord[int(point)]
        vec_a_1 = np.array([(cart[a,0]-cart[b,0]),(cart[a,1]-cart[b,1]),(cart[a,2]-cart[b,2])])
        vec_a_2 = np.array([(cart[c,0]-cart[b,0]),(cart[c,1]-cart[b,1]),(cart[c,2]-cart[b,2])])
        vec_b_1 = np.array([(cart[b,0]-cart[c,0]),(cart[b,1]-cart[c,1]),(cart[b,2]-cart[c,2])])
        vec_b_2 = np.array([(cart[d,0]-cart[c,0]),(cart[d,1]-cart[c,1]),(cart[d,2]-cart[c,2])])
        norm_a = np.cross(vec_a_1,vec_a_2)
        norm_b = np.cross(vec_b_1,vec_b_2)
        cos_theta = norm_a.dot(norm_b)/(np.sqrt(norm_a.dot(norm_a)*norm_b.dot(norm_b)))
        theta = np.arccos(cos_theta)
        theta = 180*theta/np.pi
        
        return theta
    
    def generate_mass_list(self):
        
        mass_list = np.zeros((len(self.atom),1))
        
        for i,each_atom in enumerate(self.atom):
            
            mass_list[i] = mass_dict[each_atom]
            
        return mass_list
    
    def velocity(self,n):
        
        n -= 1
        velocity_list = []
        vel_n = []
        kinetic_energy_list = []
        vel_vec_n = []
        
        for pts in range(self.time-1):
            
            velocity = np.zeros((len(self.atom),3))
            for atm_num,atom in enumerate(self.coord[pts]):
                
                next_pts = self.coord[pts+1]
                velocity[atm_num,0] = next_pts[atm_num,0] - atom[0]
                velocity[atm_num,1] = next_pts[atm_num,1] - atom[1]
                velocity[atm_num,2] = next_pts[atm_num,2] - atom[2]
           
            velocity_list.append(velocity)
            vel_n.append((np.sqrt(velocity[n,0]**2 + velocity[n,1]**2 + velocity[n,2]**2)))
            vel_vec_n.append(velocity[n,])
            
        return velocity_list,vel_n,vel_vec_n
    
    def func_grp_ke(self,func_list):
        
        vel = self.velocity(0)
        mass = self.generate_mass_list()
        partial_vel = np.zeros((len(func_list),3))
        partial_atomlist = np.zeros((len(func_list),1))
        func_grp_ke_list = []
        
        for each_pts in vel[0]:
        
            for i,atom_number in enumerate(func_list):
                
                atom_number -= 1
                partial_atomlist[i,0] = mass[atom_number]
                partial_vel[i,0] = each_pts[atom_number,0]
                partial_vel[i,1] = each_pts[atom_number,1]
                partial_vel[i,2] = each_pts[atom_number,2]
                
            func_grp_ke_list.append(0.5*np.sum(partial_vel**2*partial_atomlist))
                
        return func_grp_ke_list
    
    def theta(self,n):
        
        n -= 1
        theta_list = []
        displacement_vector = self.coord[-1][n] - self.coord[0][n]
        displacement_vector_norm = np.linalg.norm(displacement_vector)
        
        for each_momentum_vector in self.velocity(n)[2]:
            
            momentum_norm = np.linalg.norm(each_momentum_vector)
            cos_theta = displacement_vector.dot(each_momentum_vector)/(displacement_vector_norm*momentum_norm)
            theta = np.arccos(cos_theta)
            theta_list.append(180*theta/np.pi)
            
        return theta_list
    
    def force_along_direction(self,n):
        
        force_prj_list = []
        force_list = self.velocity(n)[2]
        for i,each_theta in enumerate(self.theta(n)):
            
            force_magnitude = np.linalg.norm(force_list[i])
            force_prj_list.append(force_magnitude*np.cos(each_theta))
            
        return force_prj_list
    
class Uphill_trajectory(Trajectory):
    
    def traj_filter(self,traj_plot_rule,total_traj,traj_plot_amount):

        if traj_plot_rule == 'all':
            
            self.plotting_status = True
            
        elif traj_plot_rule == 'defined':
            
            if len(self.structname) != 1:
                
                self.plotting_status = True
                
        elif type(traj_plot_rule) == str:
            
            if traj_plot_rule in self.structname:
                
                self.plotting_status = True
                
class Downhill_trajectory(Trajectory):
        
    def traj_filter(self,traj_plot_rule,total_traj,traj_plot_amount):

        if traj_plot_rule == 'all':
            
            self.plotting_status = True
            
        elif traj_plot_rule == 'defined':
            
            if len(self.end_point_name) != 1:
                
                self.plotting_status = True
                
        elif traj_plot_rule == 'not_recrossing' :
            
            if len(self.end_point_name) ==  2 and 'Unknown' not in self.end_point_name:
                
                self.plotting_status = True
                
        elif type(traj_plot_rule) == str:
            
            if traj_plot_rule in self.structname:
                
                self.plotting_status = True
                
    def get_endpoint_type(self):
        
        for i in list(self.structname):
            
            self.end_point_type.add(i.type)
        
            print(self.end_point_type)
    
    def recrossing_check(self):
        
        pass
            
                
            
                
        
                
            
