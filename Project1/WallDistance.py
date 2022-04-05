# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 18:22:49 2020

@author: VISHAL SUBRAMANIASIVAM
"""
def WallDistance(x2d,y2d,point_x,point_y):
    import numpy as np
    ni = len(x2d)
    dist_LW = np.zeros(ni)
    dist_UP = np.zeros(ni)
    for i in range(ni):
        dist_LW[i] = np.sqrt(((point_x-x2d[i,0])**2) + ((point_y-y2d[i,0])**2))
        dist_UP[i] = np.sqrt(((point_x-x2d[i,-1])**2) + ((point_y-y2d[i,-1])**2))
        
    if min(dist_LW)<=min(dist_UP):
        min_dist = min(dist_LW)
    else:
        min_dist = min(dist_UP)
        
    return min_dist
