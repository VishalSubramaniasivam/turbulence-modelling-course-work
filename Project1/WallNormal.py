# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 20:06:05 2020

@author: VISHAL SUBRAMANIASIVAM
"""
def WallNormal(xf2d,yf2d):
    import numpy as np
# Bottom Wall
    sx_bottom = np.diff(xf2d[:,0],axis=0)
    sy_bottom = np.diff(yf2d[:,0],axis=0)
    
    d=np.sqrt(sx_bottom**2+sy_bottom**2)
    nx_bottom=-sy_bottom/d
    ny_bottom=sx_bottom/d
# Top Wall
    sx_top = np.diff(xf2d[:,-1],axis=0)
    sy_top = np.diff(yf2d[:,-1],axis=0)
    
    d=np.sqrt(sx_top**2+sy_top**2)
    nx_top=-sy_top/d
    ny_top=sx_top/d    
    
    nx_bottom = np.insert(nx_bottom,0,nx_bottom[0],axis=0)
    ny_bottom = np.insert(ny_bottom,0,ny_bottom[0],axis=0)
    nx_top = np.insert(nx_top,0,nx_top[0],axis=0)
    ny_top = np.insert(ny_top,0,ny_top[0],axis=0)
    nx_bottom = np.insert(nx_bottom,-1,nx_bottom[-1],axis=0)
    ny_bottom = np.insert(ny_bottom,-1,ny_bottom[-1],axis=0)
    nx_top = np.insert(nx_top,-1,nx_top[-1],axis=0)
    ny_top = np.insert(ny_top,-1,ny_top[-1],axis=0)    
    
    return nx_bottom,ny_bottom,nx_top,ny_top