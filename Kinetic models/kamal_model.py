
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 10:09:53 2022

@author: alan.tabore
"""
import numpy as np
import matplotlib.pyplot as plt

def rate_constant(A,Ea,T):
    """
    Compute the value of rate constant with the Arrhenius equation

    Parameters
    ----------
    A : Float
        Pre-exponential factor.
    Ea : Float
        Activation energy.
    T : Float
        Temperature of reaction.

    Returns
    -------
    k : Float
        rate constant for the Arrhenius equation.

    """
    old_settings=np.seterr(all="ignore")
    k=A*np.exp(-Ea/(8.314*T))
    np.seterr(**old_settings)
    return k

def rate_for_kamal(A1,E1,A2,E2,m,n,alpha,T):
    """
    Compute the value of the rate of reaction for a Kamal equation.

    Parameters
    ----------
    A1 : Float
        Pre-exponential factor of the regular reaction.
    E1 : Float
        Activation energy of the regular reaction.
    A2 : Float
        Pre-exponential factor of the autocatalyzed reaction.
    E2 : Float
        Activation energy of the autocatalyzed reaction.
    m : TYPE
        Order of reaction for the autocatalyzed reaction.
    n : TYPE
        Order of reaction for the regular reaction.
    alpha : Float
        Extent of reaction.
    T : Float
        Temperature of reaction.

    Returns
    -------
    rate : Float
        Rate of reaction for a Kamal equation.

    """
    old_settings=np.seterr(all="ignore")
    rate = ( rate_constant(A1,E1,T) + (rate_constant(A2,E2,T)*alpha**m))*(1-alpha)**n
    np.seterr(**old_settings)
    return rate

def compute_alpha(time,temperature,A1,E1,A2,E2,m,n):
    """
    Compute the extent and rate of reaction for the Kamal equation.

    Parameters
    ----------
    time : List
        List containing times during reaction.
    temperature : List
        List containing temperatures during reaction.
     A1 : Float
         Pre-exponential factor of the regular reaction.
     E1 : Float
         Activation energy of the regular reaction.
     A2 : Float
         Pre-exponential factor of the autocatalyzed reaction.
     E2 : Float
         Activation energy of the autocatalyzed reaction.
     m : TYPE
         Order of reaction for the autocatalyzed reaction.
     n : TYPE
         Order of reaction for the regular reaction.

    Returns
    -------
    extent : List
        List containing the evolution of extent during reaction.
    rate : List
        List containing the evolution of rate during reaction.

    """
    
    extent=[0] #at the beginning of reaction extent is equal to 0
    rate=[rate_for_kamal(A1,E1,A2,E2,m,n,0,temperature[0])] #rate at the beginning of reaction
    
    for i in range(len(time)-1):
        extent_for_next_step = extent[i] + rate_for_kamal(A1,E1,A2,E2,m,n,extent[i],temperature[i])*(time[i+1]-time[i])
        if extent_for_next_step > 1:
            extent.append(1)
            rate.append(0)
        else: 
            extent.append( extent_for_next_step )
            rate.append(rate_for_kamal(A1,E1,A2,E2,m,n,extent[i+1],temperature[i+1]))     
        
    return extent, rate

def compute_Ea(A1,E1,A2,E2,m,n,alpha,T):
    return (E1*rate_constant(A1, E1, T)+E2*rate_constant(A2, E2, T)*alpha**m)/(rate_constant(A1, E1, T)+rate_constant(A2, E2, T)*alpha**m)



#%%
if __name__ == "__main__": 

    
    number_of_points = 10000
    time = np.linspace(0, 1800, number_of_points)  # 30 minutes
    times = [time, time, time, time]
    temperatures = [np.linspace(293, 443, number_of_points),  # 5°C/min
                   np.linspace(293, 593, number_of_points),  # 10°C/min
                   np.linspace(293, 743, number_of_points),  # 15°C/min
                   np.linspace(293, 1093, number_of_points)]  # 20°C/min
    conversions=[]
    rates=[]

    fig, ax1 = plt.subplots(num=1)
    ax1_bis = ax1.twinx()
    ax1.set_xlabel("time")
    ax1.set_ylabel("conversion")
    ax1_bis.set_ylabel("rate")

    for i in range(len(temperatures)):
        extent,rate = compute_alpha(time,temperatures[i],1.666e8,80000,1.666e13,120000,1,0.7)
        conversions.append(extent)
        rates.append(rate)
        ax1.plot(time,extent,label="conversion for a heating rate of " + str((i+1)*5) + "°/min",linestyle='dotted')
        ax1_bis.plot(time,rate,label="rate for a heating rate of" + str((i+1)*5) + "°/min")
        ax1.legend()
        ax1_bis.legend()
        
#%%       
    my_time=np.linspace(0, 1800, number_of_points)  # 30 minutes
    my_temp=np.linspace(293, 593, number_of_points)  # 10°C/min
    
    fig, ax2 = plt.subplots(num=2)
    ax2_bis=ax2.twinx()
    ax2.set_xlabel("time (s)")
    ax2.set_ylabel("conversion")
    ax2_bis.set_ylabel("rate")
    # line,=ax2.plot([], [])
    
    points_for_m=10
    points_for_n=10
    points_for_E1=10
    points_for_E2=10
    
    for m in range(points_for_m):
        for n in range(points_for_n):
            for E1 in range(points_for_E1):
                for E2 in range(points_for_E2):
                    ax2.clear()
                    ax2_bis.clear()
                    ax2.set_xlabel("time (s)")
                    ax2.set_ylabel("conversion")
                    ax2_bis.set_ylabel("rate")
                    extent,rate = compute_alpha(my_time, my_temp, 1.666e8,\
                                                80000+E1*(120000-80000)/points_for_E1, 1.666e13,\
                                                    100000+E2*(300000-100000)/points_for_E2,\
                                                        0.1+m*(3-0.1)/points_for_m,\
                                                            0.1+n*(3-0.1)/points_for_n)

                    ax2.plot(my_time,extent,label="E1="+str((80000+E1*(120000-80000)/points_for_E1))+"\n"+ \
                             "E2="+str((100000+E2*(300000-100000)/points_for_E2))+"\n"+ \
                                 "n="+str(0.1+n*(3-0.1)/points_for_n)+"\n"+ \
                                     "m="+str(0.1+m*(3-0.1)/points_for_m))
                    ax2_bis.plot(my_time,rate,color='tab:orange')
                    ax2.legend()

                    plt.draw()
                    plt.pause(.001)
                                 
                        
                        
                    

        