# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 15:17:54 2022

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
    k=A*np.exp(-Ea/(8.314*T))
    return k

def rate_for_two_parallel_first_order_reactions(A1,E1,A2,E2,T,alpha1,alpha2):
    """
    Compute the global rate of reaction for two parallel reaction with the same contribution of both reaction to the global reaction rate.

    Parameters
    ----------
    A1 : Float
        Pre-exponential factor of reaction n°1.
    E1 : Float
        Activation energy of reaction n°1.
    A2 : Float
        Pre-exponential factor of reaction n°2.
    E2 : Float
        Activation energy of reaction n°2.
    T : Float
        Temperature of reaction.
    alpha1 : Float
        Extent of reaction n°1.
    alpha2 : Float
        Extent of reaction n°2.

    Returns
    -------
    rate : Float
        Global rate of reaction for two parallel reaction with the same contribution of both reaction to the global reaction rate.

    """
    rate = (1/2)* (rate_constant(A1,E1,T)*(1-alpha1) + rate_constant(A2,E2,T)*(1-alpha2))
    return rate
    

    
def compute_alpha(time,temperature,A1,E1,A2,E2):
    """
    Compute the extent and rate of reaction for two parallel reaction with the same contribution of both reaction to the global reaction rate.

    Parameters
    ----------
    time : List
        List containing times during reaction.
    temperature : List
        List containing temperatures during reaction.
    A1 : Float
        Pre-exponential factor of reaction n°1.
    E1 : Float
        Activation energy of reaction n°1.
    A2 : Float
        Pre-exponential factor of reaction n°2.
    E2 : Float
        Activation energy of reaction n°2.

    Returns
    -------
    extent : List
        List containing the global extent of reaction.
    extent1 : List
        List containing the extent of reaction n°1.
    extent2 : List
        List containing the extent of reaction n°2.
    rate : List
        List containing the global rate of reaction.
    rate1 : List
        List containing the rate of reaction n°1.
    rate2 : List
        List containing the rate of reaction n°2.

    """
    extent=[0]
    extent1=[0]
    extent2=[0]
    rate=[rate_for_two_parallel_first_order_reactions(A1,E1,A2,E2,temperature[0],0,0)]
    rate1=[rate_constant(A1,E1,temperature[0])]
    rate2=[rate_constant(A2,E2,temperature[0])]
    
    
    for i in range(len(time)-1):
        next_extent=extent[i] + rate_for_two_parallel_first_order_reactions(A1,E1,A2,E2,temperature[i],extent1[i],extent2[i])*(time[i+1]-time[i])
        
        if next_extent > 1:
            extent.append(1)
            extent1.append(1)
            extent2.append(1)
            rate.append(0)
            rate1.append(0)
            rate2.append(0)
        
        else :    
            extent.append( extent[i] + rate_for_two_parallel_first_order_reactions(A1,E1,A2,E2,temperature[i],extent1[i],extent2[i])*(time[i+1]-time[i]) )
            extent1.append( extent1[i] +  rate_constant(A1,E1,temperature[i])*(1-extent1[i])*(time[i+1]-time[i]) )    
            extent2.append( extent2[i] + rate_constant(A2,E2,temperature[i])*(1-extent2[i])*(time[i+1]-time[i]) )
            rate.append(rate_for_two_parallel_first_order_reactions(A1,E1,A2,E2,temperature[i],extent1[i],extent2[i]))
            rate1.append(rate_constant(A1,E1,temperature[i])*(1-extent1[i]))
            rate2.append(rate_constant(A2,E2,temperature[i])*(1-extent2[i]))
            
            
        
    return extent, extent1, extent2, rate, rate1, rate2



        
def global_energy_of_activation(temperature,alpha1,alpha2,A1,E1,A2,E2):
    """
    Return the evolution of apparent activation energy for two parallel reactions.
    It's only valid for a given couple of extent list and temperature list.
    Ea=f(extent,temperature)

    Parameters
    ----------
    temperature : List
        List of temperatures during reaction.
    alpha1 : List
        List of extent of reaction n°1 during reaction.
    alpha2 : List
        List of extent of reaction n°2 during reaction.
    A1 : Float
        Pre-exponential factor of reaction n°1.
    E1 : Float
        Activation energy of reaction n°1.
    A2 : Float
        Pre-exponential factor of reaction n°2.
    E2 : Float
        Activation energy of reaction n°2.

    Returns
    -------
    Ea_to_return : List
        List of the evolution of apparent activation energy for two parallel reactions.

    """
    Ea_to_return=[]

    for i in range(len(alpha1)):
        if alpha1[i] == 1.0 and alpha2[i]==1.0:
            Ea_to_return.append(0)
        else:
            
            Ea_to_return.append( (E1*rate_constant(A1, E1, temperature[i])*(1-alpha1[i]) + E2*rate_constant(A2, E2, temperature[i])*(1-alpha2[i]) )/
                      (rate_constant(A1, E1, temperature[i])*(1-alpha1[i]) + rate_constant(A2, E2, temperature[i])*(1-alpha2[i]) ))
    return Ea_to_return
    
    
    
    
if __name__ == "__main__":     
    
    time_list=np.linspace(0, 5000, 20000)
    
    temperature1=np.linspace(353,393,10000)
    temperature2=np.full(10000, 393)
    temperature = np.concatenate((temperature1 , temperature2))
    
    # temperature = np.linspace(353,413,2000)
    
    # temperature = np.full(2000,363)
    
    alpha,alpha1,alpha2,rate,rate1,rate2 = compute_alpha(time_list, temperature, 1.666e8, 80000, 1.666e13, 120000)
    Ea_list=global_energy_of_activation(temperature, alpha1, alpha2, 1.666e8, 80000, 1.666e13, 120000)
    
    
# =============================================================================
#  Plot 1   
# =============================================================================


    fig=plt.figure()
    ax=plt.axes()
    
    ax.plot(time_list,alpha,label="extent of global reaction")
    ax.plot(time_list,alpha1,label="extent of reaction 1")
    ax.plot(time_list,alpha2,label="extent of reaction 2")
    
    ax.set_xlabel('time (s)')
    ax.set_ylabel('extent of reaction')
    plt.legend()
    
    ax2=plt.twinx()
    ax2.plot(time_list,temperature,color='r', label="temperature")
    ax2.set_ylabel('Temperature (K)')
    plt.legend()
    
    plt.show()
# =============================================================================
#    Plot 2   
# =============================================================================
    fig=plt.figure()
    ax=plt.axes()
    
    ax.plot(alpha,Ea_list,label="Activation energy")
    ax.set_xlabel('Conversion')
    ax.set_ylabel('Ea (J/mol)')
    plt.legend()
    ax2=plt.twinx()
    ax2.plot(alpha,temperature,color='r', label="temperature")
    ax2.set_ylabel('Temperature (K)')
    plt.legend()
    
# =============================================================================
#     Plot 3
# =============================================================================
    
    fig=plt.figure()
    ax=plt.axes()
    ax.set_xlabel('time (s)')
    ax.set_ylabel('rate')
    ax.plot(time_list,rate,label="rate of global reaction")
    ax.plot(time_list,rate1,label="rate of reaction 1")
    ax.plot(time_list,rate2,label="rate of reaction 2")
    plt.legend()
    
    ax2=plt.twinx()
    ax2.plot(time_list,temperature,color='r', label="temperature")
    ax2.set_ylabel('Temperature (K)')
    plt.legend()

    
    