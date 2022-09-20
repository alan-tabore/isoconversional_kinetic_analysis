# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 17:26:08 2022

@author: alan.tabore

Based on N. Sbirrazzuoli, « Is the Friedman Method Applicable to Transformations with Temperature Dependent Reaction Heat? », Macromolecular Chemistry and Physics, vol. 208, nᵒ 14, p. 1592‑1597, 2007, doi: 10.1002/macp.200700100.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from bisect import bisect_left

def find_closest_value_of_conversion(myNumber, myList):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0, myList[0]
    if pos == len(myList):
        return len(myList), myList[-1]
    before_index=pos - 1
    before = myList[pos - 1]
    after_index=pos
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after_index, after
    else:
        return before_index, before
    
    
def friedman_method(min_conv, max_conv, number_of_points, conversions, rates, temperatures):
    """
    

    Parameters
    ----------
    min_conv : Float
        Minimum conversion at which the activation energy will be computed.
    max_conv : Float
        Maximum conversion at which the activation energy will be computed.
    number_of_points : Float
        Number of conversion points at which the activation energy will be computed.
    conversions : List
        List containing a list of evolution of conversion for various heating rates. Size (n,m) with n the number of heating rates and m the number of conversion points.
    rates : List
        List containing a list of evolution of rates for various heating rates. Size (n,m) with n the number of heating rates and m the number of time points.
    temperatures : List
        List containing a list of evolution of temperature for various heating rates. Size (n,m) with n the number of heating rates and m the number of temperature points.


    Returns
    -------
    conv_list : List
        List containing conversion points at which activation energy was assessed.
    Ea_calculated : List
        List containing assessed activation energy.
    Intercept : List
        List containing the intercept from the linear regressions 
        of ln(dx/dt) as function of 1/T.

    """
    conv_list = np.linspace(min_conv, max_conv, number_of_points)  #Creates the list of conversion points at which the activation energy will be computed
    Ea_calculated = []
    intercept=[]
    
    
    
    for i in range(len(conv_list)):
        one_over_T=[]
        rate=[]
        for k in range(len(conversions)):
            index,value=find_closest_value_of_conversion(conv_list[i], conversions[k])
            one_over_T.append(1/temperatures[k][index])
            rate.append(rates[k][index])
        result = stats.linregress(one_over_T, np.log(rate))
        Ea_calculated.append(result.slope*(-8.314))
        intercept.append(result.intercept)

    return conv_list, Ea_calculated, intercept








        
if __name__ == "__main__":

    import two_parallel_first_order_reaction_modelling as tpform
    import interpolation
    

    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
                  'tab:purple', 'tab:brown', 'tab:pink', 'tab:grey', 'tab:olive', 'tab:cyan']
    linestyle = ['dashed', ]
    marker = ['o']


    number_of_points = 10000
    time_list = np.linspace(0, 1800, number_of_points)  # 30 minutes
    temperatures = [np.linspace(293, 443, number_of_points),  # 5°C/min
                   np.linspace(293, 593, number_of_points),  # 10°C/min
                   np.linspace(293, 743, number_of_points),  # 15°C/min
                   np.linspace(293, 1093, number_of_points)]  # 20°C/min
    conversions = []
    rates=[]
    Ea_list = []
    
    fig, ax1 = plt.subplots(num=1)
    ax1_bis = ax1.twinx()

    
    fig, ax2 = plt.subplots(num=2)
    ax2_bis = ax2.twinx()


    fig, ax3 = plt.subplots(num=3)
    ax3_bis = ax3.twinx()


    for i in range(len(temperatures)):
        
        alpha, alpha1, alpha2, rate, rate1, rate2 = tpform.compute_alpha(time_list,
                                                                         temperatures[i],
                                                                         1.666e8,
                                                                         80000,
                                                                         1.666e13,
                                                                         120000)
        conversions.append(alpha)
        rates.append(rate)
        Ea_list.append(tpform.global_energy_of_activation(temperatures[i], alpha1, alpha2, 1.666e8, 80000, 1.666e13, 120000))
        
        plt.figure(1)
        ax1.set_xlabel('time')
        ax1.plot(time_list, conversions[i],label="conversion for "+str(5*i+5)+"°C/min")
        ax1.set_ylabel('conversion')
        ax1_bis.plot(time_list,rates[i],label="rate for "+str(5*i+5)+"°C/min", linestyle='dashed')
        ax1_bis.set_ylabel('rate')
        ax1.legend()
        ax1_bis.legend()
        
        plt.figure(2)
        ax2.set_xlabel('conversion')
        ax2.plot(conversions[i], Ea_list[i],label="Ea for "+str(5*i+5)+"°C/min")
        ax2.set_ylabel('Ea')
        ax2.legend()
        
        plt.figure(3)
        ax3.set_xlabel('temperature')
        ax3.plot(temperatures[i], Ea_list[i],label="Ea for "+str(5*i+5)+"°C/min",linestyle='dashed')
        ax3.set_ylabel('Ea')
        ax3_bis.plot(temperatures[i], conversions[i], label="conversion for "+str(5*i+5)+"°C/min")
        ax3_bis.set_ylabel('conversion')
        ax3.legend()
        ax3_bis.legend()

        
    
    

    
    new_conversions, new_rates, new_temperatures = interpolation.cubic_spline(conversions, rates, temperatures, 20000)
    
# =============================================================================
#     Plot 4
# =============================================================================

    fig, ax4 = plt.subplots(num=4)
    ax4_bis = ax4.twinx()
    
    for i in range(len(new_conversions)):
        plt.figure(4)
        ax4.set_xlabel('conversion')
        ax4.plot(new_conversions[i], new_rates[i], label="rate for "+str(5*i+5)+"°C/min")
        ax4.set_ylabel('rate')
        ax4_bis.plot(new_conversions[i], new_temperatures[i], label="temperature for "+str(5*i+5)+"°C/min", linestyle='dashed')
        ax4_bis.set_ylabel('temperature')
        ax4.legend()
        ax4_bis.legend()
    
    my_convs, my_Ea, my_intercept = friedman_method(0.001, 0.99, 10000, new_conversions, new_rates, new_temperatures)
    
# =============================================================================
#     Plot 5
# =============================================================================
    
    ax2.plot(my_convs, my_Ea, label="Ea from friedman method")
    ax2.legend()
    
    fig5, ax5 = plt.subplots(num=5)
    ax5.plot(my_convs, np.exp(my_intercept), label="Intercept from friedman method")
    ax5.legend()
    
    A1=np.exp(my_intercept)

    
    
                
            
        