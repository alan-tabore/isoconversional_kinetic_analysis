# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:12:12 2022

@author: alan.tabore

Based on S. Vyazovkin, « Evaluation of activation energy of thermally stimulated solid-state reactions under arbitrary variation of temperature », Journal of Computational Chemistry, vol. 18, nᵒ 3, p. 393‑402, 1997, doi: 10.1002/(SICI)1096-987X(199702)18:3<393::AID-JCC9>3.0.CO;2-P.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
import time
from bisect import bisect_left
from progress_bar import printProgressBar
import interpolation










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




def get_data_at_conversion(conversion,list_of_conversions,list_of_times,list_of_temperatures):
    new_list_of_conversions=[]
    new_list_of_times=[]
    new_list_of_temperatures=[]
    for i in range(len(list_of_conversions)):
        index,value = find_closest_value_of_conversion(conversion, list_of_conversions[i])
        new_list_of_conversions.append(list_of_conversions[i][:index+1])
        new_list_of_times.append(list_of_times[i][:index+1])
        new_list_of_temperatures.append(list_of_temperatures[i][:index+1])

    
    return new_list_of_conversions, new_list_of_times, new_list_of_temperatures
        
    

def integral_J_to_compute(Ea, time_list, temperature_list):
    """
    Compute the integral described in Vyazovkin paper

    Parameters
    ----------
    Ea : float
        Activation energy that will be optimized.
    time_list : list
        list of times for the DSC scan on which integrals are calculated.
    temperature_list : list
        list of times for the DSC scan on which integrals are calculated.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    int_sum=0
    for i in range(len(time_list)-1):
        int_sum = int_sum + np.exp( -Ea / (8.314* (temperature_list[i]+temperature_list[i+1]) /2) )*(time_list[i+1]-time_list[i])

    return int_sum



def func_to_minimize(Ea,list_of_time_list, list_of_temperature_list):
    """
    Compute the function to minimize described in Vyazovkin paper

    Parameters
    ----------
    Ea : float
        Energy of activation for which the function will be computed
    list_of_time_list : list
        List of time lists for multiple experiments carried out with different temperature program.
    list_of_temperature_list : list
        List of temperature lists for multiple experiments carried out with different temperature program.

    Returns
    -------
    func_sum : float
        Value of the function for the given energy of activation.

    """
    my_integrals_J=[]
    func_sum=0

    for i in range(len(list_of_time_list)):
        my_integrals_J.append(integral_J_to_compute(Ea, list_of_time_list[i], list_of_temperature_list[i]))

    for i in range(len(my_integrals_J)):
        for k in range(len(my_integrals_J)):
            if i!=k:
                # print("i=",i," k=",k)
                # print(integral_J_to_compute(Ea, list_of_time_list[i], list_of_temperature_list[i]))
                # print(integral_J_to_compute(Ea, list_of_time_list[k], list_of_temperature_list[k]))
                # print(integral_J_to_compute(Ea, list_of_time_list[i], list_of_temperature_list[i])/integral_J_to_compute(Ea, list_of_time_list[k], list_of_temperature_list[k])) 
                func_sum = func_sum + (my_integrals_J[i]/my_integrals_J[k])
                # print("func sum=",func_sum)   
    return func_sum


def get_min_func_as_function_of_conversion_and_Ea(min_conv,max_conv,number_of_points,conv_lists,time_lists,temperature_lists):
    conv_list=np.linspace(min_conv, max_conv,number_of_points)
    Ea_to_test=np.linspace(10000,210000,2000)
    results_list=[]
    
    for i in range(len(conv_list)):
        results_of_function_to_minimize=[]
        conv, time, temperature = get_data_at_conversion(conv_list[i], conv_lists, time_lists, temperature_lists)
        for k in range(len(Ea_to_test)):
            results_of_function_to_minimize.append(func_to_minimize(Ea_to_test[k], time, temperature))
        
        results_list.append(results_of_function_to_minimize)
        
    return conv_list, Ea_to_test,results_list

def get_Ea_as_function_of_conversion(min_conv,max_conv,number_of_points,conv_lists,time_lists,temperature_lists):
    conv_list=np.linspace(min_conv, max_conv,number_of_points)
    Ea_to_test=np.linspace(10000,210000,2000)
    
    Ea_calculated=[]
    
    for i in range(len(conv_list)):
        results_of_function_to_minimize=[]
        conv, time, temperature = get_data_at_conversion(conv_list[i], conv_lists, time_lists, temperature_lists)
        for k in range(len(Ea_to_test)):
            results_of_function_to_minimize.append(func_to_minimize(Ea_to_test[k], time, temperature))
        
        minimum=min(results_of_function_to_minimize)
        index_of_mininum=results_of_function_to_minimize.index(minimum)
        Ea_calculated.append( Ea_to_test[index_of_mininum])

        
    return conv_list, Ea_calculated


def optimize_BFGS(initial_guess,min_conv,max_conv,number_of_points,conv_lists,time_lists,temperature_lists):
    t_start = time.process_time()
    print("Progress of optimization with BFGS method:")
    conv_list=np.linspace(min_conv, max_conv,number_of_points)
    Ea_calculated=[]
    for i in range(len(conv_list)):
        conv, time_, temperature = get_data_at_conversion(conv_list[i], conv_lists, time_lists, temperature_lists)
        result=scipy.optimize.minimize(func_to_minimize, x0=initial_guess, args=(time_, temperature), method='BFGS', options={'gtol': 1e-15} )
        Ea_calculated.append(float(result.x))
        printProgressBar(i, len(conv_list), length=50)    
    t_stop = time.process_time()
    print("Done ! It took ", t_stop-t_start,"s for the BFGS optimization to be performed")
    return conv_list, Ea_calculated

def optimize_CG(initial_guess,min_conv,max_conv,number_of_points,conv_lists,time_lists,temperature_lists):
    t_start = time.process_time()
    print("Progress of optimization with CG method:")
    conv_list=np.linspace(min_conv, max_conv,number_of_points)
    Ea_calculated=[]
    for i in range(len(conv_list)):
        conv, time_, temperature = get_data_at_conversion(conv_list[i], conv_lists, time_lists, temperature_lists)
        result=scipy.optimize.minimize(func_to_minimize, x0=initial_guess, args=(time_, temperature), method='CG', options={'gtol': 1e-15} )
        Ea_calculated.append(float(result.x))
        initial_guess=float(result.x)
        # print("\r", "Algorithm completed at "+str((i+1)/len(conv_list))+"%", end="")
        printProgressBar(i, len(conv_list),length=50)
    t_stop = time.process_time()
    print("Done ! It took ", t_stop-t_start,"s for the CG optimization to be performed")
    return conv_list, Ea_calculated

        
if __name__ == "__main__":    
    import two_parallel_first_order_reaction_modelling as tpform
    

    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
                  'tab:purple', 'tab:brown', 'tab:pink', 'tab:grey', 'tab:olive', 'tab:cyan']
    linestyle = ['dashed', ]
    marker = ['o']

    issou1 = []
    issou2 = []

    number_of_points = 10000
    time_list = np.linspace(0, 1800, number_of_points)  # 30 minutes
    time_lists = [time_list, time_list, time_list, time_list]
    temperature = [np.linspace(293, 443, number_of_points),  # 5°C/min
                   np.linspace(293, 593, number_of_points),  # 10°C/min
                   np.linspace(293, 743, number_of_points),  # 15°C/min
                   np.linspace(293, 893, number_of_points)]  # 20°C/min
    conversions = []
    Ea_list = []

    ax = plt.axes()

    for i in range(len(temperature)):
        alpha, alpha1, alpha2, rate, rate1, rate2 = tpform.compute_alpha(
            time_list, temperature[i], 1.666e8, 80000, 1.666e13, 120000)
        conversions.append(alpha)
        Ea_list.append(tpform.global_energy_of_activation(
            temperature[i], alpha1, alpha2, 1.666e8, 80000, 1.666e13, 120000))
        plt.figure(1)
        plt.plot(time_list, conversions[i],
                 label="conversion for "+str(5*i+5)+"°C/min")
        plt.figure(2)
        plt.plot(conversions[i], Ea_list[i],
                 label="Ea for "+str(5*i+5)+"°C/min")
    

    
    new_conversions, new_times, new_temperatures = interpolation.linear_interpolation(conversions, time_lists, temperature, 20000)
    new_conversions2, new_times2, new_temperatures2 = interpolation.cubic_spline(conversions, time_lists, temperature, 20000)
    
    
    plt.figure(1)
    plt.legend()

    # t1_start = time.process_time()
    # conv_to_plot, Ea_to_plot = get_Ea_as_function_of_conversion(
    #     0.001, 0.99, 2500, new_conversions, new_times, new_temperatures)
    # t1_stop = time.process_time()


    # t2_start =time.process_time()
    # conv_to_plot2, Ea_to_plot2 = get_Ea_as_function_of_conversion(
    #     0.001, 0.99, 2500, new_conversions2, new_times2, new_temperatures2)
    # t2_stop =time.process_time()


    t3_start = time.process_time()
    conv_to_plot3, Ea_to_plot3 = optimize_CG(
        50000, 0.001, 0.99, 4000, new_conversions2, new_times2, new_temperatures2)
    t3_stop = time.process_time()
    
    # t4_start = time.process_time()
    # conv_to_plot4, Ea_to_plot4 = optimize_BFGS(
    #     50000, 0.001, 0.99, 2500, new_conversions2, new_times2, new_temperatures2)
    # t4_stop = time.process_time()

    plt.figure(2)

    # plt.plot(conv_to_plot, Ea_to_plot, label="Vyazovkin modified integral (linear interpolation)")
    # plt.plot(conv_to_plot2,Ea_to_plot2,label="Vyazovkin modified integral (cubic spline)")
    plt.plot(conv_to_plot3, Ea_to_plot3, label="Vyazovkin modified integral (cubic spline) + CG")
    # plt.plot(conv_to_plot4, Ea_to_plot4, label="Vyazovkin modified integral (cubic spline) + BFGS")
    plt.legend()

    # print("Elapsed time during the linear interpolation:", t1_stop-t1_start)
    # print("Elapsed time during the cubic spline in seconds:", t2_stop-t2_start)
    print("Elapsed time during the cubic spline + CG method in seconds:", t3_stop-t3_start)
    # print("Elapsed time during the cubic spline + BFGS method in seconds:", t4_stop-t4_start)
    
   
    
    
    
    

    

    
    # for i in range(len(fixed_convs)):
    #     print(fixed_convs[i])
    #     print(minimize(func_to_minimize, 60000, args=(my_times,my_convs) ,method='Newton-CG' ))
        
    
    
    # fig=plt.figure()
    # axes=plt.axes()

    # for i in range(len(my_times)):
    #     plt.plot(my_times[i],my_convs[i],color=color_list[i],label="v="+str((i+1)*5)+"°C/min considering complete reaction")
    #     plt.plot(my_times[i],my_convs2[i],color=color_list[i],linestyle=linestyle[0],label="v="+str((i+1)*5)+"°C/min considering partial reaction")
    # plt.legend()
    # plt.show()
    
    