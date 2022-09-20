# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:25:15 2022

@author: alan.tabore

Based on S. Vyazovkin, « Modification of the integral isoconversional method to account for variation in the activation energy », Journal of Computational Chemistry, vol. 22, nᵒ 2, p. 178‑183, 2001, doi: 10.1002/1096-987X(20010130)22:2<178::AID-JCC5>3.0.CO;2-#.
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy
import time
from progress_bar import printProgressBar
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
    
def check_boundaries(minimum,maximum,list_of_real_convs):
    for i in range(len(list_of_real_convs)): #for each heating rate
        if minimum < list_of_real_convs[i][0]:
            raise ValueError("the minimum conversion you want to analyze is "+str(minimum)+" but you're input data have a minimum of "+str(list_of_real_convs[i][0])+"\Please increase the value of the minimum conversion to analyze")
        if maximum > list_of_real_convs[i][-1]:
            raise ValueError("the maximumm conversion you want to analyze is "+str(maximum)+" but you're input data have a maximum of "+str(list_of_real_convs[i][-1])+"\Please lower the value of the maximum conversion to analyze")


def get_index_area_of_calculation(conversion_lists, list_of_points_for_analysis):
    """
    Return the indexes of the interval of calculation for the modified integral J (see article)

    Parameters
    ----------
    conversion_lists : List
        List with conversions at different heating rates (one list element for one heating rate)
    list_of_points_for_analysis : List
        List containing all the point of conversion at which the activation energy will be calculated

    Raises
    ------
    ValueError
        return an error if starting point is smaller than step

    Returns
    -------
    index_list : List
        Return the indexes of the interval of calculation for the modified integral J (see article)

    """
    
    index_list = []
    start_value = 2 * \
        list_of_points_for_analysis[0] - list_of_points_for_analysis[1]

    if start_value < 0:
        raise ValueError(
            "Make sure to select a starting conversion higher than the step between two conversion points\
                For example if the step is 0.1, the starting conversion should be at least 0.1\n\
                Actual value of starting conversion: "+str(list_of_points_for_analysis[0])+"\n\
                Actual value of step: "+str(list_of_points_for_analysis[1] - list_of_points_for_analysis[0]))
                

    for i in range(len(conversion_lists)):
        index_of_start_value, value_of_start_value = find_closest_value_of_conversion(
            start_value, conversion_lists[i])
        index_list.append([index_of_start_value])
        for k in range(len(list_of_points_for_analysis)):
            index, value = find_closest_value_of_conversion(
                list_of_points_for_analysis[k], conversion_lists[i])
            index_list[i].append(index)

    return index_list


def get_data_in_interval(conversion_lists, time_lists, temperature_lists, list_of_points_for_analysis):
    """
    Returns the conversion, time and temperature lists for each heating and for each interval of calculation
    1st level of the list corresponds to heating rate
    2nd level of the list corresponds to interval of calculation

    Parameters
    ----------
    conversion_lists : List
        List with conversions at different heating rates (one list element for one heating rate).
    time_lists : List
        List with times at different heating rates (one list element for one heating rate)
    temperature_lists : List
        List with temperatures at different heating rates (one list element for one heating rate)
    list_of_points_for_analysis : List
        List containing all the point of conversion at which the activation energy will be calculated

    Returns
    -------
    new_list_of_conversions : List
        List with conversions at different heating rates for different intervals
    new_list_of_times : List
        List with times at different heating rates for different intervals
    new_list_of_temperatures : List
        List with temperatures at different heating rates for different intervals

    """
    

    new_list_of_conversions = []
    new_list_of_times = []
    new_list_of_temperatures = []
    index_list = get_index_area_of_calculation(
        conversion_lists, list_of_points_for_analysis)
    for i in range(len(index_list)):
        new_list_of_conversions.append([])
        new_list_of_times.append([])
        new_list_of_temperatures.append([])
        for k in range(len(index_list[i])-1):
            new_list_of_conversions[i].append(
                conversion_lists[i][index_list[i][k]:index_list[i][k+1]])
            new_list_of_times[i].append(
                time_lists[i][index_list[i][k]:index_list[i][k+1]])
            new_list_of_temperatures[i].append(
                temperature_lists[i][index_list[i][k]:index_list[i][k+1]])
        
    

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

    int_sum = 0
    for i in range(len(time_list)-1):
        int_sum = int_sum + np.exp(-Ea / (8.314 * (temperature_list[i]+temperature_list[i+1]) / 2))*(time_list[i+1]-time_list[i])
    return int_sum


def func_to_minimize(Ea, list_of_time_list, list_of_temperature_list):
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
    
    my_integrals_J = []
    func_sum = 0

    #For each heating rate the integral is calculated
    for i in range(len(list_of_time_list)):
        my_integrals_J.append(integral_J_to_compute(
            Ea, list_of_time_list[i], list_of_temperature_list[i]))
    
    #Calculation of the result of func to minimize
    for i in range(len(my_integrals_J)):
        for k in range(len(my_integrals_J)):
            if i != k:
                    func_sum = func_sum + (my_integrals_J[i]/my_integrals_J[k])

    return func_sum





def optimize_raw(min_conv, max_conv, number_of_points, conv_lists, time_lists, temperature_lists):
    """
    Find the energy of activation for multiple conversions by using the advanced vyazovkin method by using a predefined list of Ea to test to find the optimal Ea.

    Parameters
    ----------
    min_conv : Float
        Minimum conversion at which the activation energy will be computed.
    max_conv : Float
        Maximum conversion at which the activation energy will be computed.
    number_of_points : Float
        Number of conversion points at which the activation energy will be computed.
    conv_lists : List
        List containing a list of evolution of conversion for various heating rates. Size (n,m) with n the number of heating rates and m the number of conversion points.
    time_lists : List
        List containing a list of evolution of time for various heating rates. Size (n,m) with n the number of heating rates and m the number of time points.
    temperature_lists : List
        List containing a list of evolution of temperature for various heating rates. Size (n,m) with n the number of heating rates and m the number of temperature points.

    Returns
    -------
    conv_list : List
        List containing conversion points at which activation energy was assessed.
    Ea_calculated : List
        List containing assessed activation energy.
    average_time_at_claculation_point : List
        List containing the average time of all heating programs at the conversion points at which activation energy was assessed.
    average_temperature_at_claculation_point : List
        List containing the average temperature of all heating programs at the conversion points at which activation energy was assessed.

    """
    t_start = time.process_time()
    check_boundaries(min_conv, max_conv,conv_lists)
    print("Progress of optimization with raw Ea value method:")
    conv_list = np.linspace(min_conv, max_conv, number_of_points)
    Ea_to_test = np.linspace(10000, 210000, 200)
    Ea_calculated = []
    average_time_at_claculation_point=[]
    average_temperature_at_claculation_point=[]
    conversions_for_all_heating_rates_interval_separated, times_for_all_heating_rates_interval_separated, temperatures_for_all_heating_rates_interval_separated = get_data_in_interval(conv_lists, time_lists, temperature_lists, conv_list)
    
    
    
    for i in range(len(conv_list)):
        results_of_function_to_minimize = []
        times = [heating_rate[i] for heating_rate in times_for_all_heating_rates_interval_separated]
        temperatures = [heating_rate[i] for heating_rate in temperatures_for_all_heating_rates_interval_separated]
        
        for u in range(len(Ea_to_test)):
            results_of_function_to_minimize.append(func_to_minimize(Ea_to_test[u], times, temperatures))

        # plt.plot(Ea_to_test, results_of_function_to_minimize,label="conv="+str(conv_list[i]))
        minimum = min(results_of_function_to_minimize)
        index_of_mininum = results_of_function_to_minimize.index(minimum)
        Ea_calculated.append(Ea_to_test[index_of_mininum])
        # plt.scatter(Ea_to_test[index_of_mininum],minimum,label="min for a conv of: "+str(conv_list[i]))
        # print("minimum:", minimum)
        # print("Ea found:", Ea_to_test[index_of_mininum])
            
        printProgressBar(i, len(conv_list),length=50)
    
    t_stop = time.process_time()
    
    for i in range(len(temperatures_for_all_heating_rates_interval_separated[0])):
        time_sum=0
        temp_sum=0
        for k in range(len(temperatures_for_all_heating_rates_interval_separated)):
            time_sum = time_sum + times_for_all_heating_rates_interval_separated[k][i][-1]
            temp_sum = temp_sum + temperatures_for_all_heating_rates_interval_separated[k][i][-1]
        average_time_at_claculation_point.append(time_sum/len(times_for_all_heating_rates_interval_separated))
        average_temperature_at_claculation_point.append(temp_sum/len(temperatures_for_all_heating_rates_interval_separated))
        
    print("Done ! It took ", t_stop-t_start,"s for the optimization with raw value of Ea to be performed")
    
    return conv_list, Ea_calculated, average_time_at_claculation_point, average_temperature_at_claculation_point


def optimize_BFGS(initial_guess, min_conv, max_conv, number_of_points, conv_lists, time_lists, temperature_lists):
    """
    Find the energy of activation for multiple conversions by using the advanced vyazovkin method with a Broyden-Fletcher-Goldfarb-Shanno algorithm to find the optimal Ea.

    Parameters
    ----------
    initial_guess : Float
        Initial guess of the activation energy for the first desired conversion.
    min_conv : Float
        Minimum conversion at which the activation energy will be computed.
    max_conv : Float
        Maximum conversion at which the activation energy will be computed.
    number_of_points : Float
        Number of conversion points at which the activation energy will be computed.
    conv_lists : List
        List containing a list of evolution of conversion for various heating rates. Size (n,m) with n the number of heating rates and m the number of conversion points.
    time_lists : List
        List containing a list of evolution of time for various heating rates. Size (n,m) with n the number of heating rates and m the number of time points.
    temperature_lists : List
        List containing a list of evolution of temperature for various heating rates. Size (n,m) with n the number of heating rates and m the number of temperature points.

    Returns
    -------
    conv_list : List
        List containing conversion points at which activation energy was assessed.
    Ea_calculated : List
        List containing assessed activation energy.
    average_time_at_claculation_point : List
        List containing the average time of all heating programs at the conversion points at which activation energy was assessed.
    average_temperature_at_claculation_point : List
        List containing the average temperature of all heating programs at the conversion points at which activation energy was assessed.

    """
    t_start = time.process_time()
    check_boundaries(min_conv, max_conv,conv_lists)
    print("Progress of optimization with BFGS method:")
    conv_list = np.linspace(min_conv, max_conv, number_of_points)
    Ea_calculated = []
    average_time_at_claculation_point=[]
    average_temperature_at_claculation_point=[]
    conversions_for_all_heating_rates_interval_separated, times_for_all_heating_rates_interval_separated, temperatures_for_all_heating_rates_interval_separated = get_data_in_interval(conv_lists, time_lists, temperature_lists, conv_list)
    
    for i in range(len(conv_list)):
        times = [heating_rate[i] for heating_rate in times_for_all_heating_rates_interval_separated]
        temperatures = [heating_rate[i] for heating_rate in temperatures_for_all_heating_rates_interval_separated]
        
        result = scipy.optimize.minimize(func_to_minimize, x0=initial_guess, args=(
            times, temperatures), method='BFGS', options={'gtol': 1e-15})

        Ea_calculated.append(float(result.x))
        initial_guess = float(result.x)

        printProgressBar(i, len(conv_list), length=50)
    t_stop = time.process_time()
    
    for i in range(len(temperatures_for_all_heating_rates_interval_separated[0])):
        time_sum=0
        temp_sum=0
        for k in range(len(temperatures_for_all_heating_rates_interval_separated)):
            time_sum = time_sum + times_for_all_heating_rates_interval_separated[k][i][-1]
            temp_sum = temp_sum + temperatures_for_all_heating_rates_interval_separated[k][i][-1]
        average_time_at_claculation_point.append(time_sum/len(times_for_all_heating_rates_interval_separated))
        average_temperature_at_claculation_point.append(temp_sum/len(temperatures_for_all_heating_rates_interval_separated))
    
    print("Done ! It took ", t_stop-t_start,"s for the BFGS optimization to be performed")
    return conv_list, Ea_calculated, average_time_at_claculation_point, average_temperature_at_claculation_point



def optimize_CG(initial_guess, min_conv, max_conv, number_of_points, conv_lists, time_lists, temperature_lists):
    """
    Find the energy of activation for multiple conversions by using the advanced vyazovkin method with a conjugate gradient algorithm to find the optimal Ea.

    Parameters
    ----------
    initial_guess : Float
        Initial guess of the activation energy for the first desired conversion.
    min_conv : Float
        Minimum conversion at which the activation energy will be computed.
    max_conv : Float
        Maximum conversion at which the activation energy will be computed.
    number_of_points : Float
        Number of conversion points at which the activation energy will be computed.
    conv_lists : List
        List containing a list of evolution of conversion for various heating rates. Size (n,m) with n the number of heating rates and m the number of conversion points.
    time_lists : List
        List containing a list of evolution of time for various heating rates. Size (n,m) with n the number of heating rates and m the number of time points.
    temperature_lists : List
        List containing a list of evolution of temperature for various heating rates. Size (n,m) with n the number of heating rates and m the number of temperature points.

    Returns
    -------
    conv_list : List
        List containing conversion points at which activation energy was assessed.
    Ea_calculated : List
        List containing assessed activation energy.
    average_time_at_claculation_point : List
        List containing the average time of all heating programs at the conversion points at which activation energy was assessed.
    average_temperature_at_claculation_point : List
        List containing the average temperature of all heating programs at the conversion points at which activation energy was assessed.

    """
    t_start = time.process_time()
    check_boundaries(min_conv, max_conv,conv_lists)
    print("Progress of optimization with CG method:")
    conv_list = np.linspace(min_conv, max_conv, number_of_points) #Creates the list of conversion points at which the activation energy will be computed
    Ea_calculated = []
    average_time_at_claculation_point=[]
    average_temperature_at_claculation_point=[]
    
    #Data is extracted for all small intervals of extent that will then be used to perform the isoconversional analysis
    conversions_for_all_heating_rates_interval_separated, times_for_all_heating_rates_interval_separated, temperatures_for_all_heating_rates_interval_separated = get_data_in_interval(conv_lists, time_lists, temperature_lists, conv_list)
    
    
    for i in range(len(conv_list)): #for each conversion point at which the activation energy will be computed
        times = [heating_rate[i] for heating_rate in times_for_all_heating_rates_interval_separated] #list of times on the interval of interest for all heating rates
        temperatures = [heating_rate[i] for heating_rate in temperatures_for_all_heating_rates_interval_separated] #list of temperatures on the interval of interest for all heating rates
        
        result = scipy.optimize.minimize(func_to_minimize, x0=initial_guess, args=(
            times, temperatures), method='CG', options={'gtol': 1e-15}) #Ea found for the interval analyzed

        Ea_calculated.append(float(result.x))
        initial_guess = float(result.x) #Activation energy on next interval should be quite close from the one found at this interval, hence a new initial guess to speed up the optimization process

        printProgressBar(i, len(conv_list), length=50) #Display a progress bar in the console
    t_stop = time.process_time()
    
    #Extraction of average times and temperatures at extents analyzed
    for i in range(len(temperatures_for_all_heating_rates_interval_separated[0])):
        time_sum=0
        temp_sum=0
        for k in range(len(temperatures_for_all_heating_rates_interval_separated)):
            time_sum = time_sum + times_for_all_heating_rates_interval_separated[k][i][-1] #last element is taken because interval is [extent_for_analysis-delta_extent,extent_for_analysis] (see definition in article)
            temp_sum = temp_sum + temperatures_for_all_heating_rates_interval_separated[k][i][-1] #last element is taken because interval is [extent_for_analysis-delta_extent,extent_for_analysis] (see definition in article)
        average_time_at_claculation_point.append(time_sum/len(times_for_all_heating_rates_interval_separated))
        average_temperature_at_claculation_point.append(temp_sum/len(temperatures_for_all_heating_rates_interval_separated))
                       
                       
    
    print("Done ! It took ", t_stop-t_start,"s for the CG optimization to be performed")
    return conv_list, Ea_calculated, average_time_at_claculation_point, average_temperature_at_claculation_point





if __name__ == "__main__":

    import two_parallel_first_order_reaction_modelling as tpform
    import interpolation

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
                   np.linspace(293, 1093, number_of_points)]  # 20°C/min
    conversions = []
    Ea_list = []
    

    ax = plt.axes()
    fig, ax1 = plt.subplots(num=3)
    ax2 = ax1.twinx()
    ax2.set_ylabel('conv')

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
        plt.figure(3)
        ax1.plot(temperature[i], Ea_list[i],
                 label="Ea for "+str(5*i+5)+"°C/min")
        ax2.plot(temperature[i], conversions[i],
                 label="conversion for "+str(5*i+5)+"°C/min")
        ax1.legend()
        ax2.legend()
        
    
    

    
    new_conversions, new_times, new_temperatures = interpolation.linear_interpolation(conversions, time_lists, temperature, 20000)
    new_conversions2, new_times2, new_temperatures2 = interpolation.cubic_spline(conversions, time_lists, temperature, 20000)
    
    
    
        
        
    
    plt.figure(1)
    plt.legend()

    # t1_start = time.process_time()
    # conv_to_plot, Ea_to_plot = optimize_raw(
    #     0.001, 0.99, 2500, new_conversions, new_times, new_temperatures)
    # t1_stop = time.process_time()


    # t2_start =time.process_time()
    # conv_to_plot2, Ea_to_plot2 = optimize_raw(
    #     0.001, 0.99, 2500, new_conversions2, new_times2, new_temperatures2)
    # t2_stop =time.process_time()


    t3_start = time.process_time()
    conv_to_plot3, Ea_to_plot3, average_time3, average_temp3 = optimize_CG(
        50000, 0.01, 0.99, 500, new_conversions2, new_times2, new_temperatures2)
    t3_stop = time.process_time()
    
    t4_start = time.process_time()
    conv_to_plot4, Ea_to_plot4, average_time_at_claculation_point4, average_temperature_at_claculation_point4 = optimize_BFGS(
        50000, 0.001, 0.99, 5000, new_conversions2, new_times2, new_temperatures2)
    t4_stop = time.process_time()

    plt.figure(2)

    # plt.plot(conv_to_plot, Ea_to_plot, label="Vyazovkin modified integral (linear interpolation)")
    # plt.plot(conv_to_plot2,Ea_to_plot2,label="Vyazovkin modified integral (cubic spline)")
    plt.plot(conv_to_plot3, Ea_to_plot3, label="Vyazovkin modified integral (cubic spline) + CG")
    plt.plot(conv_to_plot4, Ea_to_plot4, label="Vyazovkin modified integral (cubic spline) + BFGS")
    plt.legend()
    
    average_temp=[]
    for i in range (len(temperature[0])):
        sum_temp=0
        for k in range(len(temperature)):
            sum_temp = sum_temp + temperature[k][i]
        average_temp.append(sum_temp/len(temperature))
            

    alpha, alpha1, alpha2, rate, rate1, rate2 = tpform.compute_alpha(time_list, average_temp, 1.666e8, 80000, 1.666e13, 120000)
    Ea_theory=tpform.global_energy_of_activation(average_temp, alpha1, alpha2, 1.666e8, 80000, 1.666e13, 120000)

    plt.figure(3)
    ax1.plot(average_temp3, Ea_to_plot3,
             label="Ea measured")
    ax1.plot(average_temp, Ea_theory,
             label="Ea theory")
    ax1.legend()
    

    # print("Elapsed time during the linear interpolation:", t1_stop-t1_start)
    # print("Elapsed time during the cubic spline in seconds:", t2_stop-t2_start)
    print("Elapsed time during the cubic spline + CG method in seconds:", t3_stop-t3_start)
    print("Elapsed time during the cubic spline + BFGS method in seconds:", t4_stop-t4_start)
