# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 11:00:47 2014

Accompanying python code for the Sunspot Number Tutorial.

@author: Dr. Lisa Winter, lwinter@aer.com
"""

import numpy as np # import the NumPy or Numerical Python librar
import datetime # a library that makes it easy to work with dates
import matplotlib.pylab as plt # plotting library

# Define a function to obtain the smoothed sunspot number
def get_sunspotnumber():
    '''
    Get the sunspot number from the NASA one-month smoothed data.
    Pass back the information from the website.
    '''
    
    import urllib2 # import a library that allows us to read data from the web
    
    sunspot_number = urllib2.urlopen('http://solarscience.msfc.nasa.gov/greenwch/spot_num.txt')
    
    sunspot_dtype = [('year', int), ('month', int), ('SSN', float), ('DEV', float)] # the data type for our numpy array
    
    data = np.loadtxt(sunspot_number, skiprows=1, dtype=sunspot_dtype)   # loads the data on the webpage into a numpy array
            
    return data

def gauss(x, *p):
    A, mu, sigma = p  # These are input to the data: the normalization, mid-point, and width of the Gaussian
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


# Define a function for a Gaussian fit
def fit_gaussian(x, y):
    ''' Fit a Gaussian to the input x and y data.
    Return an array with the best-fit parameters to the data.
    '''
    
    from scipy.optimize import curve_fit # We will load from the scipy library a useful function for fitting a curve
        
    p0 = [100, 5, 5] # Define some initial conditions for first guesses on A, mu, and sigma
    
    coeff, var_matrix = curve_fit(gauss, x, y, p0=p0) # Fit the data
    
    print 'The best-fit values: ',
    print coeff
    
    return coeff # Pass back the results

# Define a function to convert datetime dates to a decimal of a year for the year in the Solar Cycle
def convert_to_solarcycleyear(dates, startdate):
    '''
    Convert calendar date to year in the solar cycle
    Input: an array of datetime objects and the startdate of the solar cycle
    Output: an array of the deminal year since the beginning of the solar cycle
    '''
    
    ssy = [((x - startdate).total_seconds())/(60.*60.*24.*365.2425) for x in dates]
    return np.array(ssy)

# Read in the data and assign to a variable
sunspotnumber_data = get_sunspotnumber()
# The data stores just year and month.  Convert this to an array of datetime objects.
sunspotnumber_dates = [datetime.datetime(x['year'], x['month'], 1) for x in sunspotnumber_data] # convert to datetime objects
sunspotnumber_dates = np.array(sunspotnumber_dates) # convert our Python list to a NumPy array


def initiallookathistoricaldata():

    # Let's explore the data.  First, let's find out what the size of our dataset is and what the first and last entries are.
    print 'The size of our dataset is: %d' % (len(sunspotnumber_data)) # prints the length of the dataset

    print 'The first entry is: ',
    print sunspotnumber_data[0]

    print 'The last entry is: ',
    print sunspotnumber_data[-1]  

    print sunspotnumber_dates[0] # print the first date

    plt.clf() # clear the plotting window
    plt.plot(sunspotnumber_dates, sunspotnumber_data['SSN'], 'k-') # plot the SSN versus date as a black (k) line
    plt.xlabel('Date')
    plt.ylabel('Sunspot Number')

    max_SSN_index = np.argmax(sunspotnumber_data['SSN']) 
    plt.plot([sunspotnumber_dates[max_SSN_index]], [sunspotnumber_data['SSN'][max_SSN_index]], 'ro', ms=10)
    print 'The Maximum Sunspot Number in the data was on: ',
    print sunspotnumber_dates[max_SSN_index], 
    print ' The number of sunspots recorded: ',
    print sunspotnumber_data['SSN'][max_SSN_index] 
    plt.savefig('HistoricalSSN.png') # save the figure


def plottherecentdata():
    # Filter the data to just plot the recent solar cycles
    ind_after1996 = sunspotnumber_dates > datetime.datetime(1996, 1, 1) # Create an index of which dates meet our condition (True) and which do not (False)
    print ind_after1996[0], ind_after1996[-1]  # We should see that the first date (in the 1700s!) is False and the most current is True

    plt.clf() # clear the plotting window
    plt.plot(sunspotnumber_dates[ind_after1996], sunspotnumber_data['SSN'][ind_after1996], 'k-') # Use the index to only plot more recent data
    plt.xlabel('Date')
    plt.ylabel('Sunspot Number')
    plt.savefig('RecentCycles.png')


def analyzesolarcycle(startdate=datetime.datetime(1996, 1, 1), enddate=datetime.datetime(2009, 1, 1), cycle=23):    
    # Filter the data to work with just the Solar Cycle 23 Data and Fit a Gaussian to find the solar maximum
    ind_solarcycle = (sunspotnumber_dates >= startdate) & (sunspotnumber_dates <= enddate)

    solarcycle_cycleyear = convert_to_solarcycleyear(sunspotnumber_dates[ind_solarcycle], startdate) # convert dates to solar cycle year
    solarcycle_SSN = sunspotnumber_data[ind_solarcycle]['SSN'] # create a new variable for sunspot number of just this cycle
    
    fit_parameters = fit_gaussian(solarcycle_cycleyear, solarcycle_SSN)
    A = fit_parameters[0]
    mu = fit_parameters[1]
    sigma = fit_parameters[2]

    fit_x = np.arange(0, 13, 0.25) # Create a new array that goes from 0 to 13 in steps of 0.25
    fit_y = A*np.exp(-(fit_x-mu)**2/(2.*sigma**2))

    expected = gauss(solarcycle_cycleyear, A, mu, sigma)
    ind = solarcycle_SSN > 0.
    chi_squared = np.sum((solarcycle_SSN[ind] - expected[ind])**2/np.abs(0.1*solarcycle_SSN[ind]))
    print 'Chi Squared : %.2lf/%d' % (chi_squared, len(solarcycle_SSN[ind])-3)


    plt.clf() # clear the plotting window
    plt.plot(solarcycle_cycleyear, solarcycle_SSN, 'k-', label='Data') # Plot the solar cycle data
    plt.plot(fit_x, fit_y, 'r--', label='Fit')
    plt.legend()
    plt.xlabel('Solar Cycle Year')
    plt.ylabel('Sunspot Number')
    plt.title('Solar Cycle %d' % (cycle))
    plt.savefig('SunspotCycle%d.png' % (cycle))
    
    days_from_startofcycle = mu*365.2425 # convert from years to days
    solarmax_date = startdate + datetime.timedelta(days=int(days_from_startofcycle))
    print 'The Estimated Solar Max Date was: %s' % (solarmax_date)
    
    peak1_ind = np.argmax(solarcycle_SSN[solarcycle_cycleyear < mu])
    peak1_date = solarmax_date = startdate + datetime.timedelta(days=int(solarcycle_cycleyear[peak1_ind]*365.2425))
    print 'Peak 1: SSN = %d, On Date = %s' % (solarcycle_SSN[peak1_ind], peak1_date)
    print solarcycle_cycleyear[peak1_ind], peak1_ind
    
    peak2_ind = np.argmax(solarcycle_SSN[solarcycle_cycleyear > mu+1./12.]) + np.argmin(np.abs(solarcycle_cycleyear - mu)) + 1
    peak2_date = solarmax_date = startdate + datetime.timedelta(days=int(solarcycle_cycleyear[peak2_ind]*365.2425))
    print 'Peak 2: SSN = %d, On Date = %s' % (solarcycle_SSN[peak2_ind], peak2_date)
    print solarcycle_cycleyear[peak2_ind], peak2_ind