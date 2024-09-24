import numpy as np
import math
import luescher.tools.kinematics as kin_tools

available_fit_params_list = [ 'ERE', 'ERE_delta', 'ERE_npi_Eq12']

# single_channel parametrizations
def output(ecm,ma,mb,fit_param, fit_params):
    if fit_param == 'ERE':
        return ere(ecm,ma,mb,*fit_params)
    elif fit_param == 'ERE_delta':
        return ere_delta(ecm,ma,mb,*fit_params)
    elif fit_param == 'ERE_npi_Eq12':
        return ere_npi12(ecm,*fit_params)

def error_output(ecm,ma,mb,fit_param, fit_params):
    if fit_param == 'ERE':
        return ere(ecm,ma,mb,*fit_params)
    elif fit_param == 'ERE_delta':
        vec = lambda e:np.array([e,e*delta_Sp(e,ma,mb)]) 
        return vec(ecm)
    elif fit_param == 'ERE_npi_Eq12':
        vec = lambda e: np.array([-e/fit_params[0]**2])
        return vec(ecm)


## parametrizations should also have some info about where they come from
def delta_Sp(ecm,ma,mb):
    return (ecm**2 - (ma+mb)**2 )/ (ma+mb)**2

def ere_delta(ecm,ma,mb,a,b):
    return (ecm)*(a+b*delta_Sp(ecm,ma,mb))

def ere_npi12(ecm,A):
    return ecm/A

# # fit ere expansion
# def ere(ecm,ma,mb,a,b):
#     return ((-1/a)+0.5*b*kin_tools.q2(ecm,ma,mb) ) #in units of reference mass, usually mpi

# define ere for range of n
def ere(ecm,ma,mb,*p):
    #qcotd = ere = p[0] + p[1] * x + ...
    # p[0]= -1/a , ....
    # x is q2
    qcotd = 0
    x = kin_tools.q2(ecm,ma,mb)
    for n in range(len(p)):
        qcotd += p[n] * x**n 
    return qcotd 

