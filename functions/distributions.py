from random import random
import numpy as np

def fn_exponencial(alpha):
    def exp_return(n):
        return np.exp(- alpha * n)

    return exp_return


def fn_exponencial_beta(beta):
    def exp_return(n):
        return beta * ((1 + beta)**(-n))

    return exp_return


def fn_constante(cte):

    def cte_return(n):
        return  cte + n - n   
    
    return cte_return



def fn_Kronecker(m):

    def Kronecker_return(n):
            if isinstance(n, np.ndarray):
                ret_array = np.zeros(len(n))
                ret_array[m - 1] = 1
                return ret_array
            else:
                if n == m:
                    return 1 
                else:
                    return 0 
    
    return Kronecker_return

def fn_random(n_max):

    random_array = np.zeros(n_max)
    for i in range(len(random_array)):
        random_array[i] = np.random.uniform(low = 0.0, high = 1.0)

    def func_return(n):
        if isinstance(n, np.ndarray):
            return random_array
        else:
            return  np.random.uniform(low = 0.0, high = 1.0)
    
    return func_return


# def fn_moño(fn, n_max):
#     def func_ret(n):
#         if isinstance(n, np.ndarray):
#             return (1-(1- gamma * sigma)**n) * fn(n) / sum((1-(1-gamma*sigma)**n) * f_n)

#         else:


#     return func_ret


# (1-(1-gamma*sigma)**n) * f_n / sum((1-(1-gamma*sigma)**n) * f_n)