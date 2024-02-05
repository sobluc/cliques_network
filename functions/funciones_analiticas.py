# -----------------------------------------------------------------------------------------------

# Archivo con funciones que calculan las expresiones analíticas para las redes.
#
# fn es la distribución de probabilidad para el tamaño de cliques (debe estar truncada y f0 = 0).
# Q es la cantidad de cliques en el sistema.
# gamma es la probabilidad de que un nodo tenga un enlace a un nodo en otro clique.

# -----------------------------------------------------------------------------------------------

import numpy as np
import math

def suma_n_de_fn(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum +=  fn(n)
        n += 1
    return sum


def suma_n_de_nfn(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += n * fn(n)
        n += 1
    return sum

def suma_n_de_n2fn(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += (n ** 2) * fn(n)
        n += 1
    return sum

def suma_n_de_n3fn(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += (n ** 3) * fn(n)
        n += 1
    return sum

def suma_n_de_n4fn(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += (n ** 4) * fn(n)
        n += 1
    return sum

def cantidad_de_nodos(fn, Q, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += n * fn(n)
        n += 1
    return sum * Q   

def distribucion_de_grado(g, fn, gamma, n_max):
    A = gamma * g * fn(g) / suma_n_de_nfn(fn, n_max)
    B = (1-gamma) * (g + 1) * fn(g + 1) / suma_n_de_nfn(fn, n_max)
    return A + B

def g_k(k, fn, gamma, n_max):
    if isinstance(k, np.ndarray):
        ret_vec = gamma * k * fn(k) / suma_n_de_nfn(fn, n_max) + (1 - gamma) * (k+1) * fn(k+1) / suma_n_de_nfn(fn, n_max)
        ret_vec[-1] = gamma * k[-1] * fn(k[-1]) / suma_n_de_nfn(fn, n_max)
        return ret_vec

    else:
        if k == n_max:
            return gamma * k * fn(k) / suma_n_de_nfn(fn, n_max) + (1 - gamma) * (k+1) * fn(k+1) / suma_n_de_nfn(fn, n_max)
        else:
            return gamma * k * fn(k) / suma_n_de_nfn(fn, n_max)

def coef_clustering_local(red, nodo):
    clique = red.clique_list[nodo[0]]
    clique_size = len(clique.nodes())    
    if red.nodes()[nodo]["ext_connection"]:
        return (clique_size - 2) / clique_size
    else:
        return 1.0

def coef_clustering_medio(fn, gamma, n_max):
    sum = 0
    n = 3
    while n <= n_max:
        sum += (n - 2 * gamma) * fn(n)
        n += 1
         
    return sum / suma_n_de_nfn(fn, n_max)


def coef_clustering_global(fn, gamma, n_max):
    sum1 = 0
    n = 1
    while n <= n_max:
        sum1 += n * (n-1) * (n-2) * fn(n)
        n += 1

    sum2 = 0
    n = 1
    while n <= n_max:
        sum2 += n * (n-1) * fn(n)
        n += 1

    return 0.5 * sum1 /(0.5 * sum1 + gamma * sum2)



def G_0_exponencial(x,alpha ,gamma):
    A = (1-np.exp(-alpha)) * (1 + gamma * (x - 1))
    B = 1 - np.exp(-alpha) * (1+gamma*(x-1))
    return A / B

def G_0_Kronecker(x, m ,gamma):
    return (1 + gamma * (x-1)) ** m


def G_0_constante(x, N_max ,gamma):
    return (((1+gamma*(x-1))**(N_max + 1)-1)/(gamma*(x-1)) - 1) / N_max


def fraccion_comp_gigante_exponencial(gamma, alpha):
    if -4*gamma + 4 * gamma * np.exp(alpha) + gamma**2 < 0:
        return 1-G_0_exponencial(1.0, alpha , gamma) 

    U = (-2 + 2 * np.exp(alpha) + gamma - np.sqrt(-4* gamma + 4* gamma *np.exp(alpha) + gamma**2))/(2 * gamma)

    if U < 1 and U >= 0:
        return 1-G_0_exponencial(U, alpha , gamma)
    else:
        return 1-G_0_exponencial(1.0, alpha , gamma)         

def fraccion_comp_gigante_Kronecker_3(gamma):

    U = (gamma - 1)**2 / gamma**2

    if U < 1 and U >= 0.0:
        return 1-G_0_Kronecker(U, 3 , gamma)
    else:
        return 1-G_0_Kronecker(1.0, 3 , gamma)

def fraccion_comp_gigante_Kronecker_4(gamma):

    if 4 - 3 * gamma < 0:
        return 1-G_0_Kronecker(1.0, 4 , gamma) 

    U = 1 - 3 / ( 2 * gamma) + (np.sqrt((4*(gamma**3) - 3 * (gamma ** 4)))) / (2 * gamma ** 3)

    if U < 1 and U >= 0.0:
        return 1 - G_0_Kronecker(U, 4 , gamma)
    else:
        return 1 - G_0_Kronecker(1.0, 4 , gamma)


def fraccion_comp_gigante_Kronecker_5(gamma):
    if 27 - 40 * gamma + 16 * gamma**2 < 0:
        return 1-G_0_Kronecker(1.0, 5 , gamma) 


    raiz = np.sqrt(27*gamma**16 - 40*gamma**17 + 16* gamma**18)
    raiz_cubica = (27*gamma**8 - 20 * gamma**9+ 3*3**(1/2)*raiz)**(1/3)
    A = - (4*gamma**3 - 3*gamma**4) / (3 * gamma **4)
    B= -((2 * 2**(1/3)) * gamma ** 2)/(3*raiz_cubica)
    C = raiz_cubica / (3*2**(1/3) * gamma**4)

    U = A + B + C 

    if U < 1 and U >= 0.0:
        return 1-G_0_Kronecker(U, 5 , gamma)
    else:
        return 1-G_0_Kronecker(1.0, 5 , gamma)


def cant_enlaces_clique_dividido_Q(k1, k2, n_moño ,  n, gamma , fn):
    
    if k1 == n and k2 == n:
        return  n_moño * (n_moño-1) * math.comb(n, n_moño) * (gamma ** n_moño) * ((1- gamma)**(n-n_moño)) * fn(n) / 2
    
    elif k1 == n-1 and k2 == n-1:
        return  (n - n_moño) * (n - n_moño - 1) * math.comb(n, n_moño) * (gamma ** n_moño) * ((1- gamma)**(n-n_moño)) * fn(n) / 2

    else:
        return (n * (n-1) - n_moño * (n_moño-1) - (n - n_moño) * (n - n_moño - 1)) * math.comb(n, n_moño) * (gamma ** n_moño) * ((1- gamma)**(n-n_moño)) * fn(n) / 2




def sumas_homofilia(fn, gamma, n_max):

    Q_e = 2 / (suma_n_de_n2fn(fn, n_max) - (1-gamma) * suma_n_de_nfn(fn, n_max))
    
    # sumo enlaces de largo alcance
    suma_1 = gamma * (suma_n_de_n2fn(fn, n_max) ** 2) / (suma_n_de_nfn(fn, n_max) * 2) 
    suma_2 = gamma  * suma_n_de_n2fn(fn, n_max)
    suma_3 = gamma  * suma_n_de_n3fn(fn, n_max)


    # sumo enlaces de cliques
    n = 1 
    while n <= n_max:
        n_moño = 0
        
        while n_moño <= n:
            suma_1 += n ** 2 * cant_enlaces_clique_dividido_Q(n, n, n_moño ,  n, gamma , fn)
            suma_1 += (n-1) ** 2 * cant_enlaces_clique_dividido_Q(n-1, n-1, n_moño ,  n, gamma , fn)
            suma_1 +=  n * (n-1) * cant_enlaces_clique_dividido_Q(n, n-1, n_moño ,  n, gamma , fn)        
            
            suma_2 += 2 *  n * cant_enlaces_clique_dividido_Q(n, n, n_moño ,  n, gamma , fn)
            suma_2 += (2 * n - 2) * cant_enlaces_clique_dividido_Q(n-1, n-1, n_moño ,  n, gamma , fn)
            suma_2 +=  (2 * n - 1) * cant_enlaces_clique_dividido_Q(n, n-1, n_moño ,  n, gamma , fn)


            suma_3 += 2 * n ** 2 * cant_enlaces_clique_dividido_Q(n, n, n_moño ,  n, gamma , fn)
            suma_3 += 2 * (n-1) ** 2 * cant_enlaces_clique_dividido_Q(n-1, n-1, n_moño ,  n, gamma , fn)
            suma_3 +=  (n ** 2 + (n-1)**2) * cant_enlaces_clique_dividido_Q(n, n-1, n_moño ,  n, gamma , fn)
 
            n_moño += 1
        
        n += 1
    

    suma_1 *= Q_e
    suma_2 *= Q_e
    suma_3 *= Q_e

    return suma_1, suma_2, suma_3


################################################################################################

def mean_val_1(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += (n ** 4 - n ** 3) * fn(n)
        n += 1
    return sum

def mean_val_2(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += (n * (n -1) ** 3) * fn(n)
        n += 1
    return sum

def mean_val_3(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += (n**2 * (n -1)**2) * fn(n)
        n += 1
    return sum


def mean_val_4(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += (n **2 *  (n -1) ** 2) * fn(n)
        n += 1
    return sum

def mean_val_5(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += (n **2 *  (n -1) ) * fn(n)
        n += 1
    return sum

def mean_val_6(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += (n  *  (n -1)**2 ) * fn(n)
        n += 1
    return sum

def mean_val_7(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += ((2  *  n - 1) * n * (n-1) ) * fn(n)
        n += 1
    return sum


def mean_val_8(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += ((2  *  n - 1) * n * (n-1) ) * fn(n)
        n += 1
    return sum

def mean_val_9(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += ( (n ** 2 +  (n-1)**2)*n*(n-1) ) * fn(n)
        n += 1
    return sum

def mean_val_10(fn, n_max):
    sum = 0
    n = 1
    while n <= n_max:
        sum += ( (n ** 2 +  (n-1)**2) * n * (n-1) ) * fn(n)
        n += 1
    return sum


def sumas_homofilia_2(fn, gamma, n_max):

    Q_e = 2 / (suma_n_de_n2fn(fn, n_max) - (1-gamma) * suma_n_de_nfn(fn, n_max))
    
    # sumo enlaces de largo alcance
    suma_1 = gamma * (suma_n_de_n2fn(fn, n_max) ** 2) / (suma_n_de_nfn(fn, n_max) * 2) 
    suma_2 = gamma  * suma_n_de_n2fn(fn, n_max)
    suma_3 = gamma  * suma_n_de_n3fn(fn, n_max)

    v1 = suma_n_de_nfn(fn, n_max)
    v2 = suma_n_de_n2fn(fn, n_max)
    v3 = suma_n_de_n3fn(fn, n_max)
    v4 = suma_n_de_n4fn(fn, n_max)

    suma_1 += (gamma**2 / 2) * mean_val_1(fn, n_max) + ((1- gamma)**2 / 2) * mean_val_2(fn, n_max) + mean_val_3(fn, n_max) / 2 - (gamma**2 / 2) * mean_val_4(fn, n_max) - ((1- gamma)**2 / 2) * mean_val_4(fn, n_max)
    suma_2 += gamma**2 * mean_val_5(fn, n_max) + (1-gamma)**2 * mean_val_6(fn, n_max) + mean_val_7(fn, n_max) / 2 - (gamma**2 / 2) * mean_val_8(fn, n_max) - ((1- gamma)**2 / 2) * mean_val_8(fn, n_max)
    suma_3 += gamma**2 * mean_val_1(fn, n_max) + (1- gamma)**2  * mean_val_2(fn, n_max) + mean_val_9(fn, n_max) / 2 - (gamma**2 / 2) * mean_val_10(fn, n_max) - ((1- gamma)**2 / 2) * mean_val_10(fn, n_max)

    suma_1 *= Q_e
    suma_2 *= Q_e
    suma_3 *= Q_e

    return suma_1, suma_2, suma_3