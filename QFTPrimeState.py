#%matplotlib inline

# general imports
import random
import math
import sys
import time

# libraries needed for the QFT
from numpy import cumsum, sum, searchsorted
from numpy.random import rand
import matplotlib.pyplot as plt
from operator import itemgetter

# libraries needed for complex numbers, fractions and gcd
from fractions import Fraction, gcd
from decimal import Decimal
import cmath


def sieveOfEratosthenes(n):
    """
    sieveOfEratosthenes(n): return the list of the primes < n.
    Code from: <dickinsm@gmail.com>, Nov 30 2006
    http://groups.google.com/group/comp.lang.python/msg/f1f10ced88c68c2d
    """
    
    sieve = range(3, n, 2)
    top = len(sieve)
    for si in sieve:
        if si:
            bottom = (si*si - 3) / 2
            if bottom >= top:
                break
            sieve[bottom::si] = [0] * -((bottom - top) / si) 

    return [2] + [el for el in sieve if el]


def modularPrimeCounting(list_Primes, mod, a):
    """ computes the modular Prime Counting Function for primes in list_Pimes equals to a modulo mod"""
    counter_a = 0
    for prime in list_Primes:
        if( prime % mod == a ):
            counter_a += 1
            #print prime
    #print ""
    return counter_a


def qftProbability(N,list_Primes,k):
    """computes the QFT probabilty P(k) = |f(k)|^2"""
    
    summation = 0
    for prime in list_Primes:
        summation += cmath.exp(2*math.pi*1j*prime*k/N)
    return abs(summation)**2


def qftPeaks(N,list_Primes, pi_N, t1):
    """computes the QFT peaks"""
    
    weights_tuple = []
    weights_tuple_normalized = []
    for k in range(0,N):
        if k % 5000 == 0:
            t2 = time.time()
            print k, "([{:.3f} s])".format(t2 - t1)
            t1 = t2
            sys.stdout.flush()
        weight = qftProbability(N, list_Primes, k)
        weights_tuple.append([k,weight])
        weights_tuple_normalized.append([k, weight/ float(N*pi_N)])    
    return weights_tuple, weights_tuple_normalized, t1


def QFT_Simulation(N, qubits, num_peaks):
    """
    - simulates the Quantum Fourier transform of the Prime state for primes less than N
    - if qubits is true, then it does so with the first power of 2 greater than N
    - num_peaks is the number of highests peaks are going to be shown numerically
    """

    
    t0 = time.time()
    t1 = time.time() 
    
    if qubits:
        n = int(math.ceil(math.log(N,2))) #qubits
        N = pow(2,n)
    list_Primes = sieveOfEratosthenes(N)
    pi_N = len(list_Primes)
    print "- Number of Primes, PI(", N, ") = ", pi_N
    print "- Normalization (divide by): N*PI(N) =", N*pi_N
    
    a, b = modularPrimeCounting(list_Primes, 3, 1), modularPrimeCounting(list_Primes, 3, 2)
    print "- Delta_{3;2,1} bias =", b, "-", a, "=", b-a
    
    a, b = modularPrimeCounting(list_Primes, 4, 1), modularPrimeCounting(list_Primes, 4, 3)
    print "- Delta_{4;3,1} Chebyshev bias =", b, "-", a, "=", b-a
    
    a, b = modularPrimeCounting(list_Primes, 6, 1), modularPrimeCounting(list_Primes, 6, 5)
    print "- Delta_{6;5,1} bias =", b, "-", a, "=", b-a

    maxWeight = []
    for i in range(0, num_peaks):
        maxWeight.append([0,0])
    
    weights_tuple = []
    t2 = time.time()
    print "\nComputing QFT: ({:.3f} s)".format(t2 - t0)
    t1 = t2
    sys.stdout.flush()
    weights_tuple, weights_tuple_normalized, t1 = qftPeaks(N,list_Primes, pi_N, t1)
    
    maxWeight_sorted = sorted( weights_tuple, key=itemgetter(1) )[N-num_peaks:]
    maxWeight_normalized_sorted = sorted( weights_tuple_normalized, key=itemgetter(1) )[N-num_peaks:]
    
    highest_peak = pi_N/float(N)
    print "\n", num_peaks, "highests peaks:"
    print "\tk \tvalue\tvalue_norm\tproportion\t(position)"
    print "\t", "-"*60
    for i in reversed(range(0,num_peaks)):
        print "\t{}\t{:.1f}\t{:.10f}\t{:.3f}\t\t({:.3f})".format(maxWeight_sorted[i][0], maxWeight_sorted[i][1], maxWeight_normalized_sorted[i][1], maxWeight_normalized_sorted[i][1] / highest_peak, maxWeight_sorted[i][0]/float(N))
    
    fraction = 8
    t2 = time.time()
    print "\nComputing peaks (normalized) of fraction {}: ({:.3f} s)".format(fraction, t2 - t0)
    t1 = t2
    sys.stdout.flush()
    for k in range(0,fraction):
        print "\tP({}*N/{}) = {}\n\t\t  ({})".format(k, fraction, weights_tuple[k*N/fraction][1], weights_tuple_normalized[k*N/fraction][1])
     
    fraction = 6
    t2 = time.time()
    print "\nComputing peaks (normalized) of fraction {}: ({:.3f} s)".format(fraction, t2 - t0)
    t1 = t2
    sys.stdout.flush()
    for k in range(0,fraction):
        print "\tP({}*N/{}) = {}\n\t\t  ({})".format(k, fraction, weights_tuple[k*N/fraction][1], weights_tuple_normalized[k*N/fraction][1])    
    
    print ""
    
    plt.close()
    plt.plot(zip(*weights_tuple_normalized)[1], 'b')

    limit_X = N / 20
    limit_Y = maxWeight_normalized_sorted[-1][1] / 20
    plt.axis( [ -limit_X , N + limit_X , -limit_Y , maxWeight_normalized_sorted[-1][1] + limit_Y ], figsize=(8, 10), dpi=500)
    plt.xlabel("Register")
    plt.ylabel("Probability")
    plt.show()
    
    t2 = time.time()
    print "\nDONE! ({:.3f} s)".format(t2 - t0)
    sys.stdout.flush()


    return

