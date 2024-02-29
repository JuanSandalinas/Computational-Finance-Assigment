"""
Black scholes model
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import norm

class Black_scholes2():
    """
    Black_scholes model class.
    Inputs:
        -   S: stock price
        -   r: risk-free interest rate
        -   vol: volatility fo the stock % in decimals
        -   T: Time period
        -   N: Number of steps/intervals
        -   K: Strike price
        -   auto: Compute euler,exact method and values. True as default
    """

    def __init__(self, S,r,vol, T, K):

        self.S = S
        self.r = r
        self.vol = vol
        self.T = T
        self.N = N
        self.K = K
        self.taos = self.T - np.arange(0,self.N+1)*self.dt


    def exact_method(self,N):
        """
        Stocks price of each interval N in period T 
        using the exact solution of Black scholes
        """
        self.dt = self.T/N
        self.N = N
        ex_St= np.zeros(self.N+1)
        ex_St[0] = self.S

        #### Begin Pre-computations
        pre_1 = (self.r-(1/2)*self.vol**2)*self.dt
        pre_2 = self.vol*np.sqrt(self.dt)
        ###### End Pre-computations

        for m in range(1,self.N+1):
            Z_m = np.random.normal(0,1,1)
            ex_St[m] = ex_St[m-1]*np.exp(pre_1 + pre_2*(Z_m))
            S_t = ex_St[m]
        
        self.ex_St = ex_St


    def euler_method(self,N):
        """
        Stocks price of each interval N in period T 
        using the euler approximation solution of Black scholes
        """
        self.dt = self.T/N
        self.N = N
        eu_St = np.zeros(self.N+1)
        eu_St[0] = self.S

        #### Begin Pre-computations        
        pre_1 = self.r*self.dt
        pre_2 = self.vol*np.sqrt(self.dt)
        #### End Pre-computations

        for m in range(1,self.N+1):
            Z_m = np.random.normal(0,1,1)
            eu_St[m] = eu_St[m-1] + eu_St[m-1]*pre_1+ eu_St[m-1]*pre_2*Z_m
        self.eu_St = eu_St

    def mc_euler_european(self,n_samples):
        """
        Does MonteCarlo method for H with euler method
        Inputs:
            - n_samples: Number of samples to draw
        """


        mc_ST = np.zeros(n_samples)
        mc_Hi = np.zeros(n_samples)

        #### Begin Pre-computations        
        pre_1 = np.exp((self.r - 0.5*self.vol**2)*self.T)
        pre_2 = self.vol*np.sqrt(self.T)
        #### End Pre-computations

        for m in range(0,n_samples):
            Z = np.random.normal(0,1,1)
            mc_ST[m] = self.S*pre_1*np.exp(pre_2*Z)
            mc_Hi[m] = max(mc_ST[m]-self.K,0)
        
        self.mc_Hi = mc_Hi
        self.mc_ST  = mc_ST
        self.mc_H = np.exp(-self.r*self.T)*(np.sum(mc_Hi)/n_samples)

def milstein_scheme(self, iota,k,sigma,p, dt, n_steps = 251, mode = "arithmetic"):
    """
    Does listen scheme for Heston model:
    Inputs:
        -    Iota: long-term variance
        -    k: Reversion rate
        -    sigma: vol of vol
        -    p: correlation
        -    dt: time step
        -    mode: If using arithmetic or geometric
    """
    self.mi_St = np.zeros(n_steps)
    self.mi_Vt = np.zeros(n_steps)
    self.mi_Vt[0] = 0
    self.mi_St[0] = self.S

    ### Stock simulation
    
    for m in range(1,n_steps+1):

        Z_1 = np.random.normal(0,1,1)
        Z_2 = np.random.normal(0,1,1)

        Zv = Z_1
        Zs = p*Z_1 + np.sqrt(1-p**2)*Z_2
        self.mi_Vt[m] = self.mi_Vt[m-1] + k*(iota - self.mi_Vt[m-1])*dt +sigma*np.sqrt(self.mi_Vt[m-1]*dt)*Zv + (1/4)*(sigma**2)*dt*(Zv-1)
        
        if mode == "arithmetic":
            self.mi_St[m] = self.mi_St[m-1] + self.r*self.mi_St[m-1]*dt + np.sqrt(self.mi_Vt[m-1]*dt)*self.mi_St[m-1]*Zs + 1/2 * self.mi_Vt[m-1]*self.mi_St[m-1]*dt*(Zs**2 - 1)
        elif mode == "geometric":
            self.mi_St[m] = self.mi_St[m-1]*np.exp((r-1/2 *self.mi_Vt[m-1])*dt + np.sqrt(self.mi_Vt[m-1]*dt)*Z_s)

def option_value_milstein_scheme(self,iota,k,sigma,p, dt, n_steps = 251, mode = "arithmetic"):
    
    """Still need to implement
    """
    



     


if __name__ == "__main__":
    T = 1
    S = 100
    K = 99
    r = 0.06
    vol = 0.2
    N = 100
    black_scholes = Black_scholes(S,r,vol,T,N, K, auto = False)
    H = []
    n_samples_list = [10,20,50,100,200,300]
    for n_samples in n_samples_list:
        black_scholes.mc_european(n_samples) 
        H += [black_scholes.mc_H]
    plt.plot(H)
    plt.show()
