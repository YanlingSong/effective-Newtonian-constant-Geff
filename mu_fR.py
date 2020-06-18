import sys
import matplotlib.pyplot as plt
import math
import numpy as np



def mu_fr(a, k):

  
    #Hubble100=0.70                                                                    
    Omega0=0.25
    OmegaLambda=0.75mu_CGEFT
    H_over_c=100./299792.458
    FR0=1.e-4


    B1   = Omega0 / pow(a, 3.) + 4. * OmegaLambda
    B2   = Omega0 + 4. * OmegaLambda
    emme = 0.5 * H_over_c * H_over_c * pow(B1, 3.) / (B2 * B2 * FR0) 
    return 1. + k * k / 3. / (k * k + a * a * emme)
#this function of mu is from arXiv:1909.06282   and the code is from PINOCCHIO-4.1,from the MG branch of https://github.com/pigimonaco/Pinocchio, 
#in the file cosmo.c 
 



mu=np.zeros((150,200))

for hhloga in range(0,150):
    for hhlogk in range(0,200):
        
        a=math.pow(10.,-hhloga/100.)
        k=math.pow(10.,-hhlogk/100.)
        mu[hhloga][hhlogk]=mu_fr(a,k)



im=plt.imshow(mu,origin='lower')
plt.xlabel("-100*logk")
plt.ylabel("-100*loga")
plt.title('fR')
plt.colorbar(im)
plt.show()











