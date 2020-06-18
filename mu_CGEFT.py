
import sys
import matplotlib.pyplot as plt
import math
import numpy as np




def mu_CGEFT(a,k):

                                 
                      
    Omega0=0.25
    OmegaLambda=0.75
    H_over_c=100./299792.458
    #Hubble100=0.70
    m0=1.
    M=math.pow(H_over_c,2./3.)         

   
    c3=1./6./math.sqrt(6.*OmegaLambda)
    ksi=math.sqrt(6.*OmegaLambda)


    H=H_over_c*math.sqrt((Omega0/a/a/a+math.sqrt(Omega0*Omega0/math.pow(a,6.)+4.*OmegaLambda))/2.)
    dHdt=(math.pow(H_over_c,4.)*OmegaLambda/H/H-H*H-H_over_c*H_over_c*Omega0/a/a/a/2.)/(1.+math.pow(H_over_c/H,4.)*OmegaLambda)
    ddHddt=(4.*dHdt*dHdt*math.pow(H_over_c/H,4.)*OmegaLambda/H-2.*dHdt*math.pow(H_over_c/H,4.)*H*OmegaLambda-2.*H*dHdt+3.*H*H_over_c*H_over_c*Omega0/a/a/a/2.)/(1.+
math.pow(H_over_c/H,4.)*OmegaLambda)


    dphidt=ksi*H_over_c*H_over_c/H
    ddphiddt=-ksi*H_over_c*H_over_c*dHdt/H/H

  
    OMEGA=1.
    dOMEGAdt=0.
    c=dphidt*dphidt*(-6.*H*c3*dphidt/M/M/M+2.*c3*ddphiddt/M/M/M-1.)/2.
    M_hat=0.
    #M1_bar=math.pow(2.*c3*dphidt*dphidt*dphidt/M/M/M,1./3.)
    M1_bar=math.pow(2.*c3,1./3.)*dphidt/M 
    M2_bar=0.
    M3_bar=0.
    m2=0.
    dM_hatdt=0.
    dM3_bardt=0.
    pdM_hat2dt=0.
    pdM1_bar3dt=0.
    dR0dt=0.


  
    A1=2.*m0*m0*OMEGA+4.*M_hat*M_hat
    A2=-m0*m0*dOMEGAdt-M1_bar*M1_bar*M1_bar+2.*H*M3_bar*M3_bar+4.*H*M_hat*M_hat
    A3=-8.*m2*m2
    B1=-1.-2.*M_hat*M_hat/m0/m0/OMEGA
    B2=1.
    #B3=-dOMEGAdt/OMEGA+M3_bar*M3_bar*(H+2.*dM3_bardt/M3_bar)/m0/m0/OMEGA  
    B3=0.
    C1=m0*m0*dOMEGAdt+2.*H*M_hat*M_hat+4.*M_hat*dM_hatdt
    C2=-m0*m0*dOMEGAdt/2.-M1_bar*M1_bar*M1_bar/2.-3.*H*M2_bar*M2_bar/2.-H*M3_bar*M3_bar/2.+2.*H*M_hat*M_hat
    C3=c-H*M1_bar*M1_bar*M1_bar/2.-pdM1_bar3dt/2.+(k*k/2./a/a-3.*dHdt)*M2_bar*M2_bar+(k*k/2./a/a-dHdt)*M3_bar*M3_bar+2.*(H*H+dHdt)*M_hat*M_hat+2.*H*pdM_hat2dt
    CPI=m0*m0*dOMEGAdt*dR0dt/4.-3.*c*dHdt+3.*(3.*H*dHdt+ddHddt)*M1_bar*M1_bar*M1_bar/2.+3.*dHdt*pdM1_bar3dt/2.+\
9.*dHdt*dHdt*M2_bar*M2_bar/2.+3.*dHdt*dHdt*M3_bar*M3_bar/2.  
  
  # /* pdM_hat2dt is the partial derivative of M_hat^2 with respect to t, pdM1_bar3dt is the partial derivative of M1_bar^3 with respect to t       */
  
    f1=B2*C3-C1*B3
    f2=B2*CPI
    f3=A1*(B3*C2-B1*C3)+A2*(B1*C1-B2*C2)+A3*(B2*C3-B3*C1)
    f4=(A3*B2-A1*B1)*CPI

    return 2.*(f1+f2*a*a/k/k)/(f3+f4*a*a/k/k)


mu=np.zeros((150,200))

for hhloga in range(0,150):
    for hhlogk in range(0,200):
        
        a=math.pow(10.,-hhloga/100.)
        k=math.pow(10.,-hhlogk/100.)
        
        mu[hhloga][hhlogk]=mu_CGEFT(a,k)


      


im=plt.imshow(mu,origin='lower')
plt.xlabel("-100*logk")
plt.ylabel("-100*loga")
plt.title('CG_EFT')
plt.colorbar(im)
plt.show()


