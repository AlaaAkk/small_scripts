#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Libraries
import  pandas as pd
import math
import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import read , write
from numpy import array, reshape, zeros, append, arange, ones
import matplotlib
import matplotlib as mpl
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt


# In[2]:



### conversions and constants ###
pi=np.pi
cmhz=29979245800.0*2*pi #cm^-1 to Hz
amu2g=1.660539e-24
ang2bohr=1.8897261
ang2cm=1e-08


p0=1013250  # atm to g/(cm s^2)

kk=1.380649e-16 # erg/k (cm^2.g/ks^2)
k=8.617333262145e-05 # ev/k
h=6.62607015e-27  # erg.s
hb=6.582119569e-16 # eV.s
hbar=1.054571817e-27 # erg.s


# In[3]:

__all__ = ['main']

# Spliting function
def split_line(lines):
    """Split input line"""
    line_array = np.array(lines.strip().split(" "))
    line_vals = line_array[line_array != ""]
    return line_vals


# In[4]:


# reading total energies in aims
def read_energy(filename):
     energy=0
     with open(filename) as out:

        for line in out:
           if line.rfind('| Total energy of the DFT / Hartree-Fock s.c.f. calculation      : ')!=-1:
              energy = float(split_line(line)[-2]) # Periodic/cluster
              return energy


# In[5]:


# reading total energies in aims
def read_omega(filename):
    """ reads frequencies in cm-1 and returns them in hz """
    ftemp=open(filename)
    fw=np.array([float(line.split()[0]) for line in ftemp])*cmhz
    ftemp.close()
    return fw


# #### $ \frac{1}{2}(E_{MoS2}-\mu^{Bulk}_{Mo})\leq \mu_{S} \leq \mu_{S8} $

# In[6]:


def main():
    def chemical_pot_1(x,y,z):
        # x is E_MoS2
        # y is \mu_Mo
        # z is \mu_S8

        mu_S0_i=0.5*(x-y)
        mu_S0_f=z
        mu_Mo_i=x-2*mu_S0_i
        mu_Mo_f=x-2*mu_S0_f

        mu_S0=arange(mu_S0_i,mu_S0_f,0.05)
        mu_Mo0=arange(mu_Mo_f,mu_Mo_i,0.05)
        var_S=mu_S0-mu_S0_f
        var_Mo=(mu_Mo0-mu_Mo_f)/2
        yield np.array(mu_S0)
        yield np.array(mu_Mo0)
        yield var_S
        yield var_Mo
        yield mu_S0_i
        yield mu_S0_f
      # return zip(np.array(mu_S0),np.array(mu_Mo0), var_S, var_Mo)
    #mu_S0,mu_Mo0,x_S,x_Mo=chemical_potential(E_MoS2,mu_Mobcc,mu_S8,var_S,var_Mo)


    # In[7]:


    def Form_S(E,E0,sign,n,mat=False):
       # chemical_potential(E_MoS2,mu_Mobcc,mu_S8)
        Ef=[]
        mu_S0,mu_Mo0,x_S,x_Mo,mu_S0_i,mu_S0_f=chemical_pot_1(E_MoS2,mu_Mobcc,mu_S8)
        if (mat==False):
           for i in mu_S0:
              Ef.append(E-E0+(sign*n*i))
        else:
           for i in mu_Mo0:
              Ef.append(E-E0+(sign*n*i))
        return Ef


    # In[8]:


    def free_energy(omega,Ti,Tf,step):
       F=[]

       omega = np.array(omega)
       for T in arange(Ti,Tf,step):

         temp3=np.array([(hb*i/2 + k*T*np.log(1-math.exp(-(hb*i)/(k*T)))) for i in omega])
         F.append(np.sum(temp3))
       return F


    # In[9]:


    def DeltaF(X,Y):
       deltaF=[]
       F1=free_energy(X,470,1700,500)
       F2=free_energy(Y,470,1700,500)
       zip_object = zip(F1, F2)
       for i, j in zip_object:
           deltaF.append(i-j)
       return deltaF


    # In[10]:


    def moment():
        molc=read('S8.in',format='aims')
        mass=molc.get_masses()
        mass=sum(mass)* amu2g
        Is=molc.get_moments_of_inertia()
        Is=Is*amu2g*ang2cm*ang2cm
        I=np.sqrt(Is[0])*np.sqrt(Is[1])*np.sqrt(Is[2])
        yield Is
        yield I
        yield mass


    # In[11]:


    def pressure(Ti,Tf,tstep):
        pS=[]
        xr=np.array(np.arange(Ti,Tf,tstep))
        xr=1/xr
        yr=(4.1879-(3209*xr))
        pS=np.exp(yr)
        return pS


    # ## $\mu=   \mu_{0}+ kT \ln\frac{p}{p_{0}}+ E_{DFT}+ \sum_{i}\frac{\hbar \omega_{i}}{2}$
    # ### $\mu_{0}=R\ln\frac{Z^{0}}{V} \frac{kT}{p_{0}} = -kT\ln [(\frac{2 \pi m }{ h^{2}})^{\frac{3}{2}}\frac{(kT)^{\frac{5}{2}}}{p_{0}}]- kT \ln(\frac{\pi^{\frac{1}{2}}}{\sigma}) \\ -kT \ln ((\frac{8 \pi kT}{h^{2}})^{\frac{3}{2}}I_{A}^{\frac{1}{2}} I_{B}^{\frac{1}{2}} I_{C}^{\frac{1}{2}})+kT \sum^{3N-6}_{i} \ln(1-exp(-\beta h \omega_{i})) $

    # ### $ A=\ln( \frac{Z^{0}_{trans} kT}{V})=\ln [(\frac{2 \pi m }{ h^{2}})^{\frac{3}{2}}\frac{(kT)^{\frac{5}{2}}}{p_{0}}]= \ln [(2 \pi m )^{\frac{3}{2}}\frac{(kT)^{\frac{5}{2}}}{h^{3} p_{0}}]$

    # ###  $ B= \ln Z^{0}_{rot}=\ln(\frac{\pi^{\frac{1}{2}}}{\sigma})+ \ln ((\frac{8 \pi kT}{h^{2}})^{\frac{3}{2}}I_{A}^{\frac{1}{2}} I_{B}^{\frac{1}{2}} I_{C}^{\frac{1}{2}})$

    # ### $ C=\ln Z^{0}_{vib}= -\sum^{3N-6}_{i} \ln(1-exp(-\frac{ \hbar \omega_{i}}{kT}))$

    # ### $ F=kT \ln\frac{p}{p_{0}}$

    # ### $ \mu_0 =- KT(A+B-C)$

    # ## $\mu=- KT(A+B-C) + E +  E_{DFT}+ \sum_{i}\frac{\hbar \omega_{i}}{2}$

    # In[12]:

    def conc(r,g,Ti,Tf,tstep):
        rho=[]
        T=arange(Ti,Tf,tstep)
        t=k*T
        for i,j in zip(r,t):
            temp=(g)*math.exp(-i/j)
            rho.append(temp)
        return rho

    def chemical_pot_2(ES,Ti,Tf,tstep,sigma,E,E0,sign,n,w,w0,wS8,mat=False):

        Is,I, m= moment()
        pS=pressure(Ti,Tf,tstep)
        Ef=np.zeros((len(pS),len(pS)))
        concent=np.zeros((len(pS),len(pS)))
        A=[]
        B=[]
        C=[]
        D=[]
        F=[]
        mu_0=[]
        mu_S=[]
        pindx=0


        for pindx,p in enumerate(list(pS*p0)):

            Tindx=0
            for Tindx,T in enumerate(list(arange(Ti,Tf,tstep))):

                A=np.log((((2*pi*m)**(3/2))*((kk*T)**(5/2)))/(p0*(h**3)))
                B=np.log(np.sqrt(pi)/sigma)+ np.log((((8*pi*kk*T)/(h**2))**(3/2))*I)
                temp=np.array([(np.log(1-math.exp(-(hbar*i)/(kk*T)))) for i in wS8])
                C=np.sum(temp)
                F.append(k*T*np.log(p/p0))
                mu_0.append(-k*T*(A+B-C))

            temp2=np.array([(hb*i)/(2) for i in wS8])
            D=np.sum(temp2)
            mu_S=(np.array(mu_0) + np.array(F) + D)/8 + ES
            mu_Mo=(E_MoS2-2*mu_S)
            if (mat==False):
               Ef[pindx,:] = np.array([E-E0+(sign*n*a)+b for a, b in zip(mu_S,DeltaF(w,w0))])
               concent[pindx,:] = np.array(conc(Ef[pindx,:],1,Ti,Tf,tstep))
            else:
               Ef[pindx,:] = np.array([E-E0+(sign*n*a)+b for a, b in zip(mu_Mo,DeltaF(w,w0))])
               concent[pindx,:] = np.array(conc(Ef[pindx,:],1,Ti,Tf,tstep))
        yield concent
        yield Ef
       # return conc, Ef


    # In[13]:


    def plotting():
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(x_S,Ef_addS,'b', label='addon S')
        ax1.plot(x_S,Ef_VS,'r', label='Mono S vacancy')
        ax1.plot(x_S,Ef_VS22, 'g', label='di S vacancy(neig)')
        ax1.plot(x_S,Ef_VS2, 'k', label='di S vacancy(up$\&$down)')
        ax1.set_xlabel(r'$\Delta \mu_{S}$', fontsize=12)
        ax1.set_ylabel('Formation energies eV', fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        ax2 = ax1.twiny()
        ax2.plot(x_Mo,Ef_Mo, 'magenta', label='Mono Mo vacancy')
        ax2.set_xlabel(r'$\frac{1}{2}\Delta \mu_{Mo}$', fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        ax1.legend(loc=4)
        ax2.legend(loc=1)

        plt.savefig('Efvsmu.png',  bbox_inches="tight",dpi=400)
        plt.savefig('Efvsmu.pdf')
        plt.show()


    # In[14]:


    def plotting_2(Ti,Tf,tstep):
        fig = plt.figure(figsize=plt.figaspect(0.5))
        #pS=pressure(Ti,Tf,tstep)
     #===============
    #  First subplot
    #===============
    # set up the axes for the first plot
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        T=arange(Ti,Tf,tstep)
        X, Y = np.meshgrid(pS, T)
        ax.plot_surface(Y, X, Ef_addS_low.T, rstride=1, cstride=1,color='red',shade=False)
        ax.plot_surface(Y, X, Ef_VS_low.T, rstride=1, cstride=1,color='blue',shade=False)
        ax.plot_surface(Y, X, Ef_VS2_low.T, rstride=1, cstride=1,color='green',shade=False)
        ax.plot_surface(Y, X, Ef_VS22_low.T, rstride=1, cstride=1,color='k',shade=False)
        ax.plot_surface(Y, X, Ef_VMo_low.T, rstride=1, cstride=1,color='orange',shade=False)
        ax.set_xlabel(' Tempreture [k]')
        ax.set_zlabel(' Formation Energy [eV]')
        ax.set_ylabel(' Pressure [atm]')
        ax.set_title('MoS2 low S environment')
    #===============
    # Second subplot
    #===============
    # set up the axes for the second plot
        ax = fig.add_subplot(1, 2, 2, projection='3d')

    # plot a 3D wireframe like in the example mplot3d/wire3d_demo

        ax.plot_surface(Y, X, Ef_addS_high.T, rstride=1, cstride=1,color='red',shade=False)
        ax.plot_surface(Y, X, Ef_VS_high.T, rstride=1, cstride=1,color='blue',shade=False)
        ax.plot_surface(Y, X, Ef_VS2_high.T, rstride=1, cstride=1,color='green',shade=False)
        ax.plot_surface(Y, X, Ef_VS22_high.T, rstride=1, cstride=1,color='k',shade=False)
        ax.plot_surface(Y, X, Ef_VMo_high.T, rstride=1, cstride=1,color='orange',shade=False)
        ax.set_xlabel(' Tempreture [k]')
        ax.set_zlabel(' Formation Energy [eV]')
        ax.set_ylabel(' Pressure [atm]')
        ax.set_title('MoS2 high S environment')

        b1 = plt.Rectangle((0, 0), 1, 1, fc="red")
        b2 = plt.Rectangle((0, 0), 1, 1, fc="blue")
        b3 = plt.Rectangle((0, 0), 1, 1, fc="green")
        b4 = plt.Rectangle((0, 0), 1, 1, fc="k")
        b5 = plt.Rectangle((0, 0), 1, 1, fc="orange")

        ax.legend([b1, b2,b3,b4,b5], ['adatom', 'Mono S', 'di up$\&$down','di neigh','Mono Mo'])
        #plt.savefig('Efvsall.png',  bbox_inches="tight",dpi=400)
        #plt.savefig('Efvsall.pdf')
        plt.show()

    def plotting_3(Ti,Tf,tstep):
        fig = plt.figure(figsize=plt.figaspect(0.5))
        #pS=pressure(Ti,Tf,tstep)
     #===============
    #  First subplot
    #===============
    # set up the axes for the first plot
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        T=arange(Ti,Tf,tstep)
        X, Y = np.meshgrid(pS, T)
        ax.plot_surface(Y, X, conc_addS_low.T, rstride=1, cstride=1,color='red',shade=False)
        ax.plot_surface(Y, X, conc_VS_low.T, rstride=1, cstride=1,color='blue',shade=False)
        ax.plot_surface(Y, X, conc_VS2_low.T, rstride=1, cstride=1,color='green',shade=False)
        ax.plot_surface(Y, X, conc_VS22_low.T, rstride=1, cstride=1,color='k',shade=False)
        ax.plot_surface(Y, X, conc_VMo_low.T, rstride=1, cstride=1,color='orange',shade=False)
        ax.set_xlabel(' Tempreture [k]')
        ax.set_zlabel(' COncentration')
        ax.set_ylabel(' Pressure [atm]')
        ax.set_title('MoS2 low S environment')
    #===============
    # Second subplot
    #===============
    # set up the axes for the second plot
        ax = fig.add_subplot(1, 2, 2, projection='3d')

    # plot a 3D wireframe like in the example mplot3d/wire3d_demo

        ax.plot_surface(Y, X, conc_addS_high.T, rstride=1, cstride=1,color='red',shade=False)
        ax.plot_surface(Y, X, conc_VS_high.T, rstride=1, cstride=1,color='blue',shade=False)
        ax.plot_surface(Y, X, conc_VS2_high.T, rstride=1, cstride=1,color='green',shade=False)
        ax.plot_surface(Y, X, conc_VS22_high.T, rstride=1, cstride=1,color='k',shade=False)
        ax.plot_surface(Y, X, conc_VMo_high.T, rstride=1, cstride=1,color='orange',shade=False)
        ax.set_xlabel(' Tempreture [k]')
        ax.set_zlabel(' concentration')
        ax.set_ylabel(' Pressure [atm]')
        ax.set_title('MoS2 high S environment')

        b1 = plt.Rectangle((0, 0), 1, 1, fc="red")
        b2 = plt.Rectangle((0, 0), 1, 1, fc="blue")
        b3 = plt.Rectangle((0, 0), 1, 1, fc="green")
        b4 = plt.Rectangle((0, 0), 1, 1, fc="k")
        b5 = plt.Rectangle((0, 0), 1, 1, fc="orange")

        ax.legend([b1, b2,b3,b4,b5], ['adatom', 'Mono S', 'di up$\&$down','di neigh','Mono Mo'])
        #plt.savefig('Efvsall.png',  bbox_inches="tight",dpi=400)
        #plt.savefig('Efvsall.pdf')
        plt.show()

# In[15]:


    """main routine"""
    # Reading Energies
    E0=read_energy('pristine')
    E1=read_energy('addX')
    E2=read_energy('VX')
    E3=read_energy('VX2')
    E4=read_energy('VX22')
    E5=read_energy('VM')
    E_MoS2=read_energy('primitive')
    ES8=read_energy('S8')
    EMo=read_energy('Mo')
    # Chemical Potential
    mu_MoS2=E0/75
    mu_S8=ES8/8
    mu_Mobcc=EMo/2
    mu_S0,mu_Mo0,x_S,x_Mo,mu_S0_i,mu_S0_f=chemical_pot_1(E_MoS2,mu_Mobcc,mu_S8)
    # Calculation Formation Energies
    Ef_addS=Form_S(E1,E0,-1,1,False)
    Ef_VS=Form_S(E2,E0,+1,1,False)
    Ef_VS2=Form_S(E3,E0,+1,2,False)
    Ef_VS22=Form_S(E4,E0,+1,2,False)
    Ef_Mo=Form_S(E5,E0,+1,1,True)
    # Plotting Formation energy versus chemical potential
    plotting()

    # Reading the frequencies
    wW=read_omega('wM')
    wS8=read_omega('wX')
    w0=read_omega('w0')
    w1=read_omega('waddX')
    w2=read_omega('wVX')
    w3=read_omega('wVX2')
    w4=read_omega('wVX22')
    w5=read_omega('wVM')
    # Calculating Chemical potential as function of temperature and pressure

    pS=pressure(470,1700,500)

    conc_addS_low,Ef_addS_low=chemical_pot_2(mu_S0_i,470,1700,500,8,E1,E0,-1,1,w1,w0,wS8,mat=False)
    conc_addS_high,Ef_addS_high=chemical_pot_2(mu_S0_f,470,1700,500,8,E1,E0,-1,1,w1,w0,wS8,mat=False)
    conc_VS_low,Ef_VS_low=chemical_pot_2(mu_S0_i,470,1700,500,8,E2,E0,1,1,w2,w0,wS8,mat=False)
    conc_VS_high,Ef_VS_high=chemical_pot_2(mu_S0_f,470,1700,500,8,E2,E0,1,1,w2,w0,wS8,mat=False)
    conc_VS2_low,Ef_VS2_low=chemical_pot_2(mu_S0_i,470,1700,500,8,E3,E0,1,2,w3,w0,wS8,mat=False)
    conc_VS2_high,Ef_VS2_high=chemical_pot_2(mu_S0_f,470,1700,500,8,E3,E0,1,2,w3,w0,wS8,mat=False)
    conc_VS22_low,Ef_VS22_low=chemical_pot_2(mu_S0_i,470,1700,500,8,E4,E0,1,2,w4,w0,wS8,mat=False)
    conc_VS22_high,Ef_VS22_high=chemical_pot_2(mu_S0_f,470,1700,500,8,E4,E0,1,2,w4,w0,wS8,mat=False)
    conc_VMo_low,Ef_VMo_low=chemical_pot_2(mu_S0_f,470,1700,500,8,E5,E0,1,1,w5,w0,wS8,mat=True)
    conc_VMo_high,Ef_VMo_high=chemical_pot_2(mu_S0_i,470,1700,500,8,E5,E0,1,1,w5,w0,wS8,mat=True)

    # Plotting
    plotting_2(470,1700,500)


    # Plotting
    plotting_3(470,1700,500)

# In[16]:


if __name__ == "__main__":
    main()

