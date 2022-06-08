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
import sys


# In[2]:



### conversions and constants ###
Pi=np.pi
cmhz=29979245800.0*2*Pi #cm^-1 to Hz
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

#__all__ = ['main']
# Spliting function
def split_line(lines):
    """Split input line"""
    line_array = np.array(lines.strip().split(" "))
    line_vals = line_array[line_array != ""]
    return line_vals


# In[4]:

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
def indices(the_list, val):
     """Always returns a list containing the indices of val in the_list"""
     retval = []
     last = 0
     while val in the_list[last:]:
          i = the_list[last:].index(val)
          retval.append(last + i)
          last += i + 1
     return retval

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


def main(filename,Ti,Tf,tstep,pi,pf,pstep):
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
         if T==0:
            temp3=np.array([(hb*i/2) for i in omega])
            F.append(np.sum(temp3))
         else:
            temp3=np.array([(hb*i/2 + k*T*np.log(1-math.exp(-(hb*i)/(k*T)))) for i in omega])
            F.append(np.sum(temp3))
       return F


    # In[9]:


    def DeltaF(X,Y,Ti,Tf,step):
       deltaF=[]
       F1=free_energy(X,Ti,Tf,step)
       F2=free_energy(Y,Ti,Tf,step)
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


    def pressure(pi,pf,pstep):
        pS=np.arange(pi,pf,pstep)
       # pS=[]
       # xr=np.array(np.arange(Ti,Tf,tstep))
       # xr=1/xr
       # yr=(4.1879-(3209*xr))
       # pS=np.exp(yr)
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

   # def conc(r,g,Ti,Tf,tstep):
   #     rho=[]
   #     T=arange(Ti,Tf,tstep)
   #     t=k*T
   #     for i,j in zip(r,t):
   #         temp=(g)*math.exp(-i/j)
   #         rho.append(temp)
   #     return rho
    def translational(Ti,Tf,tstep):
        A=[]
        trans=[]
        Is,I, m= moment()
        Tindx=0
        for Tindx,T in enumerate(list(arange(Ti,Tf,tstep))):

            if T==0:
                A=0
                trans.append(-k*T*A)
            else:
                A=np.log((((2*Pi*m)**(3/2))*((kk*T)**(5/2)))/(p0*(h**3)))
                trans.append(-k*T*A)
       # print('trans', np.shape(trans))
        return np.array(trans)
    def rotational(Ti,Tf,tstep,sigma):
        B=[]
        rot=[]
        Is,I, m= moment()
        Tindx=0
        for Tindx,T in enumerate(list(arange(Ti,Tf,tstep))):
            if T==0:
              B=0
              rot.append(-k*T*B)
            else:
              B=np.log(np.sqrt(Pi)/sigma)+ np.log((((8*Pi*kk*T)/(h**2))**(3/2))*I)
              rot.append(-k*T*B)
       # print('rot', np.shape(rot))
        return np.array(rot)
    def vibrational(Ti,Tf,tstep,elec):
        C=[]
        D=[]
        vib=[]
        Is,I, m= moment()
        Tindx=0
        for Tindx,T in enumerate(list(arange(Ti,Tf,tstep))):
            temp2=np.array([(hb*i)/(2) for i in wS8])
            D=np.sum(temp2)
            if T==0:
               vib.append(D+8*elec)
            else:
               temp=np.array([(np.log(1-math.exp(-(hbar*i)/(kk*T)))) for i in wS8])
               C=np.sum(temp)
               vib.append(k*T*C+D+8*elec)
         #   rot=rot+D
       # print('vib', np.shape(vib))
        return vib
    def mu_0(Ti,Tf,tstep,elec,sigma):
        mu_0=[]
        rot=rotational(Ti,Tf,tstep,sigma)
        vib=vibrational(Ti,Tf,tstep,elec)
        trans=translational(Ti,Tf,tstep)
        mu_0=rot+vib+trans
#        print('mu_0', np.shape(mu_0))
        return mu_0/8
    def press(pi,pf,pstep,Ti,Tf,tstep):
        pindx=0
        pS=pressure(pi,pf,pstep)
        T=np.arange(Ti,Tf,tstep)
        F=np.zeros((len(pS),len(T)))
        for pindx,p in enumerate(list(pS*p0)):
            if p==0:
                F[pindx,:]=np.zeros(len(T))
            else:
                if Ti==0:
                   F[pindx,0]=0
                   F[pindx,1:]=np.array([k*T*np.log(p/p0) for T in arange(Ti,Tf,tstep) if T!=0])
                else:
                   F[pindx,:]=np.array([k*T*np.log(p/p0) for T in arange(Ti,Tf,tstep)])
    #    print('F', F)
    #    print('F/8', F/8)
        return F/8
    def mu(Ti,Tf,tstep,elec,sigma,pi,pf,pstep):
        pindx=0
        pS=pressure(pi,pf,pstep)
        T=np.arange(Ti,Tf,tstep)
    #    mu_00=np.zeros((len(pS),len(T)))
        mu_10=np.zeros(len(T))
        mu_10=mu_0(Ti,Tf,tstep,elec,sigma)
#        print('mu_10', np.shape(mu_10))
        p_dep=np.zeros((len(pS),len(T)))
    #    mu_00=np.column_stack((mu_10,mu_10))
#        print('mu_00', np.shape(mu_00))
        p_dep=press(pi,pf,pstep,Ti,Tf,tstep)
        mu_S=np.zeros((len(pS),len(T)))
        mu_Mo=np.zeros((len(pS),len(T)))
        for pindx,p in enumerate(list(pS*p0)):
           mu_S[pindx,:]=np.array(mu_10+p_dep[pindx,:])
           mu_Mo[pindx,:]=(E_MoS2-2*mu_S[pindx,:])
        yield mu_S
        yield mu_Mo

    def formation_energy(ES,pi,pf,pstep,Ti,Tf,tstep,sigma,E,E0,sign,n,w,w0,wS8,mat=False):

        pS=pressure(pi,pf,pstep)
        T=np.arange(Ti,Tf,tstep)
        Ef=np.zeros((len(pS),len(T)))
        concent=np.zeros((len(pS),len(T)))
        mu_S,mu_Mo=mu(Ti,Tf,tstep,ES,sigma,pi,pf,pstep)
        pindx=0
        for pindx,p in enumerate(list(pS*p0)):
          # print(pindx)
           if (mat==False):
               Ef[pindx,:] = np.array([E-E0+(sign*n*a)+b for a, b in zip(mu_S[pindx,:] ,DeltaF(w,w0,Ti,Tf,tstep))])
              # concent[pindx,:] = np.array(conc(Ef[pindx,:],1,Ti,Tf,tstep))
           else:
               Ef[pindx,:] = np.array([E-E0+(sign*n*a)+b for a, b in zip(mu_Mo[pindx,:] ,DeltaF(w,w0,Ti,Tf,tstep))])
             #  concent[pindx,:] = np.array(conc(Ef[pindx,:],1,Ti,Tf,tstep))
    #    yield concent
        return Ef

    def plotting():
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(x_S,Ef_addS,'b', label='adatom Se')
        ax1.plot(x_S,Ef_VS,'r', label='Mono Se vacancy')
        ax1.plot(x_S,Ef_VS22, 'g', label='di Se vacancy(neig)')
        ax1.plot(x_S,Ef_VS2, 'k', label='di Se vacancy(up$\&$down)')
        ax1.set_xlabel(r'$\Delta \mu_{Se}$', fontsize=12)
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

       # plt.savefig('Efvsmu.png',  bbox_inches="tight",dpi=400)
        plt.savefig('Ef_mu_'+str(filename)+'.pdf')
        plt.show()


    # In[14]:


    def plotting_2(Ti,Tf,tstep):
        fig = plt.figure(figsize=plt.figaspect(0.5))
        #pS=pressure(Ti,Tf,tstep)
     #===============
    #  First subplot
    #===============
    # set up the axes for the first plot
        ax = fig.add_subplot(1, 3, 1, projection='3d')
        T=arange(Ti,Tf,tstep)
        X, Y = np.meshgrid(pS, T)
        ax.plot_surface(Y, X, Ef_addS_low.T, rstride=1, cstride=1,color='red',shade=False)
        ax.plot_surface(Y, X, Ef_VS_low.T, rstride=1, cstride=1,color='royalblue',shade=False)
        ax.plot_surface(Y, X, Ef_VS2_low.T, rstride=1, cstride=1,color='darkgreen',shade=False)
        ax.plot_surface(Y, X, Ef_VS22_low.T, rstride=1, cstride=1,color='darkmagenta',shade=False)
        ax.plot_surface(Y, X, Ef_VMo_low.T, rstride=1, cstride=1,color='darkorange',shade=False)
        ax.set_xlabel(' Tempreture [k]')
        ax.set_zlabel(' Formation Energy [eV]')
        ax.set_ylabel(' Pressure [atm]')
        ax.set_title(str(filename)+' low S environment')
    #===============
        ax = fig.add_subplot(1, 3, 2, projection='3d')

    # plot a 3D wireframe like in the example mplot3d/wire3d_demo
      #  ax.legend(bbox_to_anchor=(0.7, 0.4))
        ax.plot_surface(Y, X, Ef_addS_cross.T, rstride=1, cstride=1,color='red',shade=False)
        ax.plot_surface(Y, X, Ef_VS_cross.T, rstride=1, cstride=1,color='royalblue',shade=False)
        ax.plot_surface(Y, X, Ef_VS2_cross.T, rstride=1, cstride=1,color='darkgreen',shade=False)
        ax.plot_surface(Y, X, Ef_VS22_cross.T, rstride=1, cstride=1,color='darkmagenta',shade=False)
        ax.plot_surface(Y, X, Ef_VMo_cross.T, rstride=1, cstride=1,color='darkorange',shade=False)
        ax.set_xlabel(' Tempreture [k]')
        ax.set_zlabel(' Formation Energy [eV]')
        ax.set_ylabel(' Pressure [atm]')
        ax.set_title(str(filename)+' cross environment')
    # Second subplot
    # Second subplot
    #===============
    # set up the axes for the second plot
        ax = fig.add_subplot(1, 3, 3, projection='3d')

    # plot a 3D wireframe like in the example mplot3d/wire3d_demo

        ax.plot_surface(Y, X, Ef_addS_high.T, rstride=1, cstride=1,color='red',shade=False)
        ax.plot_surface(Y, X, Ef_VS_high.T, rstride=1, cstride=1,color='royalblue',shade=False)
        ax.plot_surface(Y, X, Ef_VS2_high.T, rstride=1, cstride=1,color='darkgreen',shade=False)
        ax.plot_surface(Y, X, Ef_VS22_high.T, rstride=1, cstride=1,color='darkmagenta',shade=False)
        ax.plot_surface(Y, X, Ef_VMo_high.T, rstride=1, cstride=1,color='darkorange',shade=False)
        ax.set_xlabel(' Tempreture [k]')
        ax.set_zlabel(' Formation Energy [eV]')
        ax.set_ylabel(' Pressure [atm]')
        ax.set_title(str(filename)+' high S environment')

        b1 = plt.Rectangle((0, 0), 1, 1, fc="red")
        b2 = plt.Rectangle((0, 0), 1, 1, fc="blue")
        b3 = plt.Rectangle((0, 0), 1, 1, fc="green")
        b4 = plt.Rectangle((0, 0), 1, 1, fc="k")
        b5 = plt.Rectangle((0, 0), 1, 1, fc="orange")

       # ax.legend([b1, b2,b3,b4,b5], ['adatom', 'Mono S', 'di up$\&$down','di neigh','Mono Mo'])
        plt.savefig('Efvsall.png',  bbox_inches="tight",dpi=400)
       # plt.savefig('Efvsall.pdf')
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
      #  plt.savefig('rhovsall.png',  bbox_inches="tight",dpi=400)
        plt.savefig('rhovsall.pdf')
        plt.show()

# In[15]:
    def plotting_4(Ti,Tf,tstep):
        fig = plt.figure(figsize=plt.figaspect(0.5))
 #       #pS=pressure(Ti,Tf,tstep)
 #    #===============
 #   #  First subplot
 #   #===============
 #   # set up the axes for the first plot
        ax = fig.add_subplot(1, 3, 1)
        T=arange(Ti,Tf,tstep)
        pindx=1
        ax.plot(T,Ef_addS_low[pindx,:].T,'b', label='adatom Se')
        ax.plot(T,Ef_VS_low[pindx,:].T,'r', label='Mono Se vacancy')
        ax.plot(T,Ef_VS22_low[pindx,:].T, 'g', label='di Se vacancy(neig)')
        ax.plot(T,Ef_VS2_low[pindx,:].T, 'k', label='di Se vacancy(up$\&$down)')
        ax.plot(T,Ef_VMo_low[pindx,:].T, 'magenta', label='Mono Mo vacancy')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        #ax.legend(loc=4)
        ax.set_xlabel(' Tempreture [k]')
        ax.set_ylabel(' Formation Energy [eV]')
        ax.set_title(str(filename)+' low S environment')
    #===============
        ax = fig.add_subplot(1, 3, 2)

    # plot a 3D wireframe like in the example mplot3d/wire3d_demo
        ax.plot(T,Ef_addS_cross[pindx,:].T,'b', label='adatom Se')
        ax.plot(T,Ef_VS_cross[pindx,:].T,'r', label='Mono Se vacancy')
        ax.plot(T,Ef_VS22_cross[pindx,:].T, 'g', label='di Se vacancy(neig)')
        ax.plot(T,Ef_VS2_cross[pindx,:].T, 'k', label='di Se vacancy(up$\&$down)')
        ax.plot(T,Ef_VMo_cross[pindx,:].T, 'magenta', label='Mono Mo vacancy')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
      #  ax.legend(bbox_to_anchor=(0.7, 0.4))
        ax.set_xlabel(' Tempreture [k]')
        ax.set_ylabel(' Formation Energy [eV]')
        ax.set_title(str(filename)+' cross environment')
    # Second subplot
    #===============
    # set up the axes for the second plot
        ax = fig.add_subplot(1, 3, 3)

    # plot a 3D wireframe like in the example mplot3d/wire3d_demo
        ax.plot(T,Ef_addS_high[pindx,:].T,'b', label='adatom Se')
        ax.plot(T,Ef_VS_high[pindx,:].T,'r', label='Mono Se vacancy')
        ax.plot(T,Ef_VS22_high[pindx,:].T, 'g', label='di Se vacancy(neig)')
        ax.plot(T,Ef_VS2_high[pindx,:].T, 'k', label='di Se vacancy(up$\&$down)')
        ax.plot(T,Ef_VMo_high[pindx,:].T, 'magenta', label='Mono Mo vacancy')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
      #  ax.legend(bbox_to_anchor=(0.7, 0.4))
        ax.set_xlabel(' Tempreture [k]')
        ax.set_ylabel(' Formation Energy [eV]')
        ax.set_title(str(filename)+' high S environment')

      #  b1 = plt.Rectangle((0, 0), 1, 1, fc="red")
      #  b2 = plt.Rectangle((0, 0), 1, 1, fc="blue")
      #  b3 = plt.Rectangle((0, 0), 1, 1, fc="green")
      #  b4 = plt.Rectangle((0, 0), 1, 1, fc="k")
      #  b5 = plt.Rectangle((0, 0), 1, 1, fc="orange")

      #  ax.legend([b1, b2,b3,b4,b5], ['adatom', 'Mono S', 'di up$\&$down','di neigh','Mono Mo'])
       # plt.savefig('Efvsall.png',  bbox_inches="tight",dpi=400)
        plt.savefig('Ef_T_'+str(filename)+'.pdf')
        plt.show()



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

    DEf=[a-b for a, b in zip(Ef_addS, Ef_VS)]
    valu=find_nearest(DEf, 0)
    indx=indices(DEf,valu)
    print('Ef_cross',Ef_addS[indx[0]])
    print('u_cross',x_S[indx[0]])
    print('u_cross',x_S[indx[0]]+mu_S0_f)
    G_cross=Ef_addS[indx[0]]
    mu_cross=x_S[indx[0]]+mu_S0_f
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
    T=[float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4])]
    pS=pressure(float(sys.argv[5]),float(sys.argv[6]),float(sys.argv[7]))
    p=[float(sys.argv[5]),float(sys.argv[6]),float(sys.argv[7])]
    Ef_addS_low=formation_energy(mu_S0_i,p[0],p[1],p[2],T[0],T[1],T[2],8,E1,E0,-1,1,w1,w0,wS8,mat=False)
    Ef_addS_high=formation_energy(mu_S0_f,p[0],p[1],p[2],T[0],T[1],T[2],8,E1,E0,-1,1,w1,w0,wS8,mat=False)
    Ef_VS_low=formation_energy(mu_S0_i,p[0],p[1],p[2],T[0],T[1],T[2],8,E2,E0,1,1,w2,w0,wS8,mat=False)
    Ef_VS_high=formation_energy(mu_S0_f,p[0],p[1],p[2],T[0],T[1],T[2],8,E2,E0,1,1,w2,w0,wS8,mat=False)
    Ef_VS2_low=formation_energy(mu_S0_i,p[0],p[1],p[2],T[0],T[1],T[2],8,E3,E0,1,2,w3,w0,wS8,mat=False)
    Ef_VS2_high=formation_energy(mu_S0_f,p[0],p[1],p[2],T[0],T[1],T[2],8,E3,E0,1,2,w3,w0,wS8,mat=False)
    Ef_VS22_low=formation_energy(mu_S0_i,p[0],p[1],p[2],T[0],T[1],T[2],8,E4,E0,1,2,w4,w0,wS8,mat=False)
    Ef_VS22_high=formation_energy(mu_S0_f,p[0],p[1],p[2],T[0],T[1],T[2],8,E4,E0,1,2,w4,w0,wS8,mat=False)
    Ef_VMo_low=formation_energy(mu_S0_f,p[0],p[1],p[2],T[0],T[1],T[2],8,E5,E0,1,1,w5,w0,wS8,mat=True)
    Ef_VMo_high=formation_energy(mu_S0_i,p[0],p[1],p[2],T[0],T[1],T[2],8,E5,E0,1,1,w5,w0,wS8,mat=True)
    Ef_addS_cross=formation_energy(mu_cross,p[0],p[1],p[2],T[0],T[1],T[2],8,E1,E0,-1,1,w1,w0,wS8,mat=False)
    Ef_VS_cross=formation_energy(mu_cross,p[0],p[1],p[2],T[0],T[1],T[2],8,E2,E0,1,1,w2,w0,wS8,mat=False)
    Ef_VS2_cross=formation_energy(mu_cross,p[0],p[1],p[2],T[0],T[1],T[2],8,E3,E0,1,2,w3,w0,wS8,mat=False)
    Ef_VS22_cross=formation_energy(mu_cross,p[0],p[1],p[2],T[0],T[1],T[2],8,E4,E0,1,2,w4,w0,wS8,mat=False)
    Ef_VMo_cross=formation_energy(mu_cross,p[0],p[1],p[2],T[0],T[1],T[2],8,E5,E0,1,1,w5,w0,wS8,mat=True)

    # Plotting
    plotting_2(float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]))


    # Plotting
   # plotting_3(float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]))

   # plotting_4(float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]))
# In[16]:
if __name__ == "__main__":
      main(str(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6]),float(sys.argv[7]))
