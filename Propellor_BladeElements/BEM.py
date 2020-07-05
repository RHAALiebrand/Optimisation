import numpy as np
import scipy.interpolate as inter
import subprocess as sp
import os
import matplotlib.pyplot as plt
# Import and reshape geometry file such that N can be adapted
propfile='prop378.dat'


def geometry(propfile,r,R):
    propgeo=np.genfromtxt(propfile,skip_header=5)
    chord=inter.interp1d( propgeo[:,0]*R , propgeo[:,1]*2*R)(r)
    pitchgeo=inter.interp1d( propgeo[:,0]*R , propgeo[:,2])(r)
    pitch=np.arctan2(pitchgeo,(np.pi*r/R))*180./np.pi
    return chord,pitch

  # Theta is in degrees

def cl_cd(alpha,Re,radius,Mach,iter):
    def Cmd(cmd):
        process.stdin.write(cmd + '\n')

    process = sp.Popen(r'xfoil.exe', stdin=sp.PIPE, stderr=sp.PIPE, stdout=sp.PIPE)
    process.stderr.close()

    Cmd('LOAD')
    Cmd('Clark-Y10.dat')
    Cmd('OPER')
    Cmd('ITER')
    Cmd('100')
    Cmd('Visc')
    Cmd(str(Re))
    #Cmd('mach'+str(Mach))
    Cmd('pacc')
    Cmd('results/'+str(radius)+'_'+str(iter)+'.dat')
    Cmd(' ')
    Cmd('Alfa'+str(alpha))
    Cmd(' ')
    Cmd('QUIT')

    process.stdout.close()
    process.stdin.close()
    process.wait()

    return

step=5.
def xfoil_fail(alpha,Re,radius,Mach,iter):
    def Cmd(cmd):
        process.stdin.write(cmd + '\n')

    process = sp.Popen(r'xfoil.exe', stdin=sp.PIPE, stderr=sp.PIPE, stdout=sp.PIPE)
    process.stderr.close()

    Cmd('LOAD')
    Cmd('Clark-Y10.dat')
    Cmd('OPER')
    Cmd('ITER')
    Cmd('100')
    Cmd('Visc')
    Cmd(str(Re))
    #Cmd('mach'+str(Mach))
    Cmd('pacc')
    Cmd('results/fail_'+str(radius)+'_'+str(iter)+'.dat')
    Cmd(' ')
    Cmd('Aseq')
    if alpha<0.0:
        Cmd(str(0))
        Cmd(str(alpha -  2))
        Cmd(str(1.0))
    else:
        Cmd(str(0))
        Cmd(str(alpha + 2))
        Cmd(str(1.0))
    Cmd(' ')
    Cmd('QUIT')

    process.stdout.close()
    process.stdin.close()
    process.wait()

    return

# Here a part of the code is copied and translated into a function to use in the second part of the task
# Here the BEM is a function
def BEM(Vinf,RPM,Nb,r,N,R,chordfraction):
    chord,theta=geometry(propfile,r,R)
    chord=chord*chordfraction
    # Tolerance for iteration
    toll=1.0e-3
    # Set test conditions
    rho=1.225   # Not exactly correct, see table
    omega=RPM/60.*2.0*np.pi
    mu=0.0000181206
    Vax_array=np.zeros(N-1)
    T=np.zeros(N)
    Q=np.zeros(N)
    for i in range(N-1):
        sigma=Nb*chord[i]/(2*np.pi*r[i])
        # Initial guess
        a=0.1
        b=0.01
        # Now start iterations
        iter=0
        finished=False
        while not finished:
            Vax=Vinf*(1+a) #Axial velocity
            Vrad=omega*r[i]*(1-b) # Radial velocity
            phi=np.arctan2(Vax,Vrad) # Inflow angle
            alpha=theta[i]-phi*180/np.pi # Local angle of attack
            Veff=np.sqrt(Vax**2+Vrad**2)
            Mach=Veff/340.294
            Re=rho*Veff*chord[i]/mu
            cl_cd(alpha, Re, r[i],Mach, iter)
            data=np.genfromtxt('results/'+str(r[i])+'_'+str(iter)+'.dat',skip_header=12)
            if len(data)==0:
                xfoil_fail(alpha, Re, r[i], Mach, iter)
                data=np.genfromtxt('results/fail_'+str(r[i])+'_'+str(iter)+'.dat',skip_header=12)
                cl=inter.interp1d(data[:,0],data[:,1],fill_value='extrapolate')(alpha)
                cd = inter.interp1d(data[:, 0], data[:, 2],fill_value='extrapolate')(alpha)
                os.remove('results/fail_'+str(r[i])+'_'+str(iter)+'.dat')
            else:
                cl=data[1]
                cd=data[2]
                os.remove('results/' + str(r[i]) + '_' + str(iter) + '.dat')

            # Determine thrust and torque
            dT_dr=0.5*rho*Veff**2*Nb*chord[i]*(cl*np.cos(phi)-cd*np.sin(phi))
            dQ_dr=0.5*rho*Veff**2*Nb*chord[i]*r[i]*(cd*np.cos(phi)+cl*np.sin(phi))
            # Now apply torque and moment balances as described in the report
            # Assume a inside the brackets to be old a, and solve for the a outside (a_iter)
            # print 'abefore',a
            a_iter=dT_dr/(4.0*np.pi*r[i]*rho*Vinf**2*(1+a))
            # print 'aiter',a_iter
            b_iter=dQ_dr/(4.0*np.pi*r[i]**3*rho*Vinf*(1+a)*omega)
            # Now take the avarge to stabalize iteration
            anew=0.5*(a+a_iter)
            bnew=0.5*(b+b_iter)
            if abs(anew-a)<toll and abs(bnew-b)<toll:
                finished=True
            a=anew
            # print 'a', a
            b=bnew
            iter=iter+1
            if iter>500:
                finished=True
            # Last update + Tip correction
        # Update when iteration is completed and correct for tip loss
        # Hub loss correction
        h = (Nb * (r[i] - 0.2*R)) / (2.0 * r[i] * np.sin(phi))
        # print 'h', np.exp(-h)
        H = 2. / np.pi * np.arccos(np.exp(-h))
        # print 'H',H
        # Tip loss correction
        f=(Nb*(R-r[i]))/(2.0*r[i]*np.sin(phi))
        F=2./np.pi*np.arccos(np.exp(-f))
        dT_dr = 0.5 * rho * Veff ** 2 * Nb * chord[i] * (cl * np.cos(phi) - cd * np.sin(phi))*F*H
        dQ_dr = 0.5 * rho * Veff ** 2 * Nb * chord[i] * r[i] * (cd * np.cos(phi) + cl * np.sin(phi))*F*H
        T[i]=dT_dr
        Q[i]=dQ_dr
        Vax_array[i]=Vax

    thrust=np.sum(T*(0.8*R/(N-1)))
    torque=np.sum(Q*(0.8*R/(N-1)))

    CT=thrust/(rho*(2.*R)**4*(RPM/60.)**2)
    CQ=torque/(rho*(2.*R)**5*(RPM/60.)**2)
    eta=Vinf/((RPM/60.)*(2*R))/(2.0*np.pi)*(CT/CQ)
    return CT,eta,CQ*2*np.pi

# Adapt rotor geometry

# Different number of blades
# Set reference
R_ref=9.5*0.3048*0.5
N=8
r_ref=np.linspace(0.2*R_ref,R_ref,N)
RPM=1400.
Vinf=25.
chordfraction=1.0
Nb_ref=2.0
CT_ref,eta_ref,CP_ref=BEM(Vinf,RPM,Nb_ref,r_ref,N,R_ref,chordfraction)

# Start adapting rotor to 3 and 4 blades, 2 baldes is reference
Nb=np.arange(3.0,4.5,1.0)
CT_array=np.zeros(len(Nb)+1)
eta_array=np.zeros(len(Nb)+1)
CT_array[0]=CT_ref
eta_array[0]=eta_ref
# Start loop over number of blades
for i in range(len(Nb)):
    finished = False
    # Iterative process
    while not finished:
        CT,eta,CP=BEM(Vinf,RPM,Nb[i],r_ref,N,R_ref,chordfraction)
        print CP
        if CP<CP_ref:
            finished=True
        Vinf=Vinf+1.0

    CT_array[i+1]=CT
    eta_array[i+1]=eta
# Plotting results
Nb=np.arange(2.0,4.5,1.0)
plt.figure()
plt.plot(Nb,CT_array)
plt.xlabel('$N_b$ [-]')
plt.ylabel('$C_T$ [-]')

plt.figure()
plt.plot(Nb,eta_array)
plt.xlabel('$N_b$ [-]')
plt.ylabel('$\eta$ [-]')
plt.show()

# Change radius
R=9.5*0.3048*0.5
# Adapt radius (0.8 is used as reference)
R_array=np.array([R,1.2*R])
CT_array=np.zeros(len(R_array))
eta_array=np.zeros(len(R_array))
for i in range(len(R_array)):
    finished=False
    while not finished:
        CT, eta, CP = BEM(Vinf, RPM, 3.0, np.linspace(0.2*R_array[i],R_array[i],N), N, R_array[i])
        print CP
        if CP<CP_ref:
            finished=True
        Vinf=Vinf+1.

        CT_array[i]=CT
        eta_array[i]=eta
# Plot results
R_plot=np.array([0.8,1.0,1.2])*100
CT_array=np.array([CT_ref,CT_array[0],CT_array[1]])
eta_array=np.array([eta_ref,eta_array[0],eta_array[1]])
plt.figure()
plt.plot(R_plot,CT_array)
plt.xlabel('$\%R_{ref}$ [-]')
plt.ylabel('$C_T$ [-]')

plt.figure()
plt.plot(R_plot,eta_array)
plt.xlabel('$\%R_{ref}$ [-]')
plt.ylabel('$\eta$ [-]')
plt.show()

# Change aspect ratio
# Setting reference to 0.9
chord_fraction=0.9
CT_ref,eta_ref,CP_ref=BEM(Vinf,RPM,3.0,r_ref,N,R_ref,chord_fraction)
chord_fraction_array=np.array([1.0,1.1])
CT_array=np.zeros(len(chord_fraction_array))
eta_array=np.zeros(len(chord_fraction_array))
for i in range(len(chord_fraction_array)):
    finished=False
    while not finished:
        CT, eta, CP = BEM(Vinf, RPM, 3.0, r_ref, N, R_ref,chord_fraction_array[i])
        print CP
        if CP<CP_ref:
            finished=True
        Vinf=Vinf+2.0
        CT_array[i]=CT
        eta_array[i]=eta
# Plot results
chord_plot=np.array([chord_fraction,chord_fraction_array[0],chord_fraction_array[1]])*100
CT_array=np.array([CT_ref,CT_array[0],CT_array[1]])
eta_array=np.array([eta_ref,eta_array[0],eta_array[1]])
plt.figure()
plt.plot(chord_plot,CT_array)
plt.xlabel('$\%chord_{ref}$ [-]')
plt.ylabel('$C_T$ [-]')

plt.figure()
plt.plot(chord_plot,eta_array)
plt.xlabel('$\%chord_{ref}$ [-]')
plt.ylabel('$\eta$ [-]')
plt.show()