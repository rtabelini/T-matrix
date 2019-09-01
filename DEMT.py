import numpy as np
shear = []
bulk=[]
dens=[]
macroporos=[]
macroaspect=[]
microporos=[]
microaspect=[]
vvp=[]
vvs=[]
a10_plug_poros = []
a10_plug_perm = []
#    arqshear = open("shear.txt",'r')
#    tshear = arqshear.readlines()
#    for i in range(len(tshear)):
#        shear.append(float(tshear[i].strip()))
    
arqa10_plug_poros = open("a10_plug_poros.txt",'r')
ta10_plug_poros = arqa10_plug_poros.readlines()
for i in range(len(ta10_plug_poros)):
    a10_plug_poros.append(float(ta10_plug_poros[i].strip()))
    
arqa10_plug_perm = open("a10_plug_perm.txt",'r')
ta10_plug_perm = arqa10_plug_perm.readlines()
for i in range(len(ta10_plug_perm)):
    a10_plug_perm.append(float(ta10_plug_perm[i].strip()))

arqbulk = open("km_irineu.txt",'r')
tbulk = arqbulk.readlines()
for i in range(len(tbulk)):
    bulk.append(float(tbulk[i].strip()))

arqshear = open("gm_irineu.txt",'r')
tshear = arqshear.readlines()
for i in range(len(tshear)):
    shear.append(float(tshear[i].strip()))

arqdens = open("dens_irineu.txt",'r')
tdens = arqdens.readlines()
for i in range(len(tdens)):
    dens.append(float(tdens[i].strip()))
    
arqmacroaspect = open("macroaspect_irineu.txt",'r')
tmacroaspect = arqmacroaspect.readlines()
for i in range(len(tmacroaspect)):
    macroaspect.append(float(tmacroaspect[i].strip()))

arqmicroaspect = open("microaspect_irineu.txt",'r')
tmicroaspect = arqmicroaspect.readlines()
for i in range(len(tmicroaspect)):
    microaspect.append(float(tmicroaspect[i].strip()))
    
arqmacroporos = open("macroporos_irineu.txt",'r')
tmacroporos = arqmacroporos.readlines()
for i in range(len(tmacroporos)):
    macroporos.append(float(tmacroporos[i].strip()))

arqmicroporos = open("microporos_irineu.txt",'r')
tmicroporos = arqmicroporos.readlines()
for i in range(len(tmicroporos)):
    microporos.append(float(tmicroporos[i].strip()))
    
arqvp = open("vp_irineu.txt",'r')
tvp = arqvp.readlines()
for i in range(len(tvp)):
    vvp.append(float(tvp[i].strip()))
    
arqvs = open("vs_irineu.txt",'r')
tvs = arqvs.readlines()
for i in range(len(tvs)):
    vvs.append(float(tvs[i].strip()))
    
def correlacao (xm, ym, xi, yi):
    sxx = 0
    syy = 0
    sxy = 0
    for i in range(len(xi)): sxx = sxx + (xi-xm)**2
    for i in range(len(xi)): syy = syy + (yi-ym)**2
    for i in range(len(xi)): sxx = sxx + (xi-xm)*(yi-ym)
    return sxy/((sxx*syy)**0.5)

def theta(alpha):
    return alpha*(np.arccos(alpha) - alpha*np.sqrt(1.0 - alpha*alpha))/(1.0 - alpha*alpha)**(3.0/2.0)

def f(alpha, theta):
    return alpha*alpha*(3.0*theta - 2.0)/(1.0 - alpha*alpha)

def PQ(A, B, R, theta, f):
    F1 = 1.0 + A*(1.5*(f + theta) - R*(1.5*f + 2.5*theta - 4.0/3.0))
    F2 = 1.0 + A*(1.0 + 1.5*(f + theta) - R*(1.5*f + 2.5*theta)) + B*(3.0 - 4.0*R) + A*(A + 3.0*B)*(1.5 - 2.0*R)*(f + theta - R*(f - theta + 2.0*theta*theta))
    F3 = 1.0 + A*(1.0 - f - 1.5*theta + R*(f + theta))
    F4 = 1.0 + (A/4.0)*(f + 3.0*theta - R*(f - theta))
    F5 = A*(-f + R*(f + theta - 4.0/3.0)) + B*theta*(3.0 - 4.0*R)
    F6 = 1.0 + A*(1.0 + f - R*(f + theta)) + B*(1.0 - theta)*(3.0 - 4.0*R)
    F7 = 2.0 + (A/4.0)*(3.0*f + 9.0*theta - R*(3.0*f + 5.0*theta)) + B*theta*(3.0 - 4.0*R)
    F8 = A*(1.0 - 2.0*R + (f/2.0)*(R - 1.0) + (theta/2.0)*(5.0*R - 3.0)) + B*(1.0 - theta)*(3.0 - 4.0*R)
    F9 = A*((R - 1.0)*f - R*theta) + B*theta*(3.0 - 4.0*R)
    
    P = 3.0*F1/F2
    Q = 2.0/F3 + 1.0/F4 + (F4*F5 + F6*F7 - F8*F9)/(F2*F4)
    return P, Q

def KG(Km, Gm, Ki, Gi, ci, theta, f):
    A = Gi/Gm - 1.0
    B = (Ki/Km - Gi/Gm)/3.0
    R = Gm/(Km + (4.0/3.0)*Gm)
    Fm = (Gm/6.0)*(9.0*Km + 8.0*Gm)/(Km + 2.0*Gm)
    
    P, Q = PQ(A, B, R, theta, f)
    
    K = Km - (Km + (4.0/3.0)*Gm)*ci*(Km - Ki)*P/3.0/(Km + (4.0/3.0)*Gm + ci*(Km - Ki)*P/3.0)
#    
#    def teste (x):
#        return (Ki - Km)*P/(1.- x)
#    from scipy.integrate import quad
#    K2 = Km + quad(teste, 0, 0.4)[0]

#    K = (Km*(Km+4./3*Gm) + 4./3*Gm*ci*(Ki - Km)*P/3.)/(Km+4./3*Gm- ci*(Ki - Km)*P/3.0)
#    K = Km + ((Km+4./3*Gm)*(ci*(Ki - Km)*P/3.0))/(Km+4./3*Gm-ci*(Ki - Km)*P/3.0)
    G = Gm - (Gm + Fm)*ci*(Gm - Gi)*Q/5.0/(Gm + Fm + ci*(Gm - Gi)*Q/5.0)
    
    return K, G
#    return G

def DEM(Km, Gm, Ki, Gi, alphai, phii, curphi=0.0, minalphaphiratio=100):
    ni = np.ceil(minalphaphiratio*phii/alphai).astype(np.int)
    ni=[100]
    dphii = phii/ni
#    print "ni...", ni
    indexes = np.repeat(np.arange(len(ni)), ni)
#    print dphii
    np.random.shuffle(indexes)
    
    thetai = theta(alphai)
    fi = f(alphai, thetai)
    
    K = np.empty(len(indexes))
    G = np.empty(len(indexes))
    phi = np.empty(len(indexes))
    
    K_ = Km
    G_ = Gm
    phi_ = curphi   

    for j in range(len(indexes)):
        i = indexes[j]
        ci = dphii[i]/(1.0 - phi_)
        phi_ += dphii[i]
        K_, G_ = KG(K_, G_, Ki[i], Gi[i], ci, thetai[i], fi[i])
        
        K[j] = K_
        G[j] = G_
        phi[j] = phi_
     
    return K, G, phi
    
def discrete_normal(zmin, zmax, n):
    from scipy.stats import norm
    z = np.linspace(zmin, zmax, n)
    x = np.empty(n + 1)
    x[0] = np.exp(-z[0]**2/2.0)/norm.cdf(z[0])
    x[-1] = np.exp(-z[-1]**2/2.0)/(1.0 - norm.cdf(z[-1]))
    x[1:-1] = (np.exp(-z[:-1]**2/2.0) - np.exp(-z[1:]**2/2.0))/(norm.cdf(z[1:]) - norm.cdf(z[:-1]))
    x /= np.sqrt(2.0*np.pi)

    p = np.empty(len(z) + 1)
    p[0] = norm.cdf(z[0])
    p[-1] = 1.0 - norm.cdf(z[-1])
    p[1:-1] = norm.cdf(z[1:]) - norm.cdf(z[:-1])
    
    return x, p

def mainsat():
    import matplotlib.pyplot as plt
    import Gassmann
    Vpm = 5.85 # km/s
    Vsm = 3.9 # km/s
    rhom = 2.65 # g/cm3
    rhof = 1.1 # g/cm3
    Gm= 32.03
    Km= 76.85
    rhom= 2.71
    Kf = 2.2
    Km= 76.85
    Gm= 32.03    
    rhom= 2.71
    for alpha in [0.01, 0.1, 0.5]:
        K, G, phi = DEM(Km, Gm, np.array([Kf]), np.array([0.0]), np.array([alpha]), np.array([0.99]))
        rho = (1.0 - phi)*rhom + phi*rhof #sat
        Ks = Gassmann.MKs(K, Km, Kf, phi) #sat
#        Ks = K #dry
        Ks = np.insert(Ks,0,Km)
        phi = np.insert(phi,0,0.)
        G = np.insert(G,0,Gm)
        rho = np.insert(rho,0,rhom)
        Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
        Vs = np.sqrt(G/rho)
        plt.plot(phi, Vp, 'b')
#    plt.grid()
#    plt.show()

def mainbiot(): #chapter 6.1 from mavko, ks diffrent from main1
    import matplotlib.pyplot as plt
    import Gassmann
    Vpm = 5.85 # km/s
    Vsm = 3.9 # km/s
    rhom = 2.65 # g/cm3
    rhof = 1.1 # g/cm3
    Gm= 32.03
    Km= 76.85
    rhom= 2.71
    Kf = 2.2
    Km= 76.85
    Gm= 32.03    
    rhom= 2.71
    
    for alpha in [0.01, 0.1, 0.5]:
        K, G, phi = DEM(Km, Gm, np.array([0.0]), np.array([0.0]), np.array([alpha]), np.array([0.99]))
        biotalpha = 1.
        biotrho0 = rhom
        biotrhofl = rhof
        Kfl = Kf
        Kfr = K
        K0 = Km
        biotmifr = G
        biotphi = phi
        a = (1.-biotphi-(Kfr/K0))
        P = ((1.-biotphi)*(a)*K0+biotphi*K0*Kfr/Kfl)/(a+biotphi*K0/Kfl)+(4.*biotmifr/3.)
        Q = (a*(biotphi*K0))/((a+(biotphi*K0/Kfl)))
        R = (K0*biotphi**2)/(a + biotphi*K0/Kfl)
        biotrho11 = (1.-biotphi)*biotrho0 - (1.-biotalpha)*biotphi*biotrhofl
        biotrho22 = biotalpha*biotphi*biotrhofl
        biotrho12 = (1.-biotalpha)*biotphi*biotrhofl
        biotrho   = biotrho0*(1.-biotphi) + biotrhofl*biotphi
        delta = P*biotrho22 + R*biotrho11 - 2.*Q*biotrho12
        Vp = ((delta + (delta**2 - 4.*(biotrho11*biotrho22 - biotrho12**2)*(P*R-Q**2))**0.5)/(2.*(biotrho11*biotrho22 - biotrho12**2)))**0.5
#        Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
        Vs = np.sqrt(G/biotrho)#        
        plt.plot(biotphi, Vp, 'b')
#        plt.plot(biotphi, Vs, 'b')
#    plt.grid()
#    plt.show()

def mainIrineuDem():
    import matplotlib.pyplot as plt
    rhof = 0.0012928 # g/cm3
    lines = []
    for i in range(len(vvp)):
        print i
        Km = bulk
        Gm = shear
        rhom = dens
        alphameso = macroaspect
        phimeso = macroporos
        phimicro = microporos
        e=10
        alphamicro=0.001
        while e>0.01:
            K, G, phi = DEM(Km[i], Gm[i], np.array([0.0]), np.array([0.0]), np.array([alphameso[i]]), np.array([phimeso[i]]))
            K, G, phi = DEM(K[-1], G[-1], np.array([0.0]), np.array([0.0]), np.array([alphamicro]), np.array([phimicro[i]]), phi[-1])

            rho = (1.0 - phi)*rhom[i]# + phi*rhof
    #        Ks = Gassmann.Ks(K, Km, Kf, phi) #sat
    #        rho = rhom
            Ks = K #dry
            Ks = np.insert(Ks,0,Km[i])
            phi = np.insert(phi,0,0.)
            G = np.insert(G,0,Gm[i])
            rho = np.insert(rho,0,rhom[i])
            Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
            Vs = np.sqrt(G/rho)
    
            alphamicro= alphamicro+0.00002
            e = np.abs(Vp[-1]-vvp[i])
        print Vp[-1], Vs[-1], alphameso, alphamicro
        newarq = open("Vp_calc_DEM.txt", 'w')
        lines.append (str(i+1)+" "+str(Vp[-1])+" "+str(Vs[-1])+ " " + str(alphameso[i]) + " "+str(alphamicro))
        newarq.writelines ('\n'.join(lines))
        plt.plot(phi, Vp, 'b')             # plot vp
        
    for i in range(len(macroporos)):
        
        if bulk[i] ==0.:
            macroporos.pop(i)
            bulk.pop(i)
        else: 
            plt.plot(macroporos[i]+microporos[i],vvp[i],'o')
#    print 'vpp', Vp[-1]
        
    plt.grid()
    plt.show()
    
    
def dem_log(rock, fluid, phit, rhob, inclusion1, inclusion2):
#    rock=(km,gm,rhom)
#    fluid= (kf,rhof)
#    phit= log
#    rhob = log
#    inclusion = (fraction, alpha)
#    rhob = (1-phi)*rhom + phi*rhof

    frac1 = inclusion1[0]
    alpha1 = inclusion1[1]
    poros1 = frac1*phit
    frac2 = inclusion2[0]
    alpha2 = inclusion2[1]
    poros2 = frac2*phit
#    fam1 = (fluid, 1./alpha1, poros1, 'cav')
#    fam2 = (fluid, 1./alpha2, poros2, 'cav')
#    alphad = 1./alphadual
    
    import matplotlib.pyplot as plt
    import Gassmann
#    rhom = rock[2] # g/cm3
#    Vpf = 1.6 # km/s
    rhof = fluid[1] # g/cm3
#    rhof = 0.2 # g/cm3
#    phimax = 0.6
    Gm= rock[1]
    Km= rock[0]
    rhom= rock[2]
    Kf = fluid[0]

    alphameso = alpha1
    phimeso = poros1
    alphamicro = alpha2
    phimicro = poros2
    
    for alpha in [0.1]:
##    for alpha in aspect:
##    for p in [0.2, 0.4, 0.6, 0.8,1]:
##        for alpha in [0.15, 0.02, 0.8]:
##        alphai = alpha
##        if alpha != 0.15:
##            alphai = p*alpha + (1-p)*0.15
#        K, G, phi = DEM(Km[i], Gm[i], np.array([0.0]), np.array([0.0]), np.array([alphameso[i]]), np.array([phimeso[i]]))
#            K, G, phi = DEM(K[-1], G[-1], np.array([0.0]), np.array([0.0]), np.array([alphamicro]), np.array([phimicro[i]]), phi[-1])
        K, G, phi = DEM(Km, Gm, np.array([0.0]), np.array([0.0]), np.array([alpha]), np.array([phit]))
##        print'kmmm \n', K,'\n\n', G, ' \n mm\n ', Km, Gm
#        
##        for i in range(len(phi)):
##            if alpha ==0.15: 
##                phi_ref = phi
##            else:
##                phi = 0.8*phi + 0.2*phi_ref
#        K, G, phi = DEM(Km, Gm, np.array([0.0]), np.array([0.0]), np.array([meso]), np.array([phimeso]))
    #    print 'KKK', Km, K,!
#        K, G, phi = DEM(K[-1], G[-1], np.array([0.0]), np.array([0.0]), np.array([micro]), np.array([phimicro]), phi[-1])
    
#        rho = (1.0 - phi)*rhom# + phi*rhof #dry
        rho = (1.0 - phi)*rhom + phi*rhof #sat
        Ks = Gassmann.Ks(K, Km, Kf, phi) #sat
#        rho = rhom
#        Ks = K #dry
        Ks = np.insert(Ks,0,Km)
        phi = np.insert(phi,0,0.)
        G = np.insert(G,0,Gm)
        rho = np.insert(rho,0,rhom)
        Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
        Vs = np.sqrt(G/rho)
        k_reuss = (phi/Kf + (1 - phi)/Km)**(-1)
        k_voigt = phi*Kf+ (1 - phi)*Km
#        print 'mai1', rho[-1], G[-1], phi[0], Vs[-1]
    
#            vp_ref.append(Vp)
#        
#            phi = p*phic + (1-p)*phir
#        if alpha != 0.15:
#        for prop in [0.2, 0.4, 0.6, 0.8]:
#            Vp = prop*vp_ref + (1-prop)*Vp
#        if alpha == ref:
#            plt.plot(phi, Vp, 'b')
#            
#        if alpha > ref:
#            plt.plot(phi, Vp, 'b--')
#        else:
#            plt.plot(phi, Vp, 'b-.')]
#        print '\nultimovp',Vp[-1]
#        plt.plot(phi, Vp, 'b')             # plot vp
#        alpha= alpha+0.0001
#        e = np.abs(Vp[-1]-3.186)
#        print e, alpha
#        plt.plot(phi, Vs, 'b')
#        plt.plot(phi, rho, 'g')
#        plt.subplot(211)
        crackdens = np.empty(len(phi))
        for i in range(len(phi)): 
            if phi[i]==0:
                crackdens[i] = 00.
            else:
                crackdens[i] = (3.*phi[i])/(4.*np.pi*alpha)
#        plt.plot(crackdens, Ks/(Ks[0]), 'g')
        plt.plot(phi, Ks/(Ks[0]), 'g')
#        plt.plot(phi, rho, 'g')
#        plt.plot(phi, Vp, 'g')
#        plt.plot(phi, Vs, 'g')
#        plt.grid()
#        print Ks
#        plt.subplot(212)
#        plt.plot(phi, Ks, 'g')
#        plt.plot(phi, G, 'g')
#        for i in range (len(Ks)):
#            if i !=0:
#                print (Ks[i] + 4.0*G[i]/3.0)-(Ks[i-1] + 4.0*G[i-1]/3.0)
        #(Ks + 4.0*G/3.0)
#        plt.plot(phi, k_reuss, 'g')
#        plt.plot(phi, (Ks + 4.0*G/3.0), 'b')       # plot c11
#        plt.plot(phi, k_voigt, 'g')

        
    for i in range(len(macroporos)):
        
        if bulk[i] ==0.:
#            i=1
#            print 'iii'
            macroporos.pop(i)
            bulk.pop(i)
            cont = cont+1
            print 'len',len(macroporos), len(bulk), cont
        else: 
            i=1
##            print '\nporos', macroporos[i]+microporos[i], vvp[i]
#            plt.plot(macroporos[i]+microporos[i],vvp[i],'o')
##            plt.plot(poros[i],vvp[i],'o')
#    print 'vpp', Vp[-1]
#    plt.xlim(0, 0.5)
    plt.grid(True)
#    plt.show()

def main1():
    import matplotlib.pyplot as plt
    import Gassmann
    Vpm = 5.85 # km/s
    Vsm = 3.9 # km/s
    rhom = 2.65 # g/cm3
#    Vpf = 1.6 # km/s
    rhof = 1.1 # g/cm3
#    rhof = 0.2 # g/cm3
    phimax = 0.6
    
#    Gm = Vsm*Vsm*rhom # GPa
#    Km = Vpm*Vpm*rhom - 4.0*Gm/3.0 # GPa
#    Kf = Vpf*Vpf*rhof # GPa
#    Gm= 30. #calcite
#    Km= 71.
    Gm= 32.03
    Km= 76.85
    rhom= 2.71
    
#    Km = 37.9#quartzo
#    Gm = 44.3
#    rhom = 2.65
    Kf = 2.2
#    Kf = 0.09
#    Gm2= 94.9 #dolomite
#    Km2= 38.8
#    rhom2= 2.92
#    Gm3= 22.9 #clay
#    Km3= 10.6
#    rhom3= 2.52
#    Gm = 0.5*(Gm1+Gm2)
#    k_reuss = (0.5/Km1 + 0.5/Km2)**(-1)
#    k_voigt = 0.5*Km1 + 0.5*Km2
#    Km = 0.5*(k_reuss + k_voigt)
#    rhom = 0.5*(rhom1+rhom2)
#    phi_ref=[]
##    vp_ref=[]
#    ref = 0.15
#    crack = 0.02
#    stiff = 0.8
#    crack2 = 0.2*crack + 0.8*ref
#    print crack2
#    crack4 = 0.4*crack + 0.6*ref
#    crack6 = 0.6*crack + 0.4*ref
#    crack8 = 0.8*crack + 0.2*ref
#    stiff2 = 0.2*stiff + 0.8*ref
#    stiff4 = 0.4*stiff + 0.6*ref
#    stiff6 = 0.6*stiff + 0.4*ref
#    stiff8 = 0.8*stiff + 0.2*ref
#    alphamix = [ref, crack, stiff, crack2, crack4, crack6, crack8, stiff2, stiff4, stiff6, stiff8]
#    for alpha in alphamix:#[0.15, 0.02, 0.8]:  #[0.01, 0.05, 0.1, 0.15] -- 15=ref, 02 = crack, 80=stiff
#    K = Km
#    G = Gm
#    plt.figure()
#    for alpha in [0.05, 0.1, 0.15, 0.2]:
#    Km =Km3
#    Gm=Gm3
#    rhom=rhom3
    phimax = 0.2303
    Km= 76.85
    Gm= 32.03    
    rhom= 2.71
    meso = 0.51
    micro = 0.04
#    porototal = 0.2303
    phimeso = 0.141
    phimicro = 0.0893
#    crack = 0.6122*meso + 0.3878*micro
#    for i in range(len(macroporos)):
#        Km= bulk[i]
#        Gm= shear[i]   
#        rhom= dens[i]
#        meso = macroaspect[i]
#        micro = microaspect[i]
#        phimeso = macroporos[i]
#        phimicro = microporos[i]
#    e=10
#    alpha=0.01
#    while e>0.001:
#        print e
#    for alpha in [0.99]:
#    for fluid in [0,1]:
#        alpha=0.99
#        if fluid ==1:
#            print fluid
#            rhof = 0.
#            Kf = 0.
#        print Kf, rhof
    for alpha in [0.1]:
##    for alpha in aspect:
##    for p in [0.2, 0.4, 0.6, 0.8,1]:
##        for alpha in [0.15, 0.02, 0.8]:
##        alphai = alpha
##        if alpha != 0.15:
##            alphai = p*alpha + (1-p)*0.15
        K, G, phi = DEM(Km, Gm, np.array([0.0]), np.array([0.0]), np.array([alpha]), np.array([0.99]))
##        print'kmmm \n', K,'\n\n', G, ' \n mm\n ', Km, Gm
#        
##        for i in range(len(phi)):
##            if alpha ==0.15: 
##                phi_ref = phi
##            else:
##                phi = 0.8*phi + 0.2*phi_ref
#        K, G, phi = DEM(Km, Gm, np.array([0.0]), np.array([0.0]), np.array([meso]), np.array([phimeso]))
    #    print 'KKK', Km, K,!
#        K, G, phi = DEM(K[-1], G[-1], np.array([0.0]), np.array([0.0]), np.array([micro]), np.array([phimicro]), phi[-1])
    
#        rho = (1.0 - phi)*rhom# + phi*rhof #dry
        rho = (1.0 - phi)*rhom + phi*rhof #sat
        Ks = Gassmann.Ks(K, Km, Kf, phi) #sat
#        rho = rhom
#        Ks = K #dry
        Ks = np.insert(Ks,0,Km)
        phi = np.insert(phi,0,0.)
        G = np.insert(G,0,Gm)
        rho = np.insert(rho,0,rhom)
        Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
        Vs = np.sqrt(G/rho)
        k_reuss = (phi/Kf + (1 - phi)/Km)**(-1)
        k_voigt = phi*Kf+ (1 - phi)*Km
#        print 'mai1', rho[-1], G[-1], phi[0], Vs[-1]
    
#            vp_ref.append(Vp)
#        
#            phi = p*phic + (1-p)*phir
#        if alpha != 0.15:
#        for prop in [0.2, 0.4, 0.6, 0.8]:
#            Vp = prop*vp_ref + (1-prop)*Vp
#        if alpha == ref:
#            plt.plot(phi, Vp, 'b')
#            
#        if alpha > ref:
#            plt.plot(phi, Vp, 'b--')
#        else:
#            plt.plot(phi, Vp, 'b-.')]
#        print '\nultimovp',Vp[-1]
#        plt.plot(phi, Vp, 'b')             # plot vp
#        alpha= alpha+0.0001
#        e = np.abs(Vp[-1]-3.186)
#        print e, alpha
#        plt.plot(phi, Vs, 'b')
#        plt.plot(phi, rho, 'g')
#        plt.subplot(211)
        crackdens = np.empty(len(phi))
        for i in range(len(phi)): 
            if phi[i]==0:
                crackdens[i] = 00.
            else:
                crackdens[i] = (3.*phi[i])/(4.*np.pi*alpha)
#        plt.plot(crackdens, Ks/(Ks[0]), 'g')
        plt.plot(phi, Ks/(Ks[0]), 'g')
#        plt.plot(phi, rho, 'g')
#        plt.plot(phi, Vp, 'g')
#        plt.plot(phi, Vs, 'g')
#        plt.grid()
#        print Ks
#        plt.subplot(212)
#        plt.plot(phi, Ks, 'g')
#        plt.plot(phi, G, 'g')
#        for i in range (len(Ks)):
#            if i !=0:
#                print (Ks[i] + 4.0*G[i]/3.0)-(Ks[i-1] + 4.0*G[i-1]/3.0)
        #(Ks + 4.0*G/3.0)
#        plt.plot(phi, k_reuss, 'g')
#        plt.plot(phi, (Ks + 4.0*G/3.0), 'b')       # plot c11
#        plt.plot(phi, k_voigt, 'g')

        
    for i in range(len(macroporos)):
        
        if bulk[i] ==0.:
#            i=1
#            print 'iii'
            macroporos.pop(i)
            bulk.pop(i)
            cont = cont+1
            print 'len',len(macroporos), len(bulk), cont
        else: 
            i=1
##            print '\nporos', macroporos[i]+microporos[i], vvp[i]
#            plt.plot(macroporos[i]+microporos[i],vvp[i],'o')
##            plt.plot(poros[i],vvp[i],'o')
#    print 'vpp', Vp[-1]
#    plt.xlim(0, 0.5)
    plt.grid(True)
#    plt.show()

def main2():
    import matplotlib.pyplot as plt
    import Gassmann
    
    Vpm = 5.85 # km/s
    Vsm = 3.9 # km/s
    rhom = 2.65 # g/cm3
    Vpf = 1.6 # km/s
    rhof = 1.1 # g/cm3
    phimax = 0.5
    
    Gm = Vsm*Vsm*rhom # GPa
    Km = Vpm*Vpm*rhom - 4.0*Gm/3.0 # GPa
    Kf = Vpf*Vpf*rhof # GPa

    # alpha_mean = [0.2] # [0.1, 0.2]
    # alpha_std = np.linspace(0.0025, 0.0525, 21) # np.linspace(0.0025, 0.0275, 11)
    alpha_mean = [0.5]
    alpha_std = np.linspace(0.0025, 0.1325, 14)    
    x, p = discrete_normal(-3.5, 3.5, 70001)
    phis = p*phimax
    for am in alpha_mean:
        
#        K, G, phi = DEM(Km, Gm, np.array([0.0]), np.array([0.0]), np.array([am]), np.array([phimax]), minalphaphiratio=10000)
        rho = (1.0 - phi)*rhom + phi*rhof
        Ks = Gassmann.Ks(K, Km, Kf, phi)
        Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
        Vs = np.sqrt(G/rho)
        plt.plot(phi, Vp, 'b')
        plt.plot(phi, Vs, 'g')
        
        for ast in alpha_std:
            alphas = x*ast + am
            K, G, phi = DEM(Km, Gm, np.zeros_like(phis), np.zeros_like(phis), alphas, phis, minalphaphiratio=10000)
            rho = (1.0 - phi)*rhom + phi*rhof
            Ks = Gassmann.Ks(K, Km, Kf, phi)
            Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
            Vs = np.sqrt(G/rho)
            plt.plot(phi, Vp, 'b')
            plt.plot(phi, Vs, 'g')
    
    plt.ylim([0.0, 6.0])
    plt.xlim([0.0, 0.5])
    plt.show()

def Poisson (Km, Gm):
    return (3.*Km-2.*Gm)/(2.*(3.*Km+Gm)) 

def Calcite ():
#    Kmcal = 71.
#    Gmcal = 30.0
    Kmcal = 76.85
    Gmcal = 32.03
    rhocal = 2.710
    Ccal = np.zeros((6,6))
    Ccal[0][0]=Ccal[1][1]=Ccal[2][2] = Kmcal +4.*Gmcal/3
    Ccal[3][3]=Ccal[4][4]=Ccal[5][5] = Gmcal
    Ccal[0][1]=Ccal[0][2]=Ccal[1][0]=Ccal[1][2]=Ccal[2][0]=Ccal[2][1] = Kmcal - 2.*Gmcal/3
#    phical = 0.2
#    alphacal = 0.2
    poiscal = Poisson(Kmcal, Gmcal)
#    (3.*Kmcal-2.*Gmcal)/(2.*(3.*Kmcal+Gmcal))    
    Ccal_ = CVoigttoKelvin(Ccal)
#    print 'antes kelvin \n', Ccal, '\n', Ccal_
#    print Ccal[0][0]**0.5
    return [Ccal_, rhocal, poiscal]

def Dolomite ():
    Kmdol = 94.9
    Gmdol = 45.0
    rhodol = 2.87
    Cdol = np.zeros((6,6))
    Cdol[0][0]=Cdol[1][1]=Cdol[2][2] = Kmdol +4.*Gmdol/3
    Cdol[3][3]=Cdol[4][4]=Cdol[5][5] = Gmdol
    Cdol[0][1]=Cdol[0][2]=Cdol[1][0]=Cdol[1][2]=Cdol[2][0]=Cdol[2][1] = Kmdol - 2.*Gmdol/3
    poisdol = Poisson(Kmdol, Gmdol)#(3.*Kmdol-2.*Gmdol)/(2.*(3.*Kmdol+Gmdol))
    Cdol_ = CVoigttoKelvin(Cdol)
    return [Cdol_, rhodol, poisdol]

def Aragonite ():
    Kmara = 44.8
    Gmara = 38.8
    rhoara = 2.92
    Cara = np.zeros((6,6))
    Cara[0][0]=Cara[1][1]=Cara[2][2] = Kmara +4.*Gmara/3
    Cara[3][3]=Cara[4][4]=Cara[5][5] = Gmara
    Cara[0][1]=Cara[0][2]=Cara[1][0]=Cara[1][2]=Cara[2][0]=Cara[2][1] = Kmara - 2.*Gmara/3
    poisara = Poisson(Kmara, Gmara)#(3.*Kmara-2.*Gmara)/(2.*(3.*Kmara+Gmara))
    Cara_ = CVoigttoKelvin(Cara)
    return [Cara_, rhoara, poisara]
    
def Quartzo ():
    Kmq = 37.9
    Gmq = 44.3
    rhoq = 2.65
    Cq = np.zeros((6,6))
    Cq[0][0]=Cq[1][1]=Cq[2][2] = Kmq +4.*Gmq/3
    Cq[3][3]=Cq[4][4]=Cq[5][5] = Gmq
    Cq[0][1]=Cq[0][2]=Cq[1][0]=Cq[1][2]=Cq[2][0]=Cq[2][1] = Kmq - 2.*Gmq/3
    poisq = Poisson(Kmq, Gmq)#(3.*Kmq-2.*Gmq)/(2.*(3.*Kmq+Gmq))
    Cq_ = CVoigttoKelvin(Cq)
    print Cq_[0][0]
    return [Cq_, rhoq, poisq]
    
def Clay ():

    rhoc = 2.52#1.48
    Cclay = np.zeros((6,6))
#    Cclay[0][0] = Cclay[1][1] = 17.15
#    Cclay[2][2] = 5.26
#    Cclay[0][1] = Cclay[1][0] = 17.15-2*6.63
#    Cclay[0][2] = Cclay[1][2] = Cclay[2][0] = Cclay[2][1] = 2.71
#    Cclay[5][5] = 6.63
#    Cclay[4][4] = Cclay[3][3] = 1.48    
#    v1 = 0.2
    Kmc = 22.9
    Gmc = 10.6
    Cclay[0][0]=Cclay[1][1]=Cclay[2][2] = Kmc +4.*Gmc/3
    Cclay[3][3]=Cclay[4][4]=Cclay[5][5] = Gmc
    Cclay[0][1]=Cclay[0][2]=Cclay[1][0]=Cclay[1][2]=Cclay[2][0]=Cclay[2][1] = Kmc - 2.*Gmc/3
    
    poisc = 0.4
    Cclay_ = CVoigttoKelvin(Cclay)
    return [Cclay_, rhoc, poisc]

def MineralMatriz (Km = 76.85, Gm = 32.03, rho = 2.710):
    Ccal = np.zeros((6,6))
    Ccal[0][0]=Ccal[1][1]=Ccal[2][2] = Km +4.*Gm/3
    Ccal[3][3]=Ccal[4][4]=Ccal[5][5] = Gm
    Ccal[0][1]=Ccal[0][2]=Ccal[1][0]=Ccal[1][2]=Ccal[2][0]=Ccal[2][1] = Km - 2.*Gm/3
    poiscal = Poisson(Km, Gm)
#    print 'minmatr', Ccal[0][0]
    Ccal_ = CVoigttoKelvin(Ccal)
#    print 'ccalal', Ccal,'\n', Ccal_
    return [Ccal_, rho, poiscal]

def FluidSat (Km = 2.2, rho = 1.1, visc=0.001):
    FluidS = np.zeros((6,6))
    FluidS[0][0]=FluidS[1][1]=FluidS[2][2] = Km
    FluidS[3][3]=FluidS[4][4]=FluidS[5][5] = 0
    FluidS[0][1]=FluidS[0][2]=FluidS[1][0]=FluidS[1][2]=FluidS[2][0]=FluidS[2][1] = Km
    poisfluid = 0.5
#    print 'minmatr', Ccal[0][0]
    FluidS_ = CVoigttoKelvin(FluidS)
#    print 'ccalal', Ccal,'\n', Ccal_
    return [FluidS_, rho, poisfluid, visc]
    
def Clayxx ():

    rhoc = 2.52#1.48
    Cclay = np.zeros((6,6))
    Cclay[0][0] = Cclay[1][1] = 17.15
    Cclay[2][2] = 5.26
    Cclay[0][1] = Cclay[1][0] = 17.15-2*6.63
    Cclay[0][2] = Cclay[1][2] = Cclay[2][0] = Cclay[2][1] = 2.71
    Cclay[5][5] = 6.63
    Cclay[4][4] = Cclay[3][3] = 1.48    
#    v1 = 0.2
#    Kmc = 22.9
#    Gmc = 10.6
#    Cclay[0][0]=Cclay[1][1]=Cclay[2][2] = Kmc +4.*Gmc/3
#    Cclay[3][3]=Cclay[4][4]=Cclay[5][5] = Gmc
#    Cclay[0][1]=Cclay[0][2]=Cclay[1][0]=Cclay[1][2]=Cclay[2][0]=Cclay[2][1] = Kmc - 2.*Gmc/3
    
    poisc = 0.4
#    Cclay_ = CVoigttoKelvin(Cclay)
    return [Cclay_, rhoc, poisc]


def Water ():
#    Kmw = 2.2
#    Gmw = 0.0
    rhow = 1.1
    Cw = np.zeros((6,6))
    Cw[0][0] = Cw[1][1] = 2.2
    Cw[2][2] = 2.2
    Cw[0][1] = Cw[1][0] = 2.2
    Cw[0][2] = Cw[1][2] = Cw[2][0] = Cw[2][1] = 2.2
    Cw[5][5] = 0
    Cw[4][4] = Cw[3][3] = 0
    viscow = 0.001 #1cp = 10^-3 Pa.s
#    phiwat = 0.2
#    alphawat = 0.2
    poisw = 0.5
    Cw_ = CVoigttoKelvin(Cw)
    return [Cw_, rhow, poisw, viscow]

def Methane ():
#    Kmw = 2.2
#    Gmw = 0.0
    rhomet = 0.2
    Cmet = np.zeros((6,6))
    Cmet[0][0] = Cmet[1][1] = 0.09
    Cmet[2][2] = 0.09
    Cmet[0][1] = Cmet[1][0] = 0.09
    Cmet[0][2] = Cmet[1][2] = Cmet[2][0] = Cmet[2][1] = 0.09
    Cmet[5][5] = 0
    Cmet[4][4] = Cmet[3][3] = 0
    viscomet = 0.029
#    phiwat = 0.2
#    alphawat = 0.2
    poismet = 0.5
    Cmet_ = CVoigttoKelvin(Cmet)
    return [Cmet_, rhomet, poismet, viscomet]

def Dry():
    rhod = 0.001
    viscod = 0.000001
    Cd = np.zeros((6,6))
    poisd = 0.5
#    Cd_ = CVoigttoKelvin(Cd)
    return [Cd, rhod, poisd, viscod]

def Ar():
    theta = np.pi/3
    b1 = 1
    b3 = 20
    alpha = b3/b1
#    A1 = (1./(4*np.pi))*((b1**2*b3)/(b1**2*np.sin(theta)**2 + b3**2*np.cos(theta)**2)**(1.5))
    A2 = (1./(4*np.pi))*((alpha)/(np.sin(theta)**2 + alpha**2*np.cos(theta)**2)**(1.5))
    print A2
    
def Dqs (C0, km, kn):
    return C0*km*kn

#def TensorGreen (alpha, C0):
#    C0 = KelvintoVoigt(C0)
#    Gtensor = E(alpha,C0)
#    Epqrs = Etensor[0]
#    Epqsr = Etensor[1]
#    Eqprs = Etensor[2]
#    Eqpsr = Etensor[3]
#    G = -0.25*(Epqrs + Epqsr + Eqprs + Eqpsr)
#    return G
#    return -Gtensor
    
#def Green (alpha, C):
##    E = [Epqrs, Epqsr, Eqprs, Eqpsr]   --> iguais !
##    G = -0.25*(E[0] + E[1] + E[2] + E[3])
#    
#    C11 = C[0][0] 
#    C44 = C[3][3]
#    C33 = C[2][2]
#    C12 = C[0][1]
#    C21 = C[1][0]
#    C13 = C[0][2]
#    Epqrs = E(alpha, C11, C44, C33, C12, C13)
#    G = -1.*Epqrs
#    return G
    
def TensorGreen (alpha, C0_):
    from scipy.integrate import quad
    C0 = CKelvintoVoigt(C0_)
    alpha = 1./alpha
    C11 = C0[0][0]
    C44 = C0[3][3]
    C33 = C0[2][2]
    C12 = C0[0][1]
    C13 = C0[0][2]
    d = C11
    f = C44
    h = C33
    e = (C11-C12)/2.
    g = C13 + C44
#    e=f
#    d=h
    alphaquad = alpha**-2
    def delta(x):
        temp=(e*(1.-x**2) + f*(alphaquad)*(x**2))*((d*(1.-x**2) + f*alphaquad*(x**2))*(f*(1.-(x**2)) + h*alphaquad*(x**2)) - (g**2)*alphaquad*(x**2)*(1.-(x**2)))
        return temp**(-1)
    
    def e1111 (x):
        return (delta(x)*(1.-(x**2)))*((f*(1.-(x**2)) + h*alphaquad*(x**2))*((3.*e+d)*(1.-(x**2)) + 4*f*alphaquad*(x**2)) - 
        (g**2)*alphaquad*(x**2)*(1.-(x**2)))
        
    def e3333 (x):
        return delta(x)*alphaquad*(x**2)*(d*(1.-(x**2)) + f*alphaquad*(x**2))*(e*(1.-(x**2))+f*alphaquad*(x**2))
#    
    def e1122 (x):
        return delta(x)*(1.-(x**2))*((f*(1.-(x**2)) + h*alphaquad*(x**2))*((e+3.*d)*(1.-(x**2)) + 4.*f*alphaquad*(x**2)) -
        3.*(g**2)*alphaquad*(x**2)*(1.-(x**2)))
        
    def e1133 (x):
        return delta(x)*alphaquad*(x**2)*(((d+e)*(1.-(x**2)) + 2.*f*alphaquad*(x**2))*(f*(1.-(x**2)) + h*alphaquad*(x**2)) -
        (g**2)*alphaquad*(x**2)*(1.-(x**2)))
    
    def e3311 (x):
        return delta(x)*(1.-(x**2))*(d*(1.-(x**2)) + f*alphaquad*(x**2))*(e*(1.-(x**2)) + f*alphaquad*(x**2))
    
    def e1212 (x):
        return delta(x)*((1.-(x**2))**2)*((g**2)*alphaquad*(x**2)-(d-e)*(f*(1.-(x**2)) + h*alphaquad*(x**2)))
    
    def e1313 (x):
        return delta(x)*g*alphaquad*(x**2)*(1.-(x**2))*(e*(1.-(x**2)) + f*alphaquad*(x**2))
#    
#    def e1212 (x):#primeira 1212=e2121
#        return delta(x)*(1-x**2)*((f*(1-x**2) + h*alphaquad*x**2)*((e+3*d)*(1-x**2) + 4*f*alphaquad*x**2) -
#        3*g**2*alphaquad*x**2*(1-x**2))
##    
#    def e1313 (x):#1133 vira 1313=2323
#        return delta(x)*alphaquad*(x**2)*(((d+e)*(1.-(x**2)) + 2.*f*alphaquad*(x**2))*(f*(1.-(x**2)) + h*alphaquad*(x**2)) -
#        (g**2)*alphaquad*(x**2)*(1.-(x**2)))
##        
#    def e3131 (x):#3311 vira 3131=3232
#        return delta(x)*(1.-(x**2))*(d*(1.-(x**2)) + f*alphaquad*(x**2))*(e*(1.-(x**2)) + f*alphaquad*(x**2))
##    
#    def e1122 (x):#1212 vira 1122
#        return delta(x)*((1.-(x**2))**2)*((g**2)*alphaquad*(x**2)-(d-e)*(f*(1.-(x**2)) + h*alphaquad*(x**2)))
##    
#    def e1133 (x):#1313 vira 1133=2233
#        return delta(x)*g*alphaquad*(x**2)*(1.-(x**2))*(e*(1.-(x**2)) + f*alphaquad*(x**2))
    
        
    E1111 = (0.5*np.pi)*quad(e1111, 0, 1)[0]                   #E1111/E2222
    E3333 = (4.*np.pi)*quad(e3333, 0, 1)[0]                     #E3333
    E1122 = (0.5*np.pi)*quad(e1122, 0, 1)[0]                   #E1122/E2211
    E1133 = (2.*np.pi)*quad(e1133, 0, 1)[0]                     #E1133/E2233
    E3311 = (2.*np.pi)*quad(e3311, 0, 1)[0]                     #E3311/E3322
    E1212 = (0.5*np.pi)*quad(e1212, 0, 1)[0]                    #E1212
    E1313 = (-2.*np.pi)*quad(e1313, 0, 1)[0]                    #E1313/E2323
    
#    E1212 = (0.5*np.pi)*quad(e1212, 0, 1)[0]                   #E1212/2121
#    E1313 = (2.*np.pi)*quad(e1313, 0, 1)[0]                     #E1313/2323
#    E3131 = (2.*np.pi)*quad(e3131, 0, 1)[0]                     #E3131/3232
#    E1122 = (0.5*np.pi)*quad(e1122, 0, 1)[0]                    #1122
#    E1133 = (-2.*np.pi)*quad(e1133, 0, 1)[0]                    #E1133/2233
    
#    Epqrs = np.zeros((6,6))
#    Epqrs[0][0] = Epqrs[1][1] = E1111
#    Epqrs[2][2] = E3333
#    Epqrs[0][1] = Epqrs[1][0] = E1122
#    Epqrs[0][2] = Epqrs[1][2] = E1133
#    Epqrs[2][0] = Epqrs[2][1] = E3311
#    Epqrs[5][5] = E1212
#    Epqrs[4][4] = Epqrs[3][3] = E1313
    
    #    Epqrs = np.zeros((6,6))
#    Epqrs[0][0] = Epqrs[1][1] = E1111
#    Epqrs[2][2] = E3333
#    Epqrs[0][1] = Epqrs[1][0] = E1212
#    Epqrs[0][2] = Epqrs[1][2] = E1313
#    Epqrs[2][0] = Epqrs[2][1] = E3131
#    Epqrs[5][5] = E1122
#    Epqrs[4][4] = E1133
    print '\neee\n', E1111,E3333,E1122,E1133,E3311,E1212,E1313,'\n\n'
    Epqrs = np.zeros((3,)*4)#np.tensordot(np.zeros((3,3)), np.zeros((3,3)),0)    
    Epqrs[0][0][0][0] = Epqrs[1][1][1][1] = E1111
    Epqrs[2][2][2][2] = E3333
    Epqrs[0][0][1][1] = Epqrs[1][1][0][0] = E1122
    Epqrs[0][0][2][2] = Epqrs[1][1][2][2] = E1133
    Epqrs[2][2][0][0] = Epqrs[2][2][1][1] = E3311
    Epqrs[0][1][0][1] = E1212
    Epqrs[0][2][0][2] = Epqrs[1][2][1][2] = E1313
#    
#    Epqrs[0][1][0][1] = Epqrs[1][0][1][0] = E1212
#    Epqrs[0][2][0][2] = Epqrs[1][2][1][2] = E1313
#    Epqrs[2][0][2][0] = Epqrs[2][1][2][1] = E3131
#    Epqrs[0][0][1][1] = E1122
#    Epqrs[0][0][2][2] = Epqrs[1][1][2][2] = E1133
#    Epqrs__ = TensortoVoigt(Epqrs)
#    Epqrs_ = CVoigttoKelvin(Epqrs__)
    Epqrs_ = CTensortoKelvin(Epqrs)
    
#    Epqsr = np.zeros((6,6))
#    Epqsr[0][0] = Epqsr[1][1] = E1111
#    Epqsr[2][2] = E3333
#    Epqsr[0][1] = Epqsr[1][0] = E1122
#    Epqsr[0][2] = Epqsr[1][2] = E1133
#    Epqsr[2][0] = Epqsr[2][1] = E3311
#    Epqsr[5][5] = E1212 
#    Epqsr[4][4] = Epqsr[3][3] = E1313 
#    
#    Eqprs = np.zeros((6,6))
#    Eqprs[0][0] = Eqprs[1][1] = E1111
#    Eqprs[2][2] = E3333
#    Eqprs[0][1] = Eqprs[1][0] = E1122
#    Eqprs[0][2] = Eqprs[1][2] = E1133
#    Eqprs[2][0] = Eqprs[2][1] = E3311
#    Eqprs[5][5] = E1212 
#    Eqprs[4][4] = Eqprs[3][3] = E1313 
#    
#    Eqpsr = np.zeros((6,6))
#    Eqpsr[0][0] = Eqpsr[1][1] = E1111
#    Eqpsr[2][2] = E3333
#    Eqpsr[0][1] = Eqpsr[1][0] = E1122
#    Eqpsr[0][2] = Eqpsr[1][2] = E1133
#    Eqpsr[2][0] = Eqpsr[2][1] = E3311
#    Eqpsr[5][5] = E1212 
#    Eqpsr[4][4] = Eqpsr[3][3] = E1313 
    
#    return -0.25*[Epqrs_, Epqrs_, Epqrs_, Epqrs_]
    return -Epqrs_

def Cinv(C0):
    """retorna S0 ao passar C0"""
    c11 = C0[0][0]
    c12 = C0[0][1]
    S0 = np.zeros((6,6))
    S0[0][0]=S0[1][1]=S0[2][2] = (c11+c12)/((c11-c12)*(c11+2*c12))
    Ccal[3][3]=Ccal[4][4]=Ccal[5][5] = Gm
    Ccal[0][1]=Ccal[0][2]=Ccal[1][0]=Ccal[1][2]=Ccal[2][0]=Ccal[2][1] = Km - 2.*Gm/3
    poiscal = Poisson(Km, Gm)
    Ccal_ = CVoigttoKelvin(Ccal)
#    print 'ccalal', Ccal,'\n', Ccal_
    return S0
    
def TensorG (C0, alpha, v):
    #alpha eh razao de aspecto menor 1
#    S0 = C0**-1
#    if alpha == 1:
#        q = 2/3
#    alpha=1./alpha
    alphaquad = alpha**(2)
#    print alphaquad, alpha
#    C0 = KelvintoVoigt(C0)
#    print C0
    S0 = np.linalg.inv(C0)
#    S0 = VoigttoKelvin(np.linalg.inv(KelvintoVoigt(C0)))
#    print 'alpaaa', alpha
#    alpha = 0.1
#    v=0.5
#    print 'qqq', S0, SVoigttoKelvin(np.linalg.inv(CKelvintoVoigt(C0)))
    if alpha <1.:
        q = ((alpha)/((1.-alphaquad)**(1.5)))*(np.arccos(alpha) - alpha*((1.-alphaquad)**0.5))
    elif alpha ==1.:
        q=2./3.
    else:
        q = ((alpha)/((alphaquad-1.)**(1.5)))*(alpha*((alphaquad-1.)**0.5) - np.arccosh(alpha))
#    print q
    S1111 = ((3.*alphaquad)/(8.*(1.-v)*(alphaquad-1.))) + (1./(4*(1-v))) * (1.-2.*v-(9./(4*(alphaquad-1.))))*q
    S3333 = (1./(2.*(1.-v)))*(1.-2.*v + (3.*alphaquad-1.)/(alphaquad-1.)-(1.-2.*v+(3.*alphaquad)/(alphaquad-1.))*q)
    S1122 = (1./(4.*(1.-v)))*(((alphaquad)/(2.*(alphaquad-1.)))-(1.-2.*v+((3.)/(4.*(alphaquad-1.))))*q)
    S1133 = (1./(2.*(1.-v)))*((-alphaquad/(alphaquad-1.))+0.5*(((3.*alphaquad)/(alphaquad-1.))-(1.-2.*v))*q)
    S3311 = (1./(2.*(1.-v)))*(2.*v-1.-(1./(alphaquad-1.))+(1.-2.*v+(3./(2*(alphaquad-1.))))*q)
    S1212 = (1./(4.*(1.-v)))*((alphaquad/(2.*(alphaquad-1)))+(1.-2.*v-(3./(4.*(alphaquad-1))))*q)
    S1313 = (1./(4.*(1.-v)))*(1.-2.*v-((alphaquad+1.)/(alphaquad-1.))-0.5*(1.-2.*v-((3*(alphaquad+1))/(alphaquad-1.)))*q)
    
    Sr = np.zeros((3,)*4)#np.tensordot(np.zeros((3,3)), np.zeros((3,3)),0)
    
    Sr[0][0][0][0] = Sr[1][1][1][1] = S1111             #S1111/S2222
    Sr[2][2][2][2] = S3333                              #S3333
    Sr[0][0][1][1] = Sr[1][1][0][0] = S1122             #S1122/S2211
    Sr[0][0][2][2] = Sr[1][1][2][2] = S1133             #S1133/S2233
    Sr[2][2][0][0] = Sr[2][2][1][1] = S3311             #S3311/S3322
    Sr[0][1][0][1] = S1212                              #S1212
    Sr[0][2][0][2] = Sr[1][2][1][2] = S1313             #S1313/S2323
    
    Sr_ = STensortoKelvin(Sr)
#    Sr_ = TensortoVoigt(Sr)
#    Sr__ = SVoigttoKelvin(Sr_)
#    Sr_[0][0]=Sr_[2][2]
#    Sr_[1][1]=Sr_[2][2]
#    Sr_[3][3]=Sr_[5][5]
#    Sr_[4][4]=Sr_[4][4]*4.
#    Sr_[5][5]=Sr_[5][5]*4.
#    Sr_[6][6]=Sr_[6][6]*4.
#    print '\ns111\n', S0, SKelvintoVoigt(S0), SVoigttoKelvin(SKelvintoVoigt(S0))
#    S0 = SVoigttoKelvin(S0)
    Gr = - np.dot(Sr_,S0)
#    Sr__ = KelvintoVoigt(Sr_)
    return Gr,Sr_ #ambo em voigt
    


def Tma (C0,  Nr, tsat=0, omega = 30.*np.pi , tal=0.00001, kfl = 30., visc = 0.001, perm = np.identity(3)):    
    alpha = Nr[0]
    poros = Nr[1]
    Cr = Nr[2]
    Gr = Nr[3]
#    I4_ = Id4()
#    I4 = Kelvin(I4_)
    S0 = np.linalg.inv(C0)
#    for i in range (len(C0)):
#        for j in range(len(C0[0])):
##            print C0[i][j]
#            if C0[i][j] != 0:
#                S0[i][j] = C0[i][j]**-1
#    I2 = 0
#    Gr = Green (Nr[0], Nr[2])   # aspect ratio rth, C rth
    
    tdry = np.dot((Cr-C0),np.linalg.inv((I4 - np.dot(Gr,(Cr-C0)))))
#    print 'tdry \n', C0, '\n', tdry
#    print tdry
    if tsat != 0:
        I2x = np.tensordot(np.identity(3), np.identity(3),0)
        I2 = Kelvin(I2x)
        kfl = Cr
        print 'tsat', kfl, 'kfl \n'
        
#        Kelvin(I2x)
        kdry = Kdry (Gr, C0, S0, I4)
        gamma = Gamma(kfl, kdry, S0)
        
        theta = Theta(kfl, S0, omega, visc, poros, gamma, tal, kdry, perm)
        zeta = Zeta(tdry, S0, poros, omega, tal, gamma, I2)
        qui = Qui(tdry, S0, I2)
#        print 'tasat', (np.complex(1, omega*tal*gamma))
        
        print 'tsat\n\n', qui
        tsat = (np.dot(theta,zeta) + np.complex(0,1)*omega*tal*np.dot(kfl,qui))/(1. + np.complex(0, 1)*omega*tal*gamma)
    
    tma = tdry + tsat  #green rth - inversa!
#    print 'tdry \n', np.real(tsat), '\n', tdry
    return tma
#    tsat = tdry + (Theta*Zeta(r) + np.complex(0, omega*tal*kfl*Qui(r)))/(1 + np.complex(0,omega*tal*Gamma(r)))
    
def Kdry (Gr, C0, S0):
    I4 = np.identity(6)
#    S0 = np.linalg.inv(C0)
    temp1 = np.dot(Gr,C0)
    temp2 = I4 - temp1
    invtemp2 = np.linalg.inv(temp2)
    kdry = np.dot(invtemp2,S0)
#    print 'printtt \n', Gr, '\n\n',C0, '\n\n', np.matrix(Gr)*np.matrix(C0), '\n\n', kdry
    return kdry
    
def Gamma(kfl, kdry, S0):
#    print 'printtt', kfl, kdry, S0
    temp1 = kdry - S0
    temp1uuvv = Muv(temp1)
#    gamma = 1. + np.dot(kfl,temp1)
    gamma = 1. + kfl*temp1uuvv
    return gamma
        
def Theta (kfl, S0, omega, visc, poros, gamma, tal, kdry, ku, kv, perm):
#    ku=omega/2650
#    kv=omega/2650    
#    ku=np.array([omega/2650.,0,0])    # ???
#    kv=np.array([0, omega/2650.,0])   # ???
#    iotg = omega*tal*gamma*np.complex(0,1)+1
#    print 'iotg \n', iotg
   # print 'prin', np.outer(np.outer(ku,kv), perm), '\n\n print', np.outer(np.identity(3),np.identity(3))
#    kuv = np.tensordot(ku, kv, 0)
#    kuvperm_ = np.tensordot(kuv, perm,0)
    
#    kuvperm = Kelvin (kuvperm_)
#    print '\n\n agora', np.dot(kfl,((poros*kdry)/(1.+ np.complex(0,1)*omega*tal*gamma))) - np.complex(0,1)*np.dot(kfl,kuvperm)/(omega*visc), '\ntheta'
#    theta = np.linalg.inv(np.dot(kfl,((1. - np.dot(kfl,S0))*((poros)/(1. + np.complex(0,1)*omega*tal*gamma)) + 
#    np.dot(kfl,((poros*kdry)/(1.+ np.complex(0,1)*omega*tal*gamma))) - np.complex(0,1)*np.dot(kfl,kuvperm)/(omega*visc))))
#    kfl = kf[0][0]
    Suuvv = Muv(S0)
    Kduuvv = Muv(kdry)
#    temp1 = 1. - np.dot(kfl,S0)
    temp1 = 1. - kfl*Suuvv #nao usada mais apos correção do modelo
    temp2 = poros/(1.+np.complex(0,1)*omega*tal*gamma)
    temp3 = (poros*Kduuvv)/(1.+np.complex(0,1)*omega*tal*gamma)
#    temp4 = np.complex(0,1)*np.dot(ku, kv)
    temp4 = np.dot(ku, perm)
#    print 'kuuu\n',ku, '\n', perm, '\n', kv, '\n', np.dot(np.dot(perm, ku), kv)
#    temp5 = np.dot(temp4, perm)
    temp5 = np.complex(0,1)*np.dot(temp4, kv)
#    temp6 = np.dot(temp5, kfl)
    temp6 = temp5*kfl
#    temp7 = np.dot(temp1, temp2)
    temp7 = temp1*temp2
#    temp8 = np.dot(kfl, temp3)
    temp8 = kfl*temp3
    temp9 = temp6/(omega*visc)
    temp10 = temp7 + temp8 - temp9
#    invtemp10 = np.linalg.inv(temp10)
    invtemp10 = (temp10)**-1
#    theta = np.dot(kfl, invtemp10)   
    theta = kfl*invtemp10
#    print theta
    
    delta = (kfl*tal)/(poros*visc)*np.dot(np.dot(perm, ku), kv)
#    print delta, theta
#    temp11 = poros*gamma/(1+np.complex(0,1)*omega*tal*gamma)
#    temp12 = np.complex(0,1)*kfl*np.dot(np.dot(perm, ku), kv)/(visc*omega*(1-delta))
#    ntheta = temp11 + temp12
#    nthetainv = (ntheta)**-1
#    nktheta = kfl*nthetainv
#    print '\n ', nktheta
    return theta*(1-delta)
    
def Zeta (tdry, S0, poros, omega, tal, gamma):
    I2 = np.identity(3)
    i2xi2 = np.tensordot(I2,I2,0)
    i2xi2m = CTensortoKelvin(i2xi2)
#    zeta = np.dot(np.dot(np.dot(tdry,S0),(I2x)),S0)*((poros*tdry)/(1. + np.complex(0,1)*omega*tal*gamma))
    temp1 = np.dot(tdry,S0)
#    temp2 = np.dot(I2, I2) ##### VERIFICAR PRODUTO DIADICO
    temp3 = np.dot(temp1, i2xi2m)
    temp4 = np.dot(temp3, S0)
    temp5 = (poros*tdry)/(1.+np.complex(0,1)*omega*tal*gamma)
    zeta = np.dot(temp4, temp5)
    return zeta
    
def Qui (tdry, S0):
#    qui = np.dot(np.dot(np.dot(np.dot(tdry,S0),(I2x)),S0),tdry)
    I2 = np.identity(3)
    i2xi2 = np.tensordot(I2,I2,0)
    i2xi2m = CTensortoKelvin(i2xi2)
    temp1 = np.dot(tdry,S0)
#    temp2 = np.dot(I2, I2) ##### VERIFICAR PRODUTO DIADICO
    temp3 = np.dot(temp1, i2xi2m)
    temp4 = np.dot(temp3, S0)
    qui = np.dot(temp4, tdry)
    return qui

def Id4 ():
    I4_ = np.tensordot(np.zeros((3,3)), np.zeros((3,3)),0)
    for j in range (len(I4_)):
        for k in range(len(I4_[0])):
            for l in range (len(I4_[0][0])):
                for m in range (len(I4_[0][0][0])):
                    if j==l and k==m:
                        I4_[j][k][l][m] = I4_[j][k][l][m] + 0.5
                    if j==m and k==l:
                        I4_[j][k][l][m] = I4_[j][k][l][m] + 0.5
    return I4_
    
def TI (c11, c12, c13, c33, c44):
    m=np.zeros((6,6))
    m[0][0] = m[1][1] = c11
    m[0][1] = m[1][0] = c12
    m[0][2] = m[1][2] = m[2][0] = m[2][1] = c13
    m[2][2] = c33
    m[3][3] = m[4][4] = c44
    m[5][5] = 0.5*(c11-c12)
    
    return m    

def TensortoVoigt(k):                       #notacao de voigt para tensor de 4a ordem
    m=np.zeros((6,6))
    m[0,0]=k[0][0][0][0]
    m[0,1]=k[0][0][1][1]
    m[0,2]=k[0][0][2][2]
    m[0,3]=k[0][0][1][2]
    m[0,4]=k[0][0][0][2]
    m[0,5]=k[0][0][0][1]
    m[1,0]=k[1][1][0][0]
    m[1,1]=k[1][1][1][1]
    m[1,2]=k[1][1][2][2]
    m[1,3]=k[1][1][1][2]
    m[1,4]=k[1][1][0][2]
    m[1,5]=k[1][1][0][1]
    m[2,0]=k[2][2][0][0]
    m[2,1]=k[2][2][1][1]
    m[2,2]=k[2][2][2][2]
    m[2,3]=k[2][2][1][2]
    m[2,4]=k[2][2][0][2]
    m[2,5]=k[2][2][1][0]
    m[3,0]=k[1][2][0][0]
    m[3,1]=k[1][2][1][1]
    m[3,2]=k[1][2][2][2]
    m[3,3]=k[1][2][1][2]
    m[3,4]=k[1][2][0][2]
    m[3,5]=k[1][2][0][1]
    m[4,0]=k[0][2][0][0]
    m[4,1]=k[0][2][1][1]
    m[4,2]=k[0][2][2][2]
    m[4,3]=k[0][2][1][2]
    m[4,4]=k[0][2][0][2]
    m[4,5]=k[0][2][0][1]
    m[5,0]=k[0][1][0][0]
    m[5,1]=k[0][1][1][1]
    m[5,2]=k[0][1][2][2]
    m[5,3]=k[0][1][1][2]
    m[5,4]=k[0][1][0][2]
    m[5,5]=k[0][1][0][1]
    return m 

def SVoigttoKelvin (v):              #conversao de notacao voigt para kelvin
    i = np.identity(6)
    i[3][3] = i[4][4] = i[5][5] = np.sqrt(2)
    i_ = np.linalg.inv(i)
    temp1 = np.dot(i_,v)
    temp2 = np.dot(temp1,i_)
    return temp2

def SKelvintoVoigt (k):              #conversao de notacao kelvin para voigt
    i = np.identity(6)
    i[3][3] = i[4][4] = i[5][5] = np.sqrt(2)
    temp1 = np.dot(i,k)
    temp2 = np.dot(temp1,i)
    return temp2
    
def CVoigttoKelvin (v):              #conversao de notacao voigt para kelvin
    i = np.identity(6)
    i[3][3] = i[4][4] = i[5][5] = np.sqrt(2)
    temp1 = np.dot(i,v)
    temp2 = np.dot(temp1,i)
    return temp2

def CKelvintoVoigt (k):              #conversao de notacao kelvin para voigt
    i = np.identity(6)
    i[3][3] = i[4][4] = i[5][5] = np.sqrt(2)
    i_ = np.linalg.inv(i)
    temp1 = np.dot(i_,k)
    temp2 = np.dot(temp1,i_)
    return temp2

def CTensortoKelvin (t):             #conversao de tensor de 4a ordem para kelvin
    temp1 = TensortoVoigt(t)
    temp2 = CVoigttoKelvin(temp1)
    return temp2

def STensortoKelvin (t):             #conversao de tensor de 4a ordem para kelvin
    temp1 = TensortoVoigt(t)
    temp2 = SVoigttoKelvin(temp1)
    return temp2

def Tuuvv (T):                      #traco do tensor na forma tensorial
    Sum = 0.
    for i in range(3):
        for j in range (3):
            Sum += T[i][i][j][j]
    return Sum

def Muv (M):                        ##traco do tensor na forma matricial Kelvin
    Sum = 0.
    for i in range(3):
        for j in range (3):
            Sum += M[i][j]
    return Sum
    
def DeltaC (ci, c0):
    return ci - c0

def T_r (deltaC, Gr):
    I4 = (np.identity(6))
#    I4 = CVoigttoKelvin((np.identity(6)))
#    I4 = KelvintoVoigt(np.identity(6))
#    print Id4(),'id4'
    temp = I4 - np.dot(Gr, deltaC)
    
#    print ' \ntrrrr1',deltaC,'\n\n', np.linalg.det(temp)
    invtemp = np.linalg.inv(temp)   
#    print 'fimter\n'
    tr = np.dot(deltaC,invtemp)
#    print I4
#    tempx1 = I4 - np.dot(deltaC, Gr)
#    tempx2 = np.linalg.inv(tempx1)
#    tempx3 = np.dot(tempx2, deltaC)
    return tr
    
def T1 (vr, tr):
    T1 = np.dot(vr,tr)
#    print 't', vr, tr[0][0], T1[0][0]
    return T1
        
def T2 (vr, tr, gdrs, vs, ts):
    T1_ = T1(vr, tr)
    temp1 = T1_
    temp2 = np.dot(temp1,gdrs)
    temp3 = np.dot(temp2, ts)
    temp4 = np.dot(temp3, vs)
    return temp4
    
def IncludeT (C0, Cr, Gr, poison, alpha, vr, k, omega, rhom, viscr):  #ni = number of inclusions
#    print 'kmoega\n', k, omega
    if omega==0.: omega=1.
    kfl = Cr[0][0]    
#    S0 = np.linalg.inv(C0)
#    print 'omegaa3', omega
#    freq = 100
#    omega = 2*np.pi*freq
#    omega = 100000000
    vpm = ((C0[0][0]+4.*C0[3][3]/3)/rhom)**0.5
    ku = np.array([omega/vpm,omega/vpm,omega/vpm]) 
#    ku = omega/vpm
    kv = np.array([omega/vpm,omega/vpm,omega/vpm])
    tal = 0.00001
#    visc = 0.001
    perm = k
    
    dry = Dry()
    Cdry = dry[0]
    dC_dry = DeltaC(Cdry, C0)
    G_dry, Sr = TensorG(C0, alpha, poison) 
    
    S0 = np.linalg.inv(C0)
#    G1_ = TensorG(C0, alpha, poisi)
    tdry = T_r(dC_dry, G_dry)
    
    kdry = Kdry (G_dry, C0, S0)
    
    gamma = Gamma(kfl, kdry, S0)
#    print 'freqqq\n',omega
    theta = Theta(kfl, S0, omega, viscr*(10**-9), vr, gamma, tal, kdry, ku/1000, kv/1000, perm*(10**-12))
#    theta = Theta(kfl, S0, omega, viscr, vr, gamma, tal, kdry, ku, kv, perm)
    
    zeta = Zeta(tdry, S0, vr, omega, tal, gamma)
    qui = Qui(tdry, S0)
    temp1 = np.dot(theta, zeta)
    temp2 = np.complex(0,1)*omega*tal*np.dot(kfl, qui)
#    print 't', theta
    temp3 = 1.+np.complex(0,1)*omega*tal*gamma
#    print 't', temp1, '\n', temp2, '\n', temp3
    tsat = tdry + (temp1 + temp2)/(temp3)
#    print  '\ntttsat\n', tsat, '\ntdry\n', (temp1 + temp2)/(temp3), vpm, viscr
    return tsat
    
def IncludeDem (mineral1, fluido2, alpha, phimax, curphi, ni):  #ni = number of inclusions
    fluido2 = Water()
#    fluido2 = Dry()
    C0 = mineral1[0] 
    rhom = mineral1[1]
    poisi = mineral1[2]
    C1 = fluido2[0]
    rhoi = fluido2[1]
#    poisi = fluido2[2]
#    G1 = TensorG(C0, alpha, poisi)
#    print G1
    step = phimax/ni
    phi_ = curphi
    ts = np.zeros((6,6))
    vs = np.zeros((6,6))
    phi = np.empty(ni)
    rho = np.empty(ni)
#    Vp = np.empty(ni)
#    Vs = np.empty(ni)
    Ct = np.empty((ni,6,6))
    for i in range (ni):
#        curphi = phi[i]
        phi_ += step
        G1 = TensorG(C0, alpha, poisi)
#        dC = DeltaC(C1, C0)
        tr = IncludeT (C0, C1, G1, alpha, phi_)        
#        print tr
#        tdry = 
#        tma = tdry + tsat
        Ct_ = Ctotal (C0, tr, step, G1, ts, vs)
        Ct[i] = Ct_
        C0 = Ct_
        ts = tr
        vs = phi_
        
        phi[i] = phi_
        
        rho[i] = (1.0 - phi[i])*rhom + phi[i]*rhoi
#        print rho[i]
#        Vp[i] = ((Ct[i][0][0]/rho)**-0.5)**-1
#        Vs[i] = ((Ct[i][3][3]/rho)**-0.5)**-1
    Gct = Ct[-1][3][3]*0.5
    Kct = Ct[-1][0][0] - 4.*Gct/3.
    poisct = (3.*Kct-2.*Gct)/(2.*(3.*Kct+Gct))
#    print Gct, Kct,poisct, Ct[-1]
    return [Ct, rho, phi, poisct]

def Ctotal (C0, tr, vr, gr, ts=0, vs=0, gd=0):
    if vr == 0:
        print "vr  = 0"
#    (C0, tr, step, G1, ts, vs)
#    I4_ = Id4()    
#    I4 = Kelvin(I4_)
#    I4 = KelvintoVoigt(np.identity(6))
    I4 = np.identity(6)
    T1_ = T1(vr, tr)
    
    T2_ = T2(T1_, gd, vs, ts, I4)    
#    print 'tt\n', T1_,'\n',T2_
    invT1 = np.linalg.inv(T1_)
    temp1 = np.dot(invT1, T2_)
    temp2 = I4 + temp1
    invtemp2 = np.linalg.inv(temp2)
    temp3 = np.dot(T1_,invtemp2)
    Ct = C0 + temp3#*2.6
#    print temp3[0][0]
#    print 'tt\n', Ct[0][0], C0[0][0], temp3[0][0], T1_[0][0]
    return Ct



#def OPA (y, v, Cr, alpha, poisi):
#    Cdem = y
##    phi = v
#    I = np.identity(6)
#    Gdem = TensorG(y, alpha, poisi)
#    temp1 = 1./(1.-v)
#    temp2 = Cr-Cdem
#    temp3 = (I-Gdem*temp2)
#    temp4 = np.linalg.inv(temp3)
#    
#    dydv = [temp1*temp2*temp4]
#    return dydv
#
##def includeOPA (mineral1, mineral2, alpha, phimax, curphi, ni):  #ni = number of inclusions
##    
##    return
#
#def main4():
#    calcita = Calcite()
#    water = Water()
#    C0 = calcita[0] 
#    rhom = calcita[1]
#    poisi = calcita[2]
#    Cr = water[0] 
#    rhoi = water[1]
#    alpha = 0.3
#    phi = 0.4
##    rhom = 2.71
##    rhoi = 0.99
#    rho = (1.0 - phi)*rhom + phi*rhoi
#    print C0
#    y0 = C0
#    v = np.linspace(0,0.4,200)
#    from scipy.integrate import odeint
#    sol = odeint(OPA, y0, v, args=(Cr, alpha, poisi))
#    print sol
    
#def IncludeTest (mineral1, mineral2, alpha, phimax, curphi, ni): 
#    C0 = mineral1[0] 
#    rhom = mineral1[1]
#    poisi = mineral1[2]
#    Cr = mineral2[0]
#    rhoi = mineral2[1]
#    step = phimax/ni
#    phi_ = curphi
#    ts = np.zeros((6,6))
#    vs = np.zeros((6,6))
#    phi = np.empty(ni)
#    rho = np.empty(ni)
#    poison = np.empty(ni)
#    poison[0] = poisi
#    Ct = np.empty((ni,6,6))
#    for i in range (ni):
#        if i==0: 
#            poison[i] = poisi
#        else:
#            Gct = (C0[3][3]/2)
#            Kct = C0[0][0] - 4.*Gct/3        
#            poisi = (3.*Kct-2.*Gct)/(2.*(3.*Kct+Gct))
#            poison[i] = poisi
#        Gr = TensorGreen(alpha, C0)
#        dC = DeltaC(Cr, C0)
#        tr = T_r (dC , Gr)        
#        alphad = alpha*step
#        Gd = TensorGreen(alphad, C0)
#        Ct_ = Ctotal (C0, tr, step, Gr, ts, vs, Gd)
#        Ct[i] = Ct_
#        C0 = Ct_
#        ts = tr
#        vs = phi_
#        phi_ += step
#        phi[i] = phi_
##        
#        rho[i] = (1.0 - phi[i])*rhom + phi[i]*rhoi
#    return [Ct, rho, phi, poison]
def Ct (C0, T1_, T2_):
    I4 = np.identity(6)
    invT1 = np.linalg.inv(T1_)
    temp1 = np.dot(invT1, T2_)
    temp2 = I4 + temp1 #ou invT1
    invtemp2 = np.linalg.inv(temp2)
    temp3 = np.dot(T1_,invtemp2)
    Ct = C0 + temp3
#    print 't3\n',temp3
    return Ct

def CSingle (C0, tr, vr, gr):
    I4 = np.identity(6)
    T1_ = T1(vr, tr)
    T2_ = I4
#    T2_ = -np.dot(np.dot(T1_,gr),T1_)
    invT1 = np.linalg.inv(T1_)
    temp1 = np.dot(invT1, T2_)
    temp2 = I4 + temp1 #ou invT1
    invtemp2 = np.linalg.inv(temp2)
    temp3 = np.dot(T1_,invtemp2)
    Ct = C0 + temp3
#    print '\ncsingle', vr, C0, T1_
#    print temp3
    return Ct

def CFam (C0, tr, vr, gr, ts, vs, gdrs):
    I4 = np.identity(6)
    T1_ = T1(vr, tr)
    T2_ = T2(T1_, gdrs, vs, ts)
    invT1 = np.linalg.inv(T1_)
#    print 'tts', tr,'\n',gr,'\n',ts
    temp1 = np.dot(invT1, T2_)
    temp2 = I4 + temp1 #ou invT1
    invtemp2 = np.linalg.inv(temp2)
    temp3 = np.dot(T1_,invtemp2)
    Ct = C0 + temp3
#    print temp3
    return Ct

def IncludeFam (matrix, gr_fam, cavity, k, omega):  #ni = number of inclusions
#    gr_fam = [fam1, fam2]
#    fam# = [water, alpha, phi]
#    print 'rhoi',gr_fam[0]
#    print 'fam', omega
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    Cr_i = []
    rho_i = []
    alpha_i = []
    phi_i = []
    type_i=[]
    for fam in gr_fam:
#        print 'fam',fam
        Cr_i.append(fam[0][0])
        rho_i.append(fam[0][1])
        alpha_i.append(fam[1])
        phi_i.append(fam[2])
        type_i.append(fam[3])
#    print 'iii\n', rho_i, alpha_i, phi_i
#    for i in range (len(rho_i)):
#    print 'rhof', rho_i, phi_i, (1-np.sum(phi_i))*rhom+np.dot(rho_i,phi_i)
#        rhof = (1.-rho_i)*rhom + rho_i*phi_i
#    print "rhoi", rho_i, alpha_i
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]
#    vs=0.
#    ts = np.identity(6)
    poison = poisi
    C1=0.
    C2=0.
    vf = 0.
    vs=0.
    vr=0.
    tr=[]
    ts=[]
    vfin = 0.0
    cavity=0
#    alphad = 0.99
    for i in range(len(gr_fam)):
        vr = phi_i[i]
#        print vr,'vrrrrrr'
        alpha = alpha_i[i]
        Cr = Cr_i[i]
        rhoi = rho_i[i]
        
        Gr, er = TensorG(C0, alpha, poison) 
        dC = DeltaC(Cr, C0)
        if cavity == 0: 
            tr = T_r (dC , Gr) # nao considera frequencia
#            ts = T_r (dC, Gr)
        else: 
            tr = IncludeT (C0, Cr, Gr, poison, alpha, vr, k, omega) # considera frequencia
#        alphad = alpha/vr # alpha|min = alphad * vr
        alphad = alpha
        gdrs, ers = TensorG(C0, alphad, poison) 
        C1 = C1+T1(vr,tr)
        
        for j in range (len(gr_fam)):
            
            if j != i:
#            if j > i:
#                print 'iii', len(gr_fam),i, j
                vs = phi_i[j]
                Cs = Cr_i[j]
                dCs = DeltaC(Cs, C0)
                alphas = alpha_i[j]
                Gs, Es = TensorG(C0, alphas, poison) 
                ts = T_r (dCs , Gs)
#                print 'gdrs111\n',alpha, alphas, alpha_i
                C2 = C2+T2(vr, tr, gdrs, vs, ts)
#                C2=np.identity(6)
        vf = vf+vr
        if type_i[i]=='cav':
            vfin += vr
#        print 'vvff', vr, rhoi
    

#    C2 = np.identity(6)
#    print 'CC11\n', C1, '\n',vr*tr+vs*ts
    Ct_ = Ct (C0, C1, C2)
    rho_ = (1-np.sum(phi_i))*rhom+np.dot(rho_i,phi_i)
    lct.append(Ct_)
    lrho.append(rho_)
    
#    print 'vvrvrvrvr\n', Gr,'\n',Gs
#    for type_ in type_i: 
#        if type_=='cav': 
#            vfin = vfin+vr
#            print vfin,type_, phi_i, 'vfin'
#    if rtype == 'cav': vfin = vfin+vr
#    if stype == 'cav': vfin = vfin+vs
#    lphi.append(vr+vs)
    lphi.append(vfin)
#    lphi.append(np.sum(phi_i))
#    print 'vr', vf, lphi
#    vs = vr
#    ts = tr
    return [lct, lrho, lpoison], lphi

def IncludeSingle (matrix, sfam, cavity, omega=100):  #ni = number of inclusions
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    inclusion = sfam[0]
    alpha = sfam[1]
    phi_ = sfam[2]
    type_ = sfam[3]
    Cr = inclusion[0]
    rhoi = inclusion[1]
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]
    
    poison = poisi
    vr = phi_
    Gr,S0 = TensorG(C0, alpha, poison) 
    dC = DeltaC(Cr, C0)
    cavity = 1
    if cavity == 0: 
        tr = T_r (dC , Gr) # nao considera frequencia
#        print 'tr1', C0, Cr, dC
    else: 
        tr = IncludeT (C0, Cr, Gr, poison, alpha, vr, omega) #considera frequencia
#        print tr
#        print 'tr2', tr
    Ct_ = CSingle (C0, tr, vr, Gr)
    rho_ = (1.0 - vr)*rhom + vr*rhoi
    
    lct.append(Ct_)
    lrho.append(rho_)
    vfin = 0.0
    if type_ == 'cav': vfin = vfin+vr
    lphi.append(vfin)
#    lphi.append(vr)
#    print 'rhoo', len(lct), rho_, len(lphi), len(lpoison)
#    print Ct_
    poiss =poisi
    if type_ == 'agg': poiss = Poisson(Ct_[0][0]-4./3.*Ct_[3][3],Ct_[3][3])
    lpoison.append(poiss)
#    return [Ct_, rho_, poisi], vfin
    return [lct, lrho, lpoison], lphi

def IncludeDual (matrix, rfam, sfam, cavity, alphad):  #ni = number of inclusions
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    rinclusion = rfam[0]
    ralpha = rfam[1]
    rphi_ = rfam[2]
    rtype = rfam[3]
    sinclusion = sfam[0]
    salpha = sfam[1]
    sphi_ = sfam[2]
    stype = sfam[3]
    
    Cr = rinclusion[0]
    Cs = sinclusion[0]
    
    rhor = rinclusion[1]
    rhos = sinclusion[1]
    
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]
    
    poison = poisi
#    vr = rphi_
    vs = sphi_
    for vr in np.linspace(0,rphi_, 500):
        Gr,E = TensorG(C0, ralpha, poison)
        Gs,E = TensorG(C0, salpha, poison)
        
        dCr = DeltaC(Cr, C0)
        dCs = DeltaC(Cs, C0)
        cavity = 0
        if cavity == 0: 
            tr = T_r (dCr , Gr) # nao considera frequencia
            ts = T_r (dCs , Gs) # nao considera frequencia
    #        print 'tr1', tr
        else: 
            tr = IncludeT (C0, Cr, Gr, poison, ralpha, vr) #considera frequencia
            ts = IncludeT (C0, Cs, Gs, poison, salpha, vs) #considera frequencia
    #        print tr
    #    print 'tr2',cavity, ts
    #    alphad = 0.99
        gdrs,E = TensorG(C0, alphad, poison) 
        if vr==vs and tr.all() ==ts.all():
            C1 = T1(vr, tr)
            C2 = T2(vr, tr, gdrs, vr, tr)
            rho_ = (1.0 - vr)*rhom + vr*rhor
            print 'igual1'
        else: 
            C1 = T1(vr, tr)+T1(vs,ts)
            C2 = T2(vr, tr, gdrs, vs, ts)
            rho_ = (1.0 - vr - vs)*rhom + vr*rhor + vs*rhos
    #    print 'gdrs222\n',vr, tr, vs, ts
    #    print 'CC11', C1, C2
        Ct_ = Ct (C0, C1, C2)
    #    Ct_ = CFam (C0, tr, vr, Gr, ts, vs, gdrs)
    #    print Ct_
    #    rho_ = (1.0 - vr)*rhom + vr*rhor
        
    #    rho_ = (1.0 - vr - vs)*rhom + vr*rhor + vs*rhos
        lct.append(Ct_)
        lrho.append(rho_)
    #    lphi.append(vr)
    #    print 'innn\n\n', vr, vs
        vfin = 0.0
        if rtype == 'cav': vfin = vfin+vr
    #    if stype == 'cav': vfin = vfin+vs
        
        lphi.append(vr)
#        lphi.append(vfin)
#    print 'rhoo', len(lct), rho_, len(lphi), len(lpoison)
#    print Ct_
#    lpoison.append(poisi)
    return [lct, lrho, lpoison], lphi

def main3():
    import matplotlib.pyplot as plt
#    Vpm = 5.85 # km/s
#    Vsm = 3.9 # km/s
#    rhom = 2.650 # g/cm3
#    Vpf = 1.6 # km/s
#    rhof = 1.1 # g/cm3
#    phimax = 0.4    
#    Gm = Vsm*Vsm*rhom # GPa
#    Km = Vpm*Vpm*rhom - 4.0*Gm/3.0 # GPa
#    Kfl = Vpf*Vpf*rhof # GPa
#    print Km, Gm

    #inclusion
    quartzo = Quartzo()
    calcita = Calcite()
#    aragonita = Aragonite()
#    minmatrix = MineralMatriz()
    dolomita = Dolomite()
    clay = Clay()
    dry = Dry()
    water = Water()
    #cavity = 0(isolated)... 1 (exchange)
    cavity = 1
#    ni = 4
#    C0 = quartzo[0]
#    Ci = calcita[0]
#    alphai = 0.3
#    poisi = calcita[2]
#    i=0.
#    cont=0
#    composito = quartzo
#    while cont <10:
#        
#        fam10 = [clay, 0.05, 0.1, 'agg']
#        composito, phiii = IncludeSingle (composito, fam10, cavity)
#        Ctt = composito[0]
#        rhoo = composito[1]
#        poiscc = composito[2]
##        i = i+0.1
#        cont +=1
#        ni = len(phiii)
#        print "ii", cont#, '\n\nqtzo', quartzo, '\n\ncomp', composito, '\n\nsep',Ctt[ni-1], rhoo[ni-1], poiscc[ni-1]
#        composito = [Ctt[ni-1], rhoo[ni-1], poiscc[ni-1]]
#    vpp = np.empty(ni)
#    print 'teste', Ctt[-1], rhoo
#    vpp = (((Ctt[-1][0][0]/rhoo[-1]).real)**0.5)
#    print 'teste', vpp, rhoo[-1]

#    fam = [element, alpha, phi, 'type']
#    for j in range(len(macroporos)):#, 0.5]:
#    e=10
#    alpha=0.01
#    while e>0.01:
    vpo=[]
    omega = np.linspace(1,10000000,10000)
    for j in [1]:
#        Km= bulk[j]
#        Gm= shear[j]   
#        rhom= dens[j]
#        meso = macroaspect[j]
#        micro = microaspect[j]
#        phimeso = macroporos[j]
#        phimicro = microporos[j]
        print 'jjj', j
        phimax = 0.2303
        Km= 76.85
        Gm= 32.03    
        rhom= 2.71
        meso = 0.51
        micro = 0.04
#    porototal = 0.2303
        phimeso = 0.141
        phimicro = 0.0893
        minmatrix = MineralMatriz(Km, Gm, rhom)        
        fam1 = [water, 0.05, 0.04, 'cav']
#        fam1 = [dry, meso, phimeso, 'cav']
#        print meso, phimeso, micro, phimicro
#        fam1 = [dry, alpha, phimax, 'cav']
        fam2 = [dry, micro, phimicro, 'cav']
#        fam2 = [water, 0.1, 0.0893, 'cav']
#        fam3 = [dry, a, 0.05]
#        fam4 = [dry, a, 0.05]
#        fam5 = [dry, a, 0.05]
#        fam6 = [dry, a, 0.05]
#        gr_fam = [fam1, fam2]#, fam3, fam4, fam5]
#        fam1 = [dry, a, 0.25]
#        gr_fam = [fam1, fam2]
        gr_fam = [fam1,fam1]
#        for i in range(1): gr_fam.append(fam2)
#        gr_fam = [fam1, fam2]
    #    phi_fam = [0.1, 0.2]
    #    print water
    #    [Ct, rho, phi, poisct] = IncludeIso (calcita, water, alphai, porosf, porosi, ni)
    #    print 'len', len(gr_fam)
#        crackdensity = 3.*phi/(4.*np.pi*a)
        if len(gr_fam) == 1:
        #    for inc in gr_fam:
        #        [Ct, rho, phi, poisct] = IncludeSingle (calcita, water, alphai, porosf, cavity)
            composite, phi = IncludeSingle (calcita, gr_fam[0], cavity, j)
#            print 'n'
        else:
#            composite, phi = IncludeDual (quartzo, fam1, fam2, cavity)
            composite, phi = IncludeFam (minmatrix, gr_fam, cavity)
#        print 'result', phi,'\n', calcita,'\n',dry
        Ct = composite[0]
        rho = composite[1]
        
        poisct = composite[2]
    #    print VoigttoKelvin(Ct[-1])!
        ni = len(phi)
        Vp = np.empty(ni)
        Vs = np.empty(ni)
        Qp = np.empty(ni)
        ks= np.empty(ni)
        kk= np.empty(ni)
        ksg= np.empty(ni)
        g= np.empty(ni)
    #    print len(phi), len(ksg), len(Ct)
        for i in range(len(phi)):
            print phi[i], Ct[i][0][0]
    ##        print i
    ##        rho = (1.0 - phi[i])*rhom + phi[i]*rhoi
    #        
    #        print '\n',i, Ct[i][0][0], Ct[i][3][3]
            g[i] = Ct[i][3][3]/2
    #        g[i] = Ct[i][3][3]
            ks[i] = Ct[i][0][0] - 4./3.*g[i]
#            ksg[i] = Ct[i][0][0]
#            print Ct[i][0][0], phi, rho[i]    
            kk[i] = (1.0 - phi[i])*76.8 + phi[i]*2.2 + 4./3.*g[i]  
#            print 'vpppss', Ct[i][0][0]
            Vp[i] = (((Ct[i][0][0]/rho[i]).real)**0.5)
#        vpo.append(Vp[-1])
#        alpha= alpha+0.0001
#        e = np.abs(Vp[-1]-3.186)
#        if np.isnan(e): e=10
#        print e, alpha
#            Vs[i] = (((Ct[i][3][3]/(2.*rho[i])).real)**-0.5)**-1
    #        Qp[i] = ((Ct[i][0][0]).real)/((Ct[i][0][0]).img)
    #    print ks, rho, g, '\n', Ct[0], '\n\n', Ct[1]
        
    #    print rho
            
        
    #    print k
#        plt.subplot(211)
#        plt.plot(phi, ksg, 'r') # plot c11
#        plt.plot(phi, ks, 'r')
        plt.grid()
        plt.plot(phi, kk, 'y')
#        plt.plot(phi, g, 'r')
#    plt.plot(phi, rho, 'r')
#    plt.grid()
#    plt.show()
#        print 'ploto',Vp[-1], vpo
#        plt.plot(phi, Vp, 'r') # plot vp
#        plt.plot(phi, Vs, 'r') # plot vp
#    print 'main1 ultimo vp', Vp[-1], len(vpo), len(omega)
#    plt.show()
#    plt.semilogx(omega, vpo, 'r')
#    plt.xlim(1,10000000)
    plt.show()
        
    
def main4():
    import matplotlib.pyplot as plt
    #inclusion
    quartzo = Quartzo()
    calcita = Calcite()
#    aragonita = Aragonite()
#    minmatrix = MineralMatriz()
    dolomita = Dolomite()
    clay = Clay()
    dry = Dry()
    water = Water()
    #cavity = 0(isolated)... 1 (exchange)
    cavity = 1
#    ni = 4
#    C0 = quartzo[0]
#    Ci = calcita[0]
#    alphai = 0.3
#    poisi = calcita[2]
#    i=0.
#    cont=0
#    composito = quartzo
#    while cont <10:
#        
#        fam10 = [clay, 0.05, 0.1, 'agg']
#        composito, phiii = IncludeSingle (composito, fam10, cavity)
#        Ctt = composito[0]
#        rhoo = composito[1]
#        poiscc = composito[2]
##        i = i+0.1
#        cont +=1
#        ni = len(phiii)
#        print "ii", cont#, '\n\nqtzo', quartzo, '\n\ncomp', composito, '\n\nsep',Ctt[ni-1], rhoo[ni-1], poiscc[ni-1]
#        composito = [Ctt[ni-1], rhoo[ni-1], poiscc[ni-1]]
#    vpp = np.empty(ni)
#    print 'teste', Ctt[-1], rhoo
#    vpp = (((Ctt[-1][0][0]/rhoo[-1]).real)**0.5)
#    print 'teste', vpp, rhoo[-1]

#    fam = [element, alpha, phi, 'type']
#    for j in range(len(macroporos)):#, 0.5]:
#    e=10
#    alpha=0.01
#    while e>0.01:
    vpo=[]
    omega = np.linspace(1,10000000,10000)
    for j in [0.01, 0.1, 0.999]:
#    for j in [1]:
#        Km= bulk[j]
#        Gm= shear[j]   
#        rhom= dens[j]
#        meso = macroaspect[j]
#        micro = microaspect[j]
#        phimeso = macroporos[j]
#        phimicro = microporos[j]
#        print 'jjj', j
        phimax = 0.2303
        Km= 76.85
        Gm= 32.03    
        rhom= 2.71
        meso = 0.51
        micro = 0.04
#    porototal = 0.2303
        phimeso = 0.141
        phimicro = 0.0893
        alphad=1./0.99
        minmatrix = MineralMatriz(Km, Gm, rhom)        
        fam1 = [dry, 0.99, 0.99, 'cav'] # first order hu
#        fam1 = [dry, 0.99999, 0.99999, 'cav'] # second order hu
        print 'crack density', 3.*fam1[2]/(4.*np.pi*fam1[1])
        gr_fam = [fam1,fam1]
#        gr_fam = [fam1]
        if len(gr_fam) == 1:
            composite, phi = IncludeSingle (minmatrix, gr_fam[0], cavity, j)
            print 'n', "single"
        else:
            k = np.identity(3)*1000
            omega = 30
            fam1 = [dry, 1./j, 0.5, 'cav'] # second order hu
            fam2 = [dry, 1./j, 0.5, 'cav'] # second order hu
            composite, phi = DualCom (minmatrix, fam1, fam2, cavity, alphad, k, omega)
#            composite, phi = TMAFIRST (minmatrix, fam1, cavity)
#            composite, phi = TMASECOND (minmatrix, fam1, fam1, cavity, alphad)
#            composite, phi = IncludeDual (minmatrix, fam1, fam1, cavity, alphad)
            print '\n', "moree", 
        Ct = composite[0]
        rho = composite[1]
#        print len(Ct), len(rho), len(phi), '\n', Ct[-1], phi[-1], rho[-1]
        poisct = composite[2]
    #    print VoigttoKelvin(Ct[-1])!
        ni = len(phi)
        Vp = np.empty(ni)
        Vs = np.empty(ni)
        Qp = np.empty(ni)
        ks= np.empty(ni)
        kk= np.empty(ni)
        c11= np.empty(ni)
        c12= np.empty(ni)
        c22= np.empty(ni)
        c33= np.empty(ni)
        c44= np.empty(ni)
        c66= np.empty(ni)
        g= np.empty(ni)
    #    print len(phi), len(ksg), len(Ct)
        for i in range(len(phi)):
#            print len(phi), phi[i], Ct[i], '\n'
    ##        print i
    ##        rho = (1.0 - phi[i])*rhom + phi[i]*rhoi
    #        
    #        print '\n',i, Ct[i][0][0], Ct[i][3][3]
            g[i] = Ct[i][3][3]/2
    #        g[i] = Ct[i][3][3]
            ks[i] = Ct[i][0][0] - 4./3.*g[i]
            c11[i] = Ct[i][0][0] #c11
            c12[i] = Ct[i][0][1]
            c22[i] = Ct[i][1][1]
            c33[i] = Ct[i][2][2]
            c44[i] = Ct[i][3][3]/2.
            c66[i] = Ct[i][5][5]/2.
#            print Ct[i][0][0], phi, rho[i]    
            kk[i] = (1.0 - phi[i])*76.8 + phi[i]*2.2 + 4./3.*g[i]  
            
            Vp[i] = (((Ct[i][0][0]/rho[i]).real)**0.5)
#            print Ct[i][0][0]
#        vpo.append(Vp[-1])
#        alpha= alpha+0.0001
#        e = np.abs(Vp[-1]-3.186)
#        if np.isnan(e): e=10
#        print e, alpha
#            Vs[i] = (((Ct[i][3][3]/(2.*rho[i])).real)**-0.5)**-1
    #        Qp[i] = ((Ct[i][0][0]).real)/((Ct[i][0][0]).img)
    #    print ks, rho, g, '\n', Ct[0], '\n\n', Ct[1]
        
    #    print rho
            
        
    #    print k
#        plt.subplot(211)
#        print c44
        crackdens = np.empty(ni)
        for i in range(len(phi)): 
            if phi[i]==0:
                crackdens[i] = 00.
            else:
                crackdens[i] = (3.*phi[i])/(4.*np.pi*alphad)
#        print crackdens[0], ksg[0], minmatrix[0][0][0], len(ksg/minmatrix[0][0][0]), crackdens
#        plt.plot(crackdens, c11/c11[0], 'b') # plot c11
#        plt.plot(phi, c11/c11[0], 'r') # plot c11
        plt.plot(phi, Vp, 'r') # plot c11
#        plt.plot(phi, c12/c12[0], 'b') # plot c122
#        plt.plot(phi, c22/c22[0], '.') # plot c22
#        plt.plot(phi, c33/c33[0], '.') # plot c33
#        plt.plot(phi, c44/c44[0], 'b') # plot c44
#        plt.plot(phi, c66/c66[0], 'g') # plot c66
        plt.xlim(0,0.5)
#        plt.ylim(0,1)
#        plt.plot(phi, ks, 'r')
#        plt.grid()
#        plt.plot(phi, kk, 'y')
#        plt.plot(phi, g, 'r')
#    plt.plot(phi, rho, 'r')
#    plt.grid()
#    plt.show()
#        print 'ploto',Vp[-1], vpo
#        plt.plot(phi, Vp, 'r') # plot vp
#        plt.plot(phi, Vs, 'r') # plot vp
#    print 'main1 ultimo vp', Vp[-1], len(vpo), len(omega)
#    plt.show()
#    plt.semilogx(omega, vpo, 'r')
#    plt.xlim(1,10000000)
#    plt.show()
#    minmatrix, gr_fam, cavity
        
def main5():
    import matplotlib.pyplot as plt
    #inclusion
    quartzo = Quartzo()
    calcita = Calcite()
#    aragonita = Aragonite()
#    minmatrix = MineralMatriz()
    dolomita = Dolomite()
    clay = Clay()
    dry = Dry()
    water = Water()
    #cavity = 0(isolated)... 1 (exchange)
    cavity = 1
    vpo=[]
    omega = np.linspace(1,10000000,10000)
    for j in [2]:
#        print 'jjj', j
        phimax = 0.2303
        Km= 76.85
        Gm= 32.03    
        rhom= 2.71
        meso = 0.51
        micro = 0.04
#    porototal = 0.2303
        phimeso = 0.141
        phimicro = 0.0893
        alphad=1./0.99
        minmatrix = MineralMatriz(Km, Gm, rhom)   
#        minmatrix = calcita
#        fam1 = [dry, 0.05, 0.99, 'cav'] # first order hu
#        fam2 = [water, 0.99, 0.99, 'cav'] # second order hu
#        fam1 = [water, 0.99, 0.99, 'cav'] # second order hu
#        gr_fam = [fam1,fam2]
        gr_fam = [0,1]
#        print 'min',minmatrix[0][0][0]-4./3.*minmatrix[0][3][3], minmatrix[0][3][3]
        
        if len(gr_fam) == 1:
            composite, phi = IncludeSingle (minmatrix, gr_fam[0], cavity, j)
            
#            print 'n', "single"
        else:
#            if j==1:
#                k = np.ones((3,3))*10**5
##                k = np.identity(3)
#                fam2 = [water, 0.99, 0.49, 'cav']
#                omega = 10
#                composite, phi = SingleComPerm (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
            if j==1:
#                k = np.ones((3,3))*10**10
                k = np.identity(3)
#                fam2 = [water, 0.99, 0.69, 'cav']
#                omega = 10**8
#                composite, phi = SingleComFreq (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
#            if j==1:
##                k = np.ones((3,3))*100000000000000
#                k = np.identity(3)*0.0001
#                omega = 10
#                fam1 = [water, 0.99, 0.5, 'cav'] # second order hu
#                fam2 = [water, 0.99, 0.5, 'cav'] # second order hu
#                composite, phi = TMAFIRST (minmatrix, fam1, cavity)
#                composite, phi = DualCom (minmatrix, fam1, fam2, cavity, alphad)#, k, omega)
#                composite, phi = SingleCom (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
            else:
                k = np.identity(3)
                omega = 10
                fam1 = [water, 1./0.1, 0.7, 'cav'] # second order hu
                fam2 = [water, 1./0.1, 0.3, 'cav'] # second order hu
#                gr_fam = [fam1,fam2]
#                composite, phi = DualCom (minmatrix, fam1, fam2, cavity, alphad, k, omega)
#                composite, phi = IncludeFam (minmatrix, gr_fam, cavity, alphad, omega)#, k, omega)
#                composite, phi = TMAFIRST (minmatrix, fam1, cavity)
                composite, phi = TMASECOND (minmatrix, fam1, fam2, cavity, alphad)
#                composite, phi = SingleCom (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
            
#            print '\n', "moree"
        Ct = composite[0]
        rho = composite[1]
#        print '\n', (Ct[-1][0][0]/rho[-1])**0.5#, phi[-1], rho[-1]
        poisct = composite[2]
#        plt.plot (phi, rho)
    #    print VoigttoKelvin(Ct[-1])!
        ni = len(phi)
        Vp = np.empty(ni)
        Vs = np.empty(ni)
        Qp = np.empty(ni)
        ks= np.empty(ni)
        kk= np.empty(ni)
        ksg= np.empty(ni)
        c12= np.empty(ni)
        c22= np.empty(ni)
        c33= np.empty(ni)
        c44= np.empty(ni)
        c66= np.empty(ni)
        g= np.empty(ni)
    #    print len(phi), len(ksg), len(Ct)
        for i in range(len(phi)):
#            print len(phi), phi[i], Ct[i], '\n'
    ##        print i
    ##        rho = (1.0 - phi[i])*rhom + phi[i]*rhoi
    #        
    #        print '\n',i, Ct[i][0][0], Ct[i][3][3]
            g[i] = Ct[i][3][3].real/2
    #        g[i] = Ct[i][3][3]
            ks[i] = Ct[i][0][0].real- 4./3.*g[i]
            ksg[i] = Ct[i][0][0].real #c11
            c12[i] = Ct[i][0][1].real
            c22[i] = Ct[i][1][1].real
            c33[i] = Ct[i][2][2].real
            c44[i] = Ct[i][3][3].real
            c66[i] = Ct[i][5][5].real
#            print Ct[i][0][0], phi, rho[i]    
            kk[i] = (1.0 - phi[i])*76.8 + phi[i]*2.2 + 4./3.*g[i]  
#            print 'vpppss', Ct[i][0][0]
            Vp[i] = (((Ct[i][0][0]/rho[i]).real)**0.5)

        crackdens = np.empty(ni)
        for i in range(len(phi)): 
            if phi[i]==0:
                crackdens[i] = 00.
            else:
                crackdens[i] = (3.*phi[i])/(4.*np.pi*.99)
#        print crackdens[0], ksg[0], minmatrix[0][0][0], len(ksg/minmatrix[0][0][0]), crackdens
#        plt.plot(crackdens, ksg/minmatrix[0][0][0], 'g') # plot c11
        if j==1:
            plt.plot(phi[1:], ksg[1:]/minmatrix[0][0][0], 'b') # plot c11
            plt.xscale('symlog')
#        if j==1:
#            plt.plot(phi, ksg/minmatrix[0][0][0], 'b') # plot c11
#            plt.xscale('symlog')
#        if j==2: 
#            plt.plot(phi, ksg/minmatrix[0][0][0], 'b') # plot c11
        else:
            plt.plot(phi[1:], ksg[1:]/minmatrix[0][0][0], 'r') # plot c11
#        plt.plot(phi, c12/minmatrix[0][0][1], 'b') # plot c122
#        plt.plot(phi, c22, 'o') # plot c22
#            plt.plot(phi[1:], c33[1:]/minmatrix[0][2][2], 'g') # plot c33
#        plt.plot(phi[1:], c44[1:]/minmatrix[0][3][3], 'g') # plot c44
#        plt.plot(phi, c66/minmatrix[0][5][5], 'g') # plot c66
        plt.xlim(0,1)
        plt.ylim(0,1)
    plt.grid()
#    plt.show()
#    minmatrix, gr_fam, cavity


def main6():
    import matplotlib.pyplot as plt
    #inclusion
    quartzo = Quartzo()
    calcita = Calcite()
#    aragonita = Aragonite()
#    minmatrix = MineralMatriz()
    dolomita = Dolomite()
    clay = Clay()
    dry = Dry()
    water = Water()
    methane = Methane()
    #cavity = 0(isolated)... 1 (exchange)
    cavity = 1
    vpo=[]
    omega = np.linspace(1,10000000,10000)
    for j in [2]:
        print 'jjj', j
        phimax = 0.2303
        Km= 76.85
        Gm= 32.03    
        rhom= 2.71
        meso = 0.51
        micro = 0.04
#    porototal = 0.2303
        phimeso = 0.141
        phimicro = 0.0893
        alphad=1./0.1
        minmatrix = MineralMatriz(Km, Gm, rhom)
#        print minmatrix[0][0]        
#        fam = [dry, alpha, phi, 'cav']
#        fam1 = [dry, 0.05, 0.99, 'cav'] # first order hu
        fam2 = [water, 0.99, 0.99, 'cav'] # second order hu
        fam1 = [dry, 0.99, 0.99, 'cav'] # second order hu
        gr_fam = [fam1,fam2]
        
        
        if len(gr_fam) == 1:
            composite, phi = IncludeSingle (minmatrix, gr_fam[0], cavity, j)
            
#            print 'n', "single"
        else:
#            if j==1:
#                k = np.ones((3,3))*10**3
###                k = np.identity(3)
#                fam2 = [methane, 0.99, 0.49, 'cav']
#                omega = 100
#                print 'j=1'#abaixo phi = perm (retorna array de perm)
#                composite, phi = SingleComPerm (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
                
            if j==1:
#                k = np.ones((3,3))*10**10
                k = np.identity(3)
#                fam2 = [water, 0.09, 0.09, 'cav']
                fam1 = [water, 1./0.1, 0.3, 'cav'] 
                omega = 10**8
                composite, phi = SingleComFreq (minmatrix, fam1, cavity, alphad, k, omega)#single cavity comunicating
#            if j==1:
##                k = np.ones((3,3))*100000000000000
#                k = np.identity(3)*0.0001
#                omega = 10
#                fam1 = [water, 0.99, 0.5, 'cav'] # second order hu
#                fam2 = [water, 0.99, 0.5, 'cav'] # second order hu
#                composite, phi = TMAFIRST (minmatrix, fam1, cavity)
#                composite, phi = DualCom (minmatrix, fam1, fam2, cavity, alphad)#, k, omega)
#                composite, phi = SingleCom (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
            else:
                k = np.identity(3)
                omega = 100
                fam1 = [water, 1./0.1, 0.5, 'cav'] # second order hu
                fam2 = [water, 1./0.1, 0.5, 'cav'] # second order hu
#                gr_fam = [fam1,fam2]
#                composite, phi = DualCom (minmatrix, fam1, fam2, cavity, alphad, k, omega)
#                composite, phi = IncludeFam (minmatrix, gr_fam, cavity, alphad, omega)#, k, omega)
#                composite, phi = mainTMADem(minmatrix,fam1, fam2, cavity,alphad)
#                composite, phi = TMAFIRST (minmatrix, fam1, cavity)
                composite, phi = TMASECOND (minmatrix, fam1, fam2, cavity, alphad, k, omega)
#                composite, phi = SingleCom (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
            
#            print '\n', "moree"
        Ct = composite[0]
        rho = composite[1]
#        plt.plot (phi, rho)
#        print '\n', (Ct[-1][0][0]/rho[-1])**0.5#, phi[-1], rho[-1]
        poisct = composite[2]
    #    print VoigttoKelvin(Ct[-1])!
        ni = len(phi)
        Vp = np.empty(ni)
        Vs = np.empty(ni)
        Qp = np.empty(ni)
        ks= np.empty(ni)
        kk= np.empty(ni)
        ksg= np.empty(ni)
        c12= np.empty(ni)
        c22= np.empty(ni)
        c33= np.empty(ni)
        c44= np.empty(ni)
        c66= np.empty(ni)
        g= np.empty(ni)
    #    print len(phi), len(ksg), len(Ct)
        for i in range(len(phi)):
#            print len(phi), phi[i], Ct[i], '\n'
    ##        print i
    ##        rho = (1.0 - phi[i])*rhom + phi[i]*rhoi
    #        
    #        print '\n',i, Ct[i][0][0], Ct[i][3][3]
            g[i] = Ct[i][3][3]/2
    #        g[i] = Ct[i][3][3]
            ks[i] = Ct[i][0][0] - 4./3.*g[i]
            ksg[i] = Ct[i][0][0] #c11
            c12[i] = Ct[i][0][1]
            c22[i] = Ct[i][1][1]
            c33[i] = Ct[i][2][2]
            c44[i] = Ct[i][3][3]
            c66[i] = Ct[i][5][5]
#            print Ct[i][0][0], phi, rho[i]    
            kk[i] = (1.0 - phi[i])*76.8 + phi[i]*2.2 + 4./3.*g[i]  
#            print 'vpppss', Ct[i][0][0]
            Vp[i] = (((Ct[i][0][0]/rho[i]).real)**0.5)

        crackdens = np.empty(ni)
        for i in range(len(phi)): 
            if phi[i]==0:
                crackdens[i] = 00.
            else:
                crackdens[i] = (3.*phi[i])/(4.*np.pi*.99)
#        print crackdens[0], ksg[0], minmatrix[0][0][0], len(ksg/minmatrix[0][0][0]), crackdens
#        plt.plot(crackdens, ksg/minmatrix[0][0][0], 'g') # plot c11
        if j==1:
            plt.plot(phi[1:], ksg[1:], 'b') # plot c11
            plt.xscale('symlog')
#        if j==1:
#            plt.plot(phi, ksg/minmatrix[0][0][0], 'b') # plot c11
#            plt.xscale('symlog')
#        if j==2: 
#            plt.plot(phi, ksg/minmatrix[0][0][0], 'b') # plot c11
        else:
#            plt.plot(phi, g, 'b') # plot c11
#            plt.plot(phi, ks, 'b') # plot c11
            plt.plot(phi, ks/ks[0], 'b') # plot c11
#            plt.plot(phi[1:], Vp[1:], 'b') # plot c11
#            plt.plot(phi[1:], rho[1:], 'b') # plot c11
#        plt.plot(phi, c12/minmatrix[0][0][1], 'b') # plot c122
#        plt.plot(phi, c22, 'o') # plot c22
#        plt.plot(phi, c33/minmatrix[0][2][2], 'b') # plot c33
#        plt.plot(phi, c44/minmatrix[0][3][3], 'b') # plot c44
#        plt.plot(phi, c66/minmatrix[0][5][5], 'g') # plot c66
#        plt.xlim(0,0.5)
#        plt.ylim(0,80)
    plt.grid(True)
#    plt.show()
#    minmatrix, gr_fam, cavity
    
def main7():
    import matplotlib.pyplot as plt
    #inclusion
    quartzo = Quartzo()
    calcita = Calcite()
#    aragonita = Aragonite()
#    minmatrix = MineralMatriz()
    dolomita = Dolomite()
    clay = Clay()
    dry = Dry()
    water = Water()
    methane = Methane()
    #cavity = 0(isolated)... 1 (exchange)
    cavity = 1
    vpo=[]
    omega = np.linspace(1,10000000,10000)
    for j in [2]:
        print 'jjj', j
        phimax = 0.2303
        Km= 76.85
        Gm= 32.03    
        rhom= 2.71
        meso = 0.51
        micro = 0.04
#    porototal = 0.2303
        phimeso = 0.141
        phimicro = 0.0893
        alphad=1./0.1
        minmatrix = MineralMatriz(Km, Gm, rhom)
#        print minmatrix[0][0]        
#        fam1 = [dry, 0.05, 0.99, 'cav'] # first order hu
        fam2 = [water, 0.99, 0.99, 'cav'] # second order hu
        fam1 = [dry, 0.99, 0.99, 'cav'] # second order hu
        gr_fam = [fam1,fam2]
        
        
        if len(gr_fam) == 1:
            composite, phi = IncludeSingle (minmatrix, gr_fam[0], cavity, j)
            
#            print 'n', "single"
        else:
#            if j==1:
#                k = np.ones((3,3))*10**3
###                k = np.identity(3)
#                fam2 = [methane, 0.99, 0.49, 'cav']
#                omega = 100
#                print 'j=1'#abaixo phi = perm (retorna array de perm)
#                composite, phi = SingleComPerm (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
                
            if j==1:
                k = np.ones((3,3))*10**10
#                k = np.identity(3)
                fam2 = [water, 0.09, 0.09, 'cav']
                omega = 10**10
                composite, phi = SingleComFreq (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
#            if j==1:
##                k = np.ones((3,3))*100000000000000
#                k = np.identity(3)*0.0001
#                omega = 10
#                fam1 = [water, 0.99, 0.5, 'cav'] # second order hu
#                fam2 = [water, 0.99, 0.5, 'cav'] # second order hu
#                composite, phi = TMAFIRST (minmatrix, fam1, cavity)
#                composite, phi = DualCom (minmatrix, fam1, fam2, cavity, alphad)#, k, omega)
#                composite, phi = SingleCom (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
            else:
                k = np.identity(3)*100000
                omega = 100
                fam1 = [water, 1./0.1, 0.7, 'cav'] # second order hu
                fam2 = [water, 1./0.1, 0.3, 'cav'] # second order hu
#                gr_fam = [fam1,fam2]
                composite, phi = DualCom (minmatrix, fam1, fam2, cavity, alphad, k, omega)
#                composite, phi = IncludeFam (minmatrix, gr_fam, cavity, alphad, omega)#, k, omega)
#                composite, phi = mainTMADem(minmatrix,fam1, fam2, cavity,alphad)
#                composite, phi = TMAFIRST (minmatrix, fam1, cavity)
#                composite, phi = TMASECOND (minmatrix, fam1, fam2, cavity, alphad, k , omega)
#                composite, phi = SingleCom (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
            
#            print '\n', "moree"
        Ct = composite[0]
        rho = composite[1]
#        plt.plot (phi, rho)
#        print '\n', (Ct[-1][0][0]/rho[-1])**0.5#, phi[-1], rho[-1]
        poisct = composite[2]
    #    print VoigttoKelvin(Ct[-1])!
        ni = len(phi)
        Vp = np.empty(ni)
        Vs = np.empty(ni)
        Qp = np.empty(ni)
        ks= np.empty(ni)
        kk= np.empty(ni)
        ksg= np.empty(ni)
        c12= np.empty(ni)
        c22= np.empty(ni)
        c33= np.empty(ni)
        c44= np.empty(ni)
        c66= np.empty(ni)
        g= np.empty(ni)
    #    print len(phi), len(ksg), len(Ct)
        for i in range(len(phi)):
#            print len(phi), phi[i], Ct[i], '\n'
    ##        print i
    ##        rho = (1.0 - phi[i])*rhom + phi[i]*rhoi
    #        
    #        print '\n',i, Ct[i][0][0], Ct[i][3][3]
            g[i] = Ct[i][3][3]/2
    #        g[i] = Ct[i][3][3]
            ks[i] = Ct[i][0][0] - 4./3.*g[i]
            ksg[i] = Ct[i][0][0] #c11
            c12[i] = Ct[i][0][1]
            c22[i] = Ct[i][1][1]
            c33[i] = Ct[i][2][2]
            c44[i] = Ct[i][3][3]
            c66[i] = Ct[i][5][5]
#            print Ct[i][0][0], phi, rho[i]    
            kk[i] = (1.0 - phi[i])*76.8 + phi[i]*2.2 + 4./3.*g[i]  
#            print 'vpppss', Ct[i][0][0]
            Vp[i] = (((Ct[i][0][0]/rho[i]).real)**0.5)

        crackdens = np.empty(ni)
        for i in range(len(phi)): 
            if phi[i]==0:
                crackdens[i] = 00.
            else:
                crackdens[i] = (3.*phi[i])/(4.*np.pi*.99)
#        print crackdens[0], ksg[0], minmatrix[0][0][0], len(ksg/minmatrix[0][0][0]), crackdens
#        plt.plot(crackdens, ksg/minmatrix[0][0][0], 'g') # plot c11
        if j==1:
            plt.plot(phi[1:], ksg[1:]/minmatrix[0][0][0], 'b') # plot c11
            plt.xscale('symlog')
#        if j==1:
#            plt.plot(phi, ksg/minmatrix[0][0][0], 'b') # plot c11
#            plt.xscale('symlog')
#        if j==2: 
#            plt.plot(phi, ksg/minmatrix[0][0][0], 'b') # plot c11
        else:
#            plt.plot(phi, g, 'b') # plot c11
#            plt.plot(phi, ks, 'r') # plot c11
            plt.plot(phi, ks/ks[0], 'r') # plot c11
#            plt.plot(phi[1:], Vp[1:], 'b') # plot c11
#            plt.plot(phi[1:], rho[1:], 'b') # plot c11
#        plt.plot(phi, c12/minmatrix[0][0][1], 'b') # plot c122
#        plt.plot(phi, c22, 'o') # plot c22
#        plt.plot(phi, c33/minmatrix[0][2][2], 'b') # plot c33
#        plt.plot(phi, c44/minmatrix[0][3][3], 'b') # plot c44
#        plt.plot(phi, c66/minmatrix[0][5][5], 'g') # plot c66
#        plt.xlim(0,0.5)
#        plt.ylim(0,80)
    plt.grid(True)
#    plt.show()
#    minmatrix, gr_fam, cavity
    

def tmatrix_log(rock, fluid, phit, rhob, perm, inclusion1, inclusion2, alphadual, method):
#    rock=(km,gm,rhom)
#    fluid= (kf,rhof)
#    phit= log
#    rhob = log
#    perm = core
#    method = {first, higher}
#    inclusion = (fraction, alpha)
    
    omega = 4000 #well log    
    minmatrix = MineralMatriz(rock[0], rock[1], rock[2]) 
    fluidsat = FluidSat(fluid[0], fluid[1])
    frac1 = inclusion1[0]
    alpha1 = inclusion1[1]
    poros1 = frac1*phit
    frac2 = inclusion2[0]
    alpha2 = inclusion2[1]
    poros2 = frac2*phit
    fam1 = (fluidsat, 1./alpha1, poros1, 'cav')
    fam2 = (fluidsat, 1./alpha2, poros2, 'cav')
    alphad = 1./alphadual
    import matplotlib.pyplot as plt
    #inclusion
    #cavity = 0(isolated)... 1 (exchange)
    cavity = 1
    for j in [2]:
        print 'jjj', j
        Km= 76.85
        Gm= 32.03    
        rhom= 2.71

               
        if method == 'first':
            composite, phi = IncludeSingle (minmatrix, fam1, cavity, omega)
        else:
#            fam1 = [water, 1./0.1, 0.7, 'cav'] # second order hu
#            fam2 = [water, 1./0.1, 0.3, 'cav'] # second order hu
#                gr_fam = [fam1,fam2]
            composite, phi = DualCom (minmatrix, fam1, fam2, cavity, alphad, perm, omega)
#                composite, phi = IncludeFam (minmatrix, gr_fam, cavity, alphad, omega)#, k, omega)
#                composite, phi = TMAFIRST (minmatrix, fam1, cavity)
#                composite, phi = TMASECOND (minmatrix, fam1, fam2, cavity, alphad, k , omega)
        
#            print '\n', "moree"
        Ct = composite[0]
        rho = composite[1]
#        plt.plot (phi, rho)
#        print '\n', (Ct[-1][0][0]/rho[-1])**0.5#, phi[-1], rho[-1]
        poisct = composite[2]
    #    print VoigttoKelvin(Ct[-1])!
        ni = len(phi)
        Vp = np.empty(ni)
        Vs = np.empty(ni)
        Qp = np.empty(ni)
        ks= np.empty(ni)
        kk= np.empty(ni)
        ksg= np.empty(ni)
        c12= np.empty(ni)
        c22= np.empty(ni)
        c33= np.empty(ni)
        c44= np.empty(ni)
        c66= np.empty(ni)
        g= np.empty(ni)
    #    print len(phi), len(ksg), len(Ct)
        for i in range(len(phi)):
#            print len(phi), phi[i], Ct[i], '\n'
    ##        print i
    ##        rho = (1.0 - phi[i])*rhom + phi[i]*rhoi
    #        
    #        print '\n',i, Ct[i][0][0], Ct[i][3][3]
            g[i] = Ct[i][3][3]/2
    #        g[i] = Ct[i][3][3]
            ks[i] = Ct[i][0][0] - 4./3.*g[i]
            ksg[i] = Ct[i][0][0] #c11
            c12[i] = Ct[i][0][1]
            c22[i] = Ct[i][1][1]
            c33[i] = Ct[i][2][2]
            c44[i] = Ct[i][3][3]
            c66[i] = Ct[i][5][5]
#            print Ct[i][0][0], phi, rho[i]    
            kk[i] = (1.0 - phi[i])*76.8 + phi[i]*2.2 + 4./3.*g[i]  
#            print 'vpppss', Ct[i][0][0]
            Vp[i] = (((Ct[i][0][0]/rho[i]).real)**0.5)

        crackdens = np.empty(ni)
        for i in range(len(phi)): 
            if phi[i]==0:
                crackdens[i] = 00.
            else:
                crackdens[i] = (3.*phi[i])/(4.*np.pi*.99)
#        print crackdens[0], ksg[0], minmatrix[0][0][0], len(ksg/minmatrix[0][0][0]), crackdens
#        plt.plot(crackdens, ksg/minmatrix[0][0][0], 'g') # plot c11
        if j==1:
            plt.plot(phi[1:], ksg[1:]/minmatrix[0][0][0], 'b') # plot c11
            plt.xscale('symlog')
#        if j==1:
#            plt.plot(phi, ksg/minmatrix[0][0][0], 'b') # plot c11
#            plt.xscale('symlog')
#        if j==2: 
#            plt.plot(phi, ksg/minmatrix[0][0][0], 'b') # plot c11
        else:
#            plt.plot(phi, g, 'b') # plot c11
#            plt.plot(phi, ks, 'r') # plot c11
            plt.plot(phi, ks/ks[0], 'r') # plot c11
#            plt.plot(phi[1:], Vp[1:], 'b') # plot c11
#            plt.plot(phi[1:], rho[1:], 'b') # plot c11
#        plt.plot(phi, c12/minmatrix[0][0][1], 'b') # plot c122
#        plt.plot(phi, c22, 'o') # plot c22
#        plt.plot(phi, c33/minmatrix[0][2][2], 'b') # plot c33
#        plt.plot(phi, c44/minmatrix[0][3][3], 'b') # plot c44
#        plt.plot(phi, c66/minmatrix[0][5][5], 'g') # plot c66
#        plt.xlim(0,0.5)
#        plt.ylim(0,80)
    plt.grid(True)
#    plt.show()
#    minmatrix, gr_fam, cavity

    
def DualCom (matrix, rfam, sfam, cavity,alphad, k=None, omega=None):  #ni = number of inclusions
    C0 =  matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    Rinclusion = rfam[0]
    Ralpha = rfam[1]
    Rphi_ = rfam[2]
    Sinclusion = sfam[0]
    Salpha = sfam[1]
    Sphi_ = sfam[2]
    Cr = Rinclusion[0]
    rhor = Rinclusion[1]
    viscr = Rinclusion[3]
    Cs = Sinclusion[0]
    rhos = Sinclusion[1]
    viscs = Sinclusion[3]
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]
#    Rphi_ = 0.5
#    Sphi_ = 0.5
#    Ralpha = 0.99
#    Salpha = 0.01
    poison = poisi
    Gr, E1 = TensorG(C0, Ralpha, poison) 
    Gs, Es = TensorG(C0, Salpha, poison) 
#    print 'grr\n', Gr, Gs
    dCr = Cr-C0
    dCs = Cs-C0
    tr = T_r (dCr , Gr)
    ts = T_r (dCs , Gs)    
    Gd, E2 = TensorG (C0, alphad, poison)
    I6 = np.identity(6)
    ni =500
#    print k, omega
#    if omega:
#        arqtestes = open("print.txt",'w')
#    else: 
#        arqtestes = open("printdry.txt",'w')
    for p in range(ni)[1:]:
        vr = p*Rphi_/ni
        vs = p*Sphi_/ni
        #novoniterp
#        alphad = Ralpha
#        Gd, E2 = TensorG (C0, alphad, poison)
        rho_ = (1.0 - vr - vs)*rhom + vr*rhor + vs*rhos
        if omega:
            tr = IncludeT (C0, Cr, Gr, poison, Ralpha, Rphi_, k, omega, rhom, viscr)
            ts = IncludeT (C0, Cs, Gs, poison, Salpha, Sphi_, k, omega, rhom, viscs)
        temp1 = vr*tr
        temp2 = vs*ts
        temp3 = temp1 + temp2
        temp4 = I6 + np.dot(Gd, temp3)
#        tempa = np.dot(temp1,Gd)
#        tempa2 = np.dot(temp2,Gd)
#        tempb = np.dot(tempa,temp2)
#        tempb2 = np.dot(tempa2,temp1)
#        tempd = tempb + tempb2
##        print 'temp3', temp3
#        tempc = np.linalg.inv(temp3)
#        
#        temp4 = I6 + np.dot(tempc,tempd)
        temp5 = np.linalg.inv(temp4)
        temp6 = np.dot(temp3, temp5)
#        print temp1[0][0], temp2[0][0],temp4[0][0]
#        arqtestes.write(str(str(tr[0][0]) +'  '+ str(ts[0][0]) +'  '+ str(temp3[0][0])+ '\n'))
        Ct = C0 + temp6
#        rho_ = (1.0 - vr - vs)*rhom + vr*rhor + vs*rhos
        lct.append(Ct)
        lphi.append((vr+vs))
        lrho.append(rho_)
#        C0 = Ct
#        print Ct[0][0]-2.*Ct[3][3]/3, Ct[3][3]/2.
#        poison = Poisson(Ct[0][0]-2.*Ct[3][3]/3, Ct[3][3]/2.).real
#        print poison.real
#    print 'lenn', temp6
    return [lct, lrho, lpoison], lphi

def SingleCom (matrix, rfam, cavity,alphad, k, omega):  #ni = number of inclusions
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    Rinclusion = rfam[0]
    Ralpha = rfam[1]
    Rphi_ = rfam[2]
    Cr = Rinclusion[0]
    rhor = Rinclusion[1]
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]    
    poison = poisi
#    print 'cr', Cr
    Gr, E1 = TensorG(C0, Ralpha, poison) 
#    Gr2 = TensorGreen(Ralpha, C0)
#    print 'grrsss\n', Gr1, '\n\n',Gr2
    dCr = Cr-C0
    tr = T_r (dCr , Gr)
#    tr = IncludeT (C0, Cr, Gr, poison, Ralpha, Rphi_)
    Gd, E2 = TensorG (C0, alphad, poison)
    I6 = np.identity(6)
#    tempE1 = I6 - E1
#    tempE2 = I6 - E2
#    invtempE1 = np.linalg.inv(tempE1)  
    ts=tr
    ni =500
    for vr in np.linspace(0,Rphi_,ni):       
        tr = IncludeT (C0, Cr, Gr, poison, Ralpha, Rphi_, k, omega)
        
#        temp4 = np.dot(E2, invtempEs)
#        temp3 = I6 + vs*temp4
#        invtemp3 = np.linalg.inv(temp3)
#        temp2 = np.dot(invtempE1, invtemp3)
#        temp1 = I6 - vr*temp2
#        Ct = np.dot(C0,temp1)
        vs=vr
        temp4 = ts*vs
        temp3 = I6 + np.dot(Gd,temp4)
        temp2 = np.linalg.inv(temp3)
        temp1 = vr*tr
        Ct = C0 + np.dot(temp1, temp2)  
        rho_ = (1.0 - vr)*rhom + vr*rhor
        lct.append(Ct)
        lphi.append(vr)
        lrho.append(rho_)
    return [lct, lrho, lpoison], lphi
        
def SingleComPerm (matrix, rfam, cavity,alphad, k, omega):  #ni = number of inclusions
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    Rinclusion = rfam[0]
    Ralpha = rfam[1]
    Rphi_ = rfam[2]
    Cr = Rinclusion[0]
    rhor = Rinclusion[1]
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]    
    poison = poisi
#    print 'cr', Cr
    Gr, E1 = TensorG(C0, Ralpha, poison) 
#    Gr2 = TensorGreen(Ralpha, C0)
#    print 'grrsss\n', Gr1, '\n\n',Gr2
    dCr = Cr-C0
    tr = T_r (dCr , Gr)
#    tr = IncludeT (C0, Cr, Gr, poison, Ralpha, Rphi_)
    Gd, E2 = TensorG (C0, alphad, poison)
    I6 = np.identity(6)
#    tempE1 = I6 - E1
#    tempE2 = I6 - E2
#    invtempE1 = np.linalg.inv(tempE1)  
    ts=tr
    ni =100
    vr=Rphi_
#    print k[0][0]
    for perm in np.linspace(0,k[0][0],ni):       
        tr = IncludeT (C0, Cr, Gr, poison, Ralpha, Rphi_, perm*np.ones((3,3)), omega)
        
#        temp4 = np.dot(E2, invtempEs)
#        temp3 = I6 + vs*temp4
#        invtemp3 = np.linalg.inv(temp3)
#        temp2 = np.dot(invtempE1, invtemp3)
#        temp1 = I6 - vr*temp2
#        Ct = np.dot(C0,temp1)
        vs=vr
        temp4 = ts*vs
        temp3 = I6 + np.dot(Gd,temp4)
        temp2 = np.linalg.inv(temp3)
        temp1 = vr*tr
        Ct = C0 + np.dot(temp1, temp2)  
        rho_ = (1.0 - vr)*rhom + vr*rhor
        print perm, Ct[0][0]
        lct.append(Ct)
        lphi.append(perm)
        lrho.append(rho_)
    return [lct, lrho, lpoison], lphi

def SingleComFreq (matrix, rfam, cavity,alphad, k, omega):  #ni = number of inclusions
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    Rinclusion = rfam[0]
    Ralpha = rfam[1]
    Rphi_ = rfam[2]
    Cr = Rinclusion[0]
    rhor = Rinclusion[1]
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]    
    poison = poisi
#    print 'cr', Cr
    Gr, E1 = TensorG(C0, Ralpha, poison) 
#    Gr2 = TensorGreen(Ralpha, C0)
#    print 'grrsss\n', Gr1, '\n\n',Gr2
    dCr = Cr-C0
    tr = T_r (dCr , Gr)
#    tr = IncludeT (C0, Cr, Gr, poison, Ralpha, Rphi_)
    Gd, E2 = TensorG (C0, alphad, poison)
    I6 = np.identity(6)
#    tempE1 = I6 - E1
#    tempE2 = I6 - E2
#    invtempE1 = np.linalg.inv(tempE1)  
    ts=tr
    ni =100
    vr=Rphi_
    potomega = np.log10(omega)
    for freq in np.logspace(0,potomega,ni):       
        tr = IncludeT (C0, Cr, Gr, poison, Ralpha, Rphi_, k, freq)
        
#        temp4 = np.dot(E2, invtempEs)
#        temp3 = I6 + vs*temp4
#        invtemp3 = np.linalg.inv(temp3)
#        temp2 = np.dot(invtempE1, invtemp3)
#        temp1 = I6 - vr*temp2
#        Ct = np.dot(C0,temp1)
        vs=vr
        temp4 = ts*vs
        temp3 = I6 + np.dot(Gd,temp4)
        temp2 = np.linalg.inv(temp3)
        temp1 = vr*tr
        Ct = C0 + np.dot(temp1, temp2)  
        rho_ = (1.0 - vr)*rhom + vr*rhor
        print freq, Ct[0][0]
        lct.append(Ct)
        lphi.append(freq)
        lrho.append(rho_)
    return [lct, lrho, lpoison], lphi

def TMAFIRST (matrix, sfam, cavity, omega=100):  #ni = number of inclusions
    
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    inclusion = sfam[0]
    alpha = sfam[1]
    phi_ = sfam[2]
#    Cr = inclusion[0]
    rhoi = 0.
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    
#    C0 = KelvintoVoigt(matrix[0])
#    Cr = KelvintoVoigt(inclusion[0])
    
    lct = [C0]
    
    
    
    poison = poisi
    vr = phi_
    Gr, E = TensorG(C0, alpha, poison) 
#    Gr = TensorGreen (alpha, C0)
#    E = -np.dot(Gr2,VoigttoKelvin(np.linalg.inv(KelvintoVoigt(C0))))
#    dC = DeltaC(Cr, C0)
#    cavity = 1
#    dC = VoigttoKelvin(DeltaC(Cr, C0))
    I6 = np.identity(6)
#    Qr = I6+Gr*(Cr - Ct)
#    Qs = I6+Gs*(Cs - Ct)
#    Ct = vr*Cr*Qr*(vs*Qs)**-1
#    Km = 76.85
#    Gm = 32.03
#    a = (3.*Km+Gm)/(9.*Km*Gm)
#    b = 1./Gm
#    c = (2*Gm-3*Km)/(18*Km*Gm)
#    print 'tfirst\n',KelvintoVoigt(C0), '\n\n', np.linalg.inv(C0), '\n',VoigttoKelvin(np.linalg.inv(KelvintoVoigt(C0))),'\n', a,b,c, '\n', E
    temp3 = I6 - E
    temp2 = np.linalg.inv(temp3)
    for vr_ in np.linspace(0,vr,500):
        temp1 = I6-vr_*temp2
        Ct = np.dot(C0,temp1)    
#        print rhom, rhoi
        rho_ = (1.0 - vr_)*rhom + vr_*rhoi
        
        lct.append(Ct)
        lphi.append(vr_)
        lrho.append(rho_)
    return [lct, lrho, lpoison], lphi

def TMASECOND (matrix, rfam, sfam, cavity,alphad, k, omega):  #ni = number of inclusions
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    Rinclusion = rfam[0]
    Ralpha = rfam[1]
    Rphi_ = rfam[2]
    Sinclusion = sfam[0]
    Salpha = sfam[1]
    Sphi_ = sfam[2]
    Cr = Rinclusion[0]
    rhor = Rinclusion[1]
    viscr = Rinclusion[3]
    Cs = Sinclusion[0]
    rhos = Sinclusion[1]
    viscs = Sinclusion[3]
#    alphad = 0.99
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]
#    Rphi_ = 0.5
#    Sphi_ = 0.5
#    Ralpha = 0.99
#    Salpha = 0.01
    poison = poisi
#    Salpha = Ralpha
    Gr, E1 = TensorG(C0, Ralpha, poison) 
#    Gr = TensorGreen (Ralpha, C0)
#    print 'grss\n', Gr, '\n\n',Gr2, C0
    Gs, Es = TensorG(C0, Salpha, poison) 
#    dCr = DeltaC(Cr, C0)
    dCr = Cr-C0
    dCs = Cs-C0
    tr = T_r (dCr , Gr)
    ts = T_r (dCs , Gs)    
    Gd, E2 = TensorG (C0, alphad, poison)
#    Gd = TensorGreen (alphad, C0)
#    I6 = SVoigttoKelvin(np.identity(6))
    I6 = np.identity(6)
#    vs=vr
#    tempE1 = I6 - E1
#    tempE2 = I6 - E2
#    tempEs = I6 - Es
#    invtempE1 = np.linalg.inv(tempE1)
#    invtempE2 = np.linalg.inv(tempE2)
#    invtempEs = np.linalg.inv(tempEs)    
#    ts=tr
    
#    Rphi_ = Rphi_*0.5
#    Sphi_ = Rphi_
#    Rphi_ = Rphi_+0.1
    ni =500
#    for vr in np.linspace(0,Rphi_,ni):
#        vs=vr
#        for vs in np.linspace(0,Sphi_,ni):
#        Gr, E = TensorG(C0, Ralpha, poison)
#        dCr = DeltaC(Cr, C0)
#        tr = T_r (dCr , Gr)
#        alphad = Ralpha/vr
#        print Ralpha/vr, alphad
#        Gd = TensorGreen (alphad, C0)
        
#        Gd, E = TensorG (C0, alphad, poison)
#        Gd = TensorGreen (alphad, C0)
#    for p in range(ni):    
#        vr = p*Rphi_/ni
#        vs = p*Sphi_/ni
#        temp4 = np.dot(E2, invtempEs)
#        temp3 = I6 + vs*temp4
#        invtemp3 = np.linalg.inv(temp3)
#        temp2 = np.dot(invtempE1, invtemp3)
#        temp1 = I6 - vr*temp2
#        Ct = np.dot(C0,temp1)
    
#        vs=vr
#        temp4 = ts*vs
#        temp3 = I6 + np.dot(Gd,temp4)
#        temp2 = np.linalg.inv(temp3)
#        temp1 = vr*tr
#        Ct = C0 + np.dot(temp1, temp2)*2.  
#    print C0[0][0]-2.*C0[3][3]/3, C0[3][3]/2.
    for p in range(ni):
        if omega:
#            print 'tem k e omega'
            tr = IncludeT (C0, Cr, Gr, poison, Ralpha, Rphi_, k, omega, rhom, viscr)
            ts = IncludeT (C0, Cs, Gs, poison, Salpha, Sphi_, k, omega, rhom, viscs)
        vr = p*Rphi_/ni
        vs = p*Sphi_/ni
#        alphad = Ralpha
#        Gd, E2 = TensorG (C0, alphad, poison)
#        temp4 = ts*vs
#        temp3 = I6 + np.dot(Gd,temp4)
#        temp2 = np.linalg.inv(temp3)
#        temp1 = vr*tr
#        
#        temp5 = I6 + np.dot(Gd,temp1)
#        temp6 = np.linalg.inv(temp5)
#        Ct = C0 + np.dot(temp1, temp2) + np.dot(temp4,temp6)
        #novoniterp
#        tr = T_r (dCr , Gr)
#        ts = T_r (dCs , Gs) 
        temp1 = vr*tr
        temp2 = vs*ts
        temp3 = temp1 + temp2
        temp4 = I6 + np.dot(Gd, temp3)
        temp5 = np.linalg.inv(temp4)
        temp6 = np.dot(temp3, temp5)
#        print temp1[0][0], temp2[0][0], temp3[0][0]
        Ct = C0 + temp6
        
        C1 = np.dot(vr,tr)
        tempg = I6 + np.dot(Gd,C1)
        tempginv = np.linalg.inv(tempg)
        tempgg = np.dot(C1,tempginv)
#        Ct = C0+tempgg
#        
#        print 'trr', tr, np.dot(C0,invtempE1)
#        rho_ = (1.0 - vr)*rhom + vr*rhor
        
#        if vr == vs: 
#            print "vr=vs"
        rho_ = (1.0 - vr - vs)*rhom + vr*rhor + vs*rhos
        lct.append(Ct)
#        lphi.append(vr)
        lphi.append((vr+vs))
        lrho.append(rho_)
#        print Ct[0][0]-2.*Ct[3][3]/3, Ct[3][3]/2.
#        poison = Poisson(Ct[0][0]-2.*Ct[3][3]/3, Ct[3][3]/2.)
#        C0=Ct
#        print 'COCt', poison
#    print 'lenn', vr, vs
    return [lct, lrho, lpoison], lphi
    
def mainIrineuTmatrix():
    import matplotlib.pyplot as plt
    #inclusion
    lines = []
    for j in range(len(vvp)):
#    for j in [9]:
        quartzo = Quartzo()
        calcita = Calcite()
    #    aragonita = Aragonite()
    #    minmatrix = MineralMatriz()
        dolomita = Dolomite()
        clay = Clay()
        dry = Dry()
        water = Water()
        #cavity = 0(isolated)... 1 (exchange)
        cavity = 1
        e=10
        alpha=0.001
        alphad=0.99
        k=np.identity(3)
        omega=10.
    
#        print j, 'jjjjjjjjj', len(vvp)
#        Km = bulk
#        Gm = shear
#        rhom = dens
#        alphameso = macroaspect
#        phimeso = macroporos
#        phimicro = microporos
#        e=10
#        alphamicro=0.001
#        while e>0.001:
        while e>0.01:
            vpo=[]
    #        omega = np.linspace(1,10000000,10000)
    #    for j in [1]:
    #        Km= bulk[j]
    #        Gm= shear[j]   
    #        rhom= dens[j]
    #        meso = macroaspect[j]
    #        micro = microaspect[j]
    #        phimeso = macroporos[j]
    #        phimicro = microporos[j]
    #        print 'jjj', j
            phimax = 0.2303
            Km= bulk[j]
            Gm= shear[j]   
            rhom= dens[j]
#            alphameso = 0.51
#            alphamicro = 0.04
    #    porototal = 0.2303
#            phimeso = 0.141
#            phimicro = 0.0893
            alphameso = macroaspect[j]
            phimeso = macroporos[j]
            phimicro = microporos[j]
            
            minmatrix = MineralMatriz(Km, Gm, rhom)        
    #        fam1 = [water, 0.1, 0.3, 'cav']
            fam1 = [dry, alphameso, phimeso, 'cav']
    #        print meso, phimeso, micro, phimicro
    #        fam1 = [dry, alpha, phimax, 'cav']
            fam2 = [dry, alpha, phimicro, 'cav']
            gr_fam = [fam1, fam2]
            
            if len(gr_fam) == 1:
                composite, phi = IncludeSingle (calcita, gr_fam[0], cavity, j)
    #            print 'n'
            else:
    #            composite, phi = IncludeDual (quartzo, fam1, fam2, cavity)
                composite, phi = IncludeFam (minmatrix, gr_fam, cavity, k, omega)
    #            composite, phi = TMAIrineu (minmatrix, fam1, fam2, cavity, alphad)#, k, omega)
    #            composite, phi = TMASECOND (minmatrix, fam1, fam2, cavity, alphad)
    #            composite, phi = SingleCom (minmatrix, fam2, cavity, alphad, k, omega)#single cavity comunicating
            Ct = composite[0]
            rho = composite[1]
            poisct = composite[2]
            ni = len(phi)
            Vp = np.empty(ni)
            Vs = np.empty(ni)
            Qp = np.empty(ni)
            ks= np.empty(ni)
            kk= np.empty(ni)
            ksg= np.empty(ni)
            g= np.empty(ni)
    #        print len(phi), len(ksg), len(Ct)
    #        print rho,'\n\n'
            for i in range(len(phi)):
                Vp[i] = (((Ct[i][0][0]/rho[i]).real)**0.5)
                Vs[i] = (((Ct[i][3][3]/(2.*rho[i])).real)**-0.5)**-1
    #        vpo.append(Vp[-1])
            alpha= alpha+0.00001
            print 'vpp', vvp[j]
            e = np.abs(Vp[-1]-vvp[j])
            if np.isnan(e): e=10
#            print e, alpha, Vp[-1]#, len(phi), len(Vp), Vp[-1], Vs[-1], Ct[i][0][0], Ct[i][3][3]
        newarq = open("Vp_calc_TMA2.txt", 'w')
        lines.append (str(j+1)+" "+str(Vp[-1])+" "+str(Vs[-1])+ " " + str(alphameso) + " "+str(alpha))
        newarq.writelines ('\n'.join(lines))
    #    plt.plot(phi, kk, 'y')
    #    plt.plot(phi, g, 'r')
    #    plt.plot(phi, rho, 'r')
        
    #    plt.show()
    #        print 'ploto',Vp[-1], vpo
        plt.plot(phi, Vp, 'r') # plot vp
    
    for i in range(len(macroporos)):
        
        if bulk[i] ==0.:
            macroporos.pop(i)
            bulk.pop(i)
        else: 
            plt.plot(macroporos[i]+microporos[i],vvp[i],'o')
    plt.grid(True)
    plt.show()
    
def TMAIrineu (matrix, rfam, sfam, cavity,alphad):  #ni = number of inclusions
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    Rinclusion = rfam[0]
    Ralpha = rfam[1]
    Rphi_ = rfam[2]
    Sinclusion = sfam[0]
    Salpha = sfam[1]
    Sphi_ = sfam[2]
    Cr = Rinclusion[0]
    rhor = Rinclusion[1]
    Cs = Sinclusion[0]
    rhos = Sinclusion[1]
#    alphad = 0.99
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]
#    Rphi_ = 0.5
#    Sphi_ = 0.5
#    Ralpha = 0.99
#    Salpha = 0.01
    poison = poisi
#    Salpha = Ralpha
    Gr, E1 = TensorG(C0, Ralpha, poison) 
#    Gr2 = TensorGreen (Ralpha, C0)
    print 'alphasss\n', Ralpha, Salpha
    Gs, Es = TensorG(C0, Salpha, poison) 
#    dCr = DeltaC(Cr, C0)
    dCr = Cr-C0
    dCs = Cs-C0
    tr = T_r (dCr , Gr)
    ts = T_r (dCs , Gs)    
    Gd, E2 = TensorG (C0, alphad, poison)
#    Gd = TensorGreen (alphad, C0)
    I6 = np.identity(6)
#    vs=vr
    tempE1 = I6 - E1
    tempE2 = I6 - E2
    tempEs = I6 - Es
    invtempE1 = np.linalg.inv(tempE1)
#    invtempE2 = np.linalg.inv(tempE2)
    invtempEs = np.linalg.inv(tempEs)    
#    ts=tr
    
    vr = Rphi_# = Rphi_*0.5
    vs = Sphi_# = Rphi_
#    Rphi_ = Rphi_+0.1
    ni =500
#    for vr in np.linspace(0,Rphi_,ni):
#        vs=vr
#        for vs in np.linspace(0,Sphi_,ni):
#        Gr, E = TensorG(C0, Ralpha, poison)
#        dCr = DeltaC(Cr, C0)
#        tr = T_r (dCr , Gr)
#        alphad = Ralpha/vr
#        print Ralpha/vr, alphad
#        Gd = TensorGreen (alphad, C0)
        
#        Gd, E = TensorG (C0, alphad, poison)
#        Gd = TensorGreen (alphad, C0)
#    for p in range(ni):    
#        vr = p*Rphi_/ni
#        vs = p*Sphi_/ni
#        temp4 = np.dot(E2, invtempEs)
#        temp3 = I6 + vs*temp4
#        invtemp3 = np.linalg.inv(temp3)
#        temp2 = np.dot(invtempE1, invtemp3)
#        temp1 = I6 - vr*temp2
#        Ct = np.dot(C0,temp1)
    
#        vs=vr
#        temp4 = ts*vs
#        temp3 = I6 + np.dot(Gd,temp4)
#        temp2 = np.linalg.inv(temp3)
#        temp1 = vr*tr
#        Ct = C0 + np.dot(temp1, temp2)*2.  
        
#    for p in range(ni):
#        vr = p*Rphi_/ni
#        vs = p*Sphi_/ni
#        temp4 = ts*vs
#        temp3 = I6 + np.dot(Gd,temp4)
#        temp2 = np.linalg.inv(temp3)
#        temp1 = vr*tr
#        
#        temp5 = I6 + np.dot(Gd,temp1)
#        temp6 = np.linalg.inv(temp5)
#        Ct = C0 + np.dot(temp1, temp2) + np.dot(temp4,temp6)
        #novoniterp
#        tr = T_r (dCr , Gr)
#        ts = T_r (dCs , Gs) 
    temp1 = vr*tr
    temp2 = vs*ts
    temp3 = temp1 + temp2
    temp4 = I6 + np.dot(Gd, temp3)
    temp5 = np.linalg.inv(temp4)
    temp6 = np.dot(temp3, temp5)
#    print 'vrrrr',vr,vs
    Ct = C0 + temp6
#        
#        print 'trr', tr, np.dot(C0,invtempE1)
#        rho_ = (1.0 - vr)*rhom + vr*rhor
    rho_ = (1.0 - vr - vs)*rhom + vr*rhor + vs*rhos
    lct.append(Ct)
#        lphi.append(vr)
    lphi.append((vr+vs))
    lrho.append(rho_)
#        C0=Ct
#        print 'COCt', C0[0][0], Ct[0][0], Cr[0][0]
#    print 'lenn', vr, vs
    return [lct, lrho, lpoison], lphi

 
def mainTMADem(matrix, rfam, sfam, cavity,alphad):
    def dCdv( v, c): 
        Gdem = Gdemdcdv
        Cr = Crdcdv
        return (Cr-c)*((np.identity(6)-Gdem*(Cr-c))**-1)/(1.-v)
    def rungeKutta(x0, y0, x, h): 
        
        # Count number of iterations using step size or 
        # step height h 
        n = np.absolute((int)((x - x0)/h))
        # Iterate for number of iterations 
        y = y0 
#        print 'rugge',y0,n
        for i in range(1, n + 1): 
#            print "Apply Runge Kutta Formulas to find next value of y"
            k1 = h * dCdv(x0, y) 
            k2 = h * dCdv(x0 + 0.5 * h, y + 0.5 * k1) 
            k3 = h * dCdv(x0 + 0.5 * h, y + 0.5 * k2) 
            k4 = h * dCdv(x0 + h, y + k3) 
      
            # Update next value of y 
            y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4) 
      
            # Update next value of x 
            x0 = x0 + h 
        return y

    
    C0 = matrix[0] 
    rhom = matrix[1]
    poisi = matrix[2]
    Rinclusion = rfam[0]
    Ralpha = rfam[1]
    Rphi_ = rfam[2]
    Sinclusion = sfam[0]
    Salpha = sfam[1]
    Sphi_ = sfam[2]
    Crdcdv = Rinclusion[0]
    rhor = Rinclusion[1]
    Cs = Sinclusion[0]
    rhos = Sinclusion[1]
#    alphad = 0.99
    lphi = [0]
    lrho = [rhom]
    lpoison = [poisi]
    lct = [C0]
#    Rphi_ = 0.5
#    Sphi_ = 0.5
#    Ralpha = 0.99
#    Salpha = 0.01
    poison = poisi
    Gdemdcdv, E2 = TensorG (C0, alphad, poison)
    
    x0 = 7.
    y = 0.9
    x = 5.
    h = 0.01
    print 'The value of y at x is:', C0[0][0], rungeKutta(x0, y, x, h)
  
# Finds value of y for a given x using step size h 
# and initial value y0 at x0. 


def plotdata ():  #ni = number of inclusions
    import matplotlib.pyplot as plt
    import numpy as np
    crackdens = np.empty(len(macroporos))
    for i in range(len(macroporos)): 
            if 1==0:
                crackdens[i] = 00.
            else:
                crackdens[i] = (3.*(macroporos[i]+microporos[i]))/(4.*np.pi*(macroaspect[i]+microaspect[i]))
#    print a10_plug_perm, type (a10_plug_perm), '\n\n',np.asarray(a10_plug_perm)
    poros = np.asarray(a10_plug_poros)
    perm = np.asarray(a10_plug_perm)
    plt.semilogy(poros,perm,'ob')
#    for i in range(len(a10_plug_perm)):        
#    for i in range(len(macroporos)):        
#        if bulk[i] ==0.:
#            macroporos.pop(i)
#            bulk.pop(i)
            
#        else: 
#            plt.plot(macroporos[i]+microporos[i],vvp[i],'o')
#        if a10_plug_perm[i]>100:
#            plt.plot(a10_plug_poros[i],a10_plug_perm[i],'ob')
#        plt.semilogy(a10_plug_poros[i],a10_plug_perm[i],'ob')
#            plt.plot(crackdens[i],vvp[i],'o')
    plt.grid()
#    plt.show()

if __name__ == '__main__':
    import numpy as np
#    main4()
#    main5()
#    main6()
#    main7()
#    plotdata()
#    main4()
    main1()
#    mainsat()
#    mainbiot()    
    Km= 76.85
    Gm= 32.03    
    rhom= 2.71
    kf = 2.2
    rhof = 1.1
    rock=[Km,Gm,rhom]
    fluid= [kf,rhof]
    phit= 0.5
    rhob = 2.1
    perm = np.identity(3)
    method = 'second'
#    omega = 4000 #well log
    inclusion1 = (0.5, 0.1)
    inclusion2 = (0.5, 0.1)
    alphad=0.1
#    tmatrix_log(rock, fluid, phit, rhob, perm, inclusion1, inclusion2, alphad, method)
#    mainIrineuDem()
#    mainIrineuTmatrix()

    dem_log(rock, fluid, phit, rhob, inclusion1, inclusion2)
#    rock=(km,gm,rhom)
#    fluid= (kf,rhof)
#    phit= log
#    rhob = log
#    inclusion = (fraction, alpha)
        
        
#    import matplotlib.pyplot as plt
#    plt.figure()
#    plt.subplot(211)
#    a= [1,3,5]
#    b= [5,5,5]
#    plt.plot(a,b)
#    plt.subplot(212)
#    a= [1,3,5]
#    b= [5,5,5]
#    plt.plot(a,a)
#    plt.show()
#    i4 = Id4()
#    print KelvintoVoigt(VoigttoKelvin(TensortoVoigt(i4)))
    
    #rho ok, G ok, K nao ok