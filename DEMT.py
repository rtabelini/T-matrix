"""
From Berryman 1980
"""
from scipy.integrate import odeint
import numpy as np

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
    dphii = phii/ni
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
#        print K_+4./3*G_
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
    
    
def main1():
    import matplotlib.pyplot as plt
    import Gassmann
#    Vpm = 5.85 # km/s
#    Vsm = 3.9 # km/s
#    rhom = 2.65 # g/cm3
#    Vpf = 1.6 # km/s
    rhof = 1.1 # g/cm3
    phimax = 0.3
    
#    Gm = Vsm*Vsm*rhom # GPa
#    Km = Vpm*Vpm*rhom - 4.0*Gm/3.0 # GPa
#    Kf = Vpf*Vpf*rhof # GPa
    Km= 76.8 #calcite
    Gm= 32.
    rhom= 2.71
    Kf = 2.2
#    Gm2= 94.9 #dolomite
#    Km2= 38.8
#    rhom2= 2.92
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
    
    for alpha in [0.2]:
#    for alpha in [0.2, 0.4, 0.6, 0.8]:
#    for p in [0.2, 0.4, 0.6, 0.8,1]:
#        for alpha in [0.15, 0.02, 0.8]:
#        alphai = alpha
#        if alpha != 0.15:
#            alphai = p*alpha + (1-p)*0.15
        K, G, phi = DEM(Km, Gm, np.array([0.0]), np.array([0.0]), np.array([alpha]), np.array([phimax]))
#        print'kmmm \n', K,'\n\n', G, ' \n mm\n ', Km, Gm
        
#        for i in range(len(phi)):
#            if alpha ==0.15: 
#                phi_ref = phi
#            else:
#                phi = 0.8*phi + 0.2*phi_ref

        rho = (1.0 - phi)*rhom + phi*rhof
        Ks = Gassmann.Ks(K, Km, Kf, phi)
#        print Ks[0]+4./3.*G[0]
        Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
        Vs = np.sqrt(G/rho)
        print rho

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
#            plt.plot(phi, Vp, 'b-.')
#        plt.plot(phi, K, 'b')
#        plt.plot(phi, G, 'g')
#        plt.plot(phi, rho, 'g')
        plt.plot(phi, Vp, 'b')
#        plt.plot(phi, Vs, 'g')
    plt.ylim([1.0, 7.0])
    plt.xlim([0.0, 0.45])
    plt.grid()
    plt.show()

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
    Km = 76.8
    Gm = 32.0
    Kf = 2.2
    # alpha_mean = [0.2] # [0.1, 0.2]
    # alpha_std = np.linspace(0.0025, 0.0525, 21) # np.linspace(0.0025, 0.0275, 11)
    alpha_mean = [0.5]
    alpha_std = np.linspace(0.0025, 0.1325, 14)
    
    x, p = discrete_normal(-3.5, 3.5, 70001)
    phis = p*phimax
    for am in alpha_mean:
        
        K, G, phi = DEM(Km, Gm, np.array([0.0]), np.array([0.0]), np.array([am]), np.array([phimax]), minalphaphiratio=10000)
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
#            Ks = Gassmann.Ks(K, Km, Kf, phi)
            Vp = np.sqrt((K + 4.0*G/3.0)/rho)
            Vs = np.sqrt(G/rho)
            plt.plot(phi, Vp, 'b')
            plt.plot(phi, Vs, 'g')
    
#    plt.ylim([3.5, 7.5])
#    plt.xlim([0.0, 0.35])
    plt.show()

def Calcite ():
    Kmcal = 76.8
    Gmcal = 32.0
    rhocal = 2.71
    Ccal = np.zeros((6,6))
    Ccal[0][0]=Ccal[1][1]=Ccal[2][2] = Kmcal +4.*Gmcal/3
    Ccal[3][3]=Ccal[4][4]=Ccal[5][5] = Gmcal
    Ccal[0][1]=Ccal[0][2]=Ccal[1][0]=Ccal[1][2]=Ccal[2][0]=Ccal[2][1] = Kmcal - 2.*Gmcal/3
#    phical = 0.2
#    alphacal = 0.2
    poiscal = (3.*Kmcal-2.*Gmcal)/(2.*(3.*Kmcal+Gmcal))    
    Ccal_ = VoigttoKelvin(Ccal)
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
    poisdol = (3.*Kmdol-2.*Gmdol)/(2.*(3.*Kmdol+Gmdol))
    Cdol_ = VoigttoKelvin(Cdol)
    return [Cdol_, rhodol, poisdol]

def Aragonite ():
    Kmara = 44.8
    Gmara = 38.8
    rhoara = 2.92
    Cara = np.zeros((6,6))
    Cara[0][0]=Cara[1][1]=Cara[2][2] = Kmara +4.*Gmara/3
    Cara[3][3]=Cara[4][4]=Cara[5][5] = Gmara
    Cara[0][1]=Cara[0][2]=Cara[1][0]=Cara[1][2]=Cara[2][0]=Cara[2][1] = Kmara - 2.*Gmara/3
    poisara = (3.*Kmara-2.*Gmara)/(2.*(3.*Kmara+Gmara))
    Cara_ = VoigttoKelvin(Cara)
    return [Cara_, rhoara, poisara]
    
def Quartzo ():
    Kmq = 37.9
    Gmq = 44.3
    rhoq = 2.65
    Cq = np.zeros((6,6))
    Cq[0][0]=Cq[1][1]=Cq[2][2] = Kmq +4.*Gmq/3
    Cq[3][3]=Cq[4][4]=Cq[5][5] = Gmq
    Cq[0][1]=Cq[0][2]=Cq[1][0]=Cq[1][2]=Cq[2][0]=Cq[2][1] = Kmq - 2.*Gmq/3
    poisq = (3.*Kmq-2.*Gmq)/(2.*(3.*Kmq+Gmq))
    Cq_ = VoigttoKelvin(Cq)
#    print Cq_[0][0]
    return [Cq_, rhoq, poisq]
    
def Clay ():
#    Kmq = 37
#    Gmq = 44
    rhoc = 1.48
    Cclay = np.zeros((6,6))
    Cclay[0][0] = Cclay[1][1] = 17.15
    Cclay[2][2] = 5.26
    Cclay[0][1] = Cclay[1][0] = 17.15-2*6.63
    Cclay[0][2] = Cclay[1][2] = Cclay[2][0] = Cclay[2][1] = 2.71
    Cclay[5][5] = 6.63
    Cclay[4][4] = Cclay[3][3] = 1.48    
#    v1 = 0.2
    poisc = 0.4
    Cclay_ = VoigttoKelvin(Cclay)
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
    viscow = 0.8*(10**-6)   # 0.8 cp * fatorconversao
#    phiwat = 0.2
#    alphawat = 0.2
    poisw = 0.5
    Cw_ = VoigttoKelvin(Cw)
    return [Cw_, rhow, poisw, viscow]

def Dry():
    rhod = 0.001
    viscod = 15*(10**-6)
    Cd = np.zeros((6,6))
    poisd = 0.5
    Cd_ = VoigttoKelvin(Cd)
    return [Cd_, rhod, poisd, viscod]
    
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
    
def Green (alpha, C):
#    E = [Epqrs, Epqsr, Eqprs, Eqpsr]   --> iguais !
#    G = -0.25*(E[0] + E[1] + E[2] + E[3])
    
    C11 = C[0][0] 
    C44 = C[3][3]
    C33 = C[2][2]
    C12 = C[0][1]
    C21 = C[1][0]
    C13 = C[0][2]
    Epqrs = E(alpha, C11, C44, C33, C12, C21, C13)
    G = -1.*Epqrs
    return G
    
def E (alpha, C11, C44, C33, C12, C21, C13):
    from scipy.integrate import quad
    d = C11
    f = C44
    h = C33
    e = (C11-C12)/2
    g = C13 + C44
    def delta(x):
        return 1./( (e*(1-x**2) + f*alpha**2*x**2)*((d*(1-x**2) + f*alpha**2*x**2)*(f*(1-x**2) + h*alpha**2*x**2) - g**2*alpha**2*x**2*(1-x**2)))
    
    def e1111 (x):
        return delta(x)*(1-x**2)*((f*(1-x**2) + h*alpha**2*x**2)*((3*e+d)*(1-x**2) + 4*f*alpha**2*x**2) - 
        g**2*alpha**2*x**2*(1-x**2))
        
    def e3333 (x):
        return delta(x)*alpha**2*x**2*(d*(1-x**2) + f*alpha**2*x**2)*(e*(1-x**2)+f*alpha**2*x**2)
    
    def e1122 (x):
        return delta(x)*(1-x**2)*((f*(1-x**2) + h*alpha**2*x**2)*((e+3*d)*(1-x**2) + 4*f*alpha**2*x**2) -
        3*g**2*alpha**2*x**2*(1-x**2))
        
    def e1133 (x):
        return delta(x)*alpha**2*x**2*(((d+e)*(1-x**2) + 2+f+alpha**2*x**2)*(f*(1-x**2) + h*alpha**2*x**2) -
        g**2*alpha**2*x**2*(1-x**2))
    
    def e3311 (x):
        return delta(x)*(1-x**2)*(d*(1-x**2) + f*alpha**2*x**2)*(e*(1-x**2) + f*alpha**2*x**2)
    
    def e1212 (x):
        return delta(x)*(1-x**2)**2*(g**2*alpha**2*x**2 - (d-e)*(f*(1-x**2) + h*alpha**2*x**2))
    
    def e1313 (x):
        return delta(x)*g*alpha**2*x**2*(1-x**2)*(e*(1-x**2) + f*alpha**2*x**2)

    E1111 = (0.5*np.pi)*quad(e1111, 0, 1)[0]                   #E2222
    E3333 = (4*np.pi)*quad(e3333, 0, 1)[0]
    E1122 = (0.5*np.pi)*quad(e1122, 0, 1)[0]                   #E2211
    E1133 = (2*np.pi)*quad(e1133, 0, 1)[0]                     #E2233
    E3311 = (2*np.pi)*quad(e3311, 0, 1)[0]                     #E3322
    E1212 = (0.5*np.pi)*quad(e1212, 0, 1)[0]
    E1313 = (-2*np.pi)*quad(e1313, 0, 1)[0]                    #E2323
    
    e = np.zeros((6,6))
    e[0][0] = e[1][1] = E1111
    e[2][2] = E3333
    e[0][1] = e[1][0] = E1122
    e[0][2] = e[1][2] = E1133
    e[2][0] = e[2][1] = E3311
    e[5][5] = E1212
    e[4][4] = e[3][3] = E1313
    
    return e
    
def TensorG (C0, alpha, v): #phd tesis annex b
    #alpha eh razao de aspecto menor 1
#    S0 = C0**-1
#    print C0, alpha, v
    S0 = np.linalg.inv(C0)
#    print '\n',S0, np.dot(C0,S0)
#    alpha = np.sqrt(alphaquad)
    alphaquad = alpha**2
    q = ((alpha)/((1.-alphaquad)**1.5))*(np.arccos(alpha) - alpha*((1.-alphaquad)**0.5))
    
    S1111 = ((3.*alphaquad)/(8.*(1.-v)*(alphaquad-1.))) + (1./(4*(1-v)))*(1-2.*v-(9./(4*(alphaquad-1.))))*q
    
    S3333 = (0.5/(1.-v))*(1.-2.*v+(3.*alphaquad-1.)/(alphaquad-1.)-(1.-2.*v+(3.*alphaquad)/(alphaquad-1.))*q)
    
    S1122 = (0.25/(1.-v))*(((alphaquad)/(2*(alphaquad-1.)))-(1.-2.*v+(0.75/(alphaquad-1.)))*q)
    
    S1133 = (0.5/(1.-v))*(-(alphaquad/(alphaquad-1.))+0.5*(((3.*alphaquad)/(alphaquad-1.))-(1.-2.*v))*q)
    
    S3311 = (0.5/(1.-v))*(2.*v-1.-(1./(alphaquad-1.)) + (1.-2.*v+((1.5)/(alphaquad-1.)))*q)
    
    S1212 = (0.25/(1.-v))*(((0.5*alphaquad)/(alphaquad-1))+(1.-2.*v-((0.75)/(alphaquad-1)))*q)
    
    S1313 = (0.25/(1.-v))*(1.-2.*v-((alphaquad+1.)/(alphaquad-1.))-0.5*(1.-2.*v-((3*(alphaquad+1))/(alphaquad-1.)))*q)
    
    Sr = np.tensordot(np.zeros((3,3)), np.zeros((3,3)),0)
    
    Sr[0][0][0][0] = Sr[1][1][1][1] = S1111
    Sr[2][2][2][2] = S3333
    Sr[0][0][1][1] = Sr[1][1][0][0] = S1122
    Sr[0][0][2][2] = Sr[1][1][2][2] = S1133
    Sr[2][2][0][0] = Sr[2][2][1][1] = S3311
    Sr[0][1][0][1] = S1212
    Sr[0][2][0][2] = Sr[1][2][1][2] = S1313
#    print Sr
    
    Sr_ = TensortoKelvin(Sr)
#    print '\ns111\n', S1111, S3333, S1122, S1133, S3311, S1212 ,S1313
    Gr = -np.dot(Sr_,S0)
    return Gr
    


#def Tma (C0,  Nr, tsat=0, omega = 30.*np.pi , tal=0.00001, kfl = 30., visc = 0.001, perm = np.identity(3)):    
#    alpha = Nr[0]
#    poros = Nr[1]
#    Cr = Nr[2]
#    Gr = Nr[3]
#    I4_ = Id4()
#    I4 = Kelvin(I4_)
#    S0 = np.linalg.inv(C0)
##    for i in range (len(C0)):
##        for j in range(len(C0[0])):
###            print C0[i][j]
##            if C0[i][j] != 0:
##                S0[i][j] = C0[i][j]**-1
##    I2 = 0
##    Gr = Green (Nr[0], Nr[2])   # aspect ratio rth, C rth
#    
#    tdry = np.dot((Cr-C0),np.linalg.inv((I4 - np.dot(Gr,(Cr-C0)))))
##    print 'tdry \n', C0, '\n', tdry
##    print tdry
#    if tsat != 0:
#        I2x = np.tensordot(np.identity(3), np.identity(3),0)
#        I2 = Kelvin(I2x)
#        kfl = Cr
#        print 'tsat', kfl, 'kfl \n'
#        
##        Kelvin(I2x)
#        kdry = Kdry (Gr, C0, S0, I4)
#        gamma = Gamma(kfl, kdry, S0)
#        
#        theta = Theta(kfl, S0, omega, visc, poros, gamma, tal, kdry, perm)
#        zeta = Zeta(tdry, S0, poros, omega, tal, gamma, I2)
#        qui = Qui(tdry, S0, I2)
##        print 'tasat', (np.complex(1, omega*tal*gamma))
#        
#        print 'tsat\n\n', qui
#        tsat = (np.dot(theta,zeta) + np.complex(0,1)*omega*tal*np.dot(kfl,qui))/(1. + np.complex(0, 1)*omega*tal*gamma)
#    
#    tma = tdry + tsat  #green rth - inversa!
##    print 'tdry \n', np.real(tsat), '\n', tdry
#    return tma
##    tsat = tdry + (Theta*Zeta(r) + np.complex(0, omega*tal*kfl*Qui(r)))/(1 + np.complex(0,omega*tal*Gamma(r)))
    
def Kdry (Gr, C0, S0):
    I4 = np.identity(6)
#    S0 = np.linalg.inv(C0)
    temp1 = np.dot(Gr,C0)
    temp2 = I4 + temp1
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
#    print gamma
    return gamma
        
def Theta (kfl, S0, omega, visc, poros, gamma, tal, kdry, perm, ku=None, kv=None):
#    ku=omega/2650
#    kv=omega/2650
#    perm = 0.2*np.ones((3,3))#*(10**(-12))
    ku=np.array([omega/2650.,0,0])    # ???
    kut = ku.reshape((1,3))
    kv=np.array([0, omega/2650.,0])   # ???
#    iotg = omega*tal*gamma*np.complex(0,1)+1
#    print 'iotg \n', iotg
   # print 'prin', np.outer(np.outer(ku,kv), perm), '\n\n print', np.outer(np.identity(3),np.identity(3))
#    kuv = np.tensordot(ku, kv, 0)
#    kuvperm_ = np.tensordot(kuv, perm,0)
    
#    kuvperm = Kelvin (kuvperm_)
#    print '\n\n agora', np.dot(kfl,((poros*kdry)/(1.+ np.complex(0,1)*omega*tal*gamma))) - np.complex(0,1)*np.dot(kfl,kuvperm)/(omega*visc), '\ntheta'
#    theta = np.linalg.inv(np.dot(kfl,((1. - np.dot(kfl,S0))*((poros)/(1. + np.complex(0,1)*omega*tal*gamma)) + 
#    np.dot(kfl,((poros*kdry)/(1.+ np.complex(0,1)*omega*tal*gamma))) - np.complex(0,1)*np.dot(kfl,kuvperm)/(omega*visc))))
    
    Suuvv = Muv(S0)
    Kduuvv = Muv(kdry)
#    print kdry, Kduuvv
#    temp1 = 1. - np.dot(kfl,S0)
#    print kfl
    temp1 = 1. - kfl*Suuvv
#    print poros
    temp2 = poros/(1.+np.complex(0,1)*omega*tal*gamma)
    temp3 = (poros*Kduuvv)/(1.+np.complex(0,1)*omega*tal*gamma)
#    temp4 = np.complex(0,1)*np.dot(ku, kv)
    
    temp4 = np.dot(kut, perm)
#    temp5 = np.dot(temp4, perm)
    temp5 = np.dot(temp4, kv)
#    print type(temp5[0])
#    temp6 = np.dot(temp5, kfl)
    temp6 = np.complex(0,1)*temp5[0]*kfl
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
    return theta
    
def Thetacor (kfl, S0, omega, visc, poros, gamma, tal, kdry, perm, ku=None, kv=None):
#    ku=omega/2650
#    kv=omega/2650
    ku=np.array([omega/2650.,0,0])    # ???
    kv=np.array([0, omega/2650.,0])   # ???
#    iotg = omega*tal*gamma*np.complex(0,1)+1
#    print 'iotg \n', iotg
   # print 'prin', np.outer(np.outer(ku,kv), perm), '\n\n print', np.outer(np.identity(3),np.identity(3))
#    kuv = np.tensordot(ku, kv, 0)
#    kuvperm_ = np.tensordot(kuv, perm,0)
    
#    kuvperm = Kelvin (kuvperm_)
#    print '\n\n agora', np.dot(kfl,((poros*kdry)/(1.+ np.complex(0,1)*omega*tal*gamma))) - np.complex(0,1)*np.dot(kfl,kuvperm)/(omega*visc), '\ntheta'
#    theta = np.linalg.inv(np.dot(kfl,((1. - np.dot(kfl,S0))*((poros)/(1. + np.complex(0,1)*omega*tal*gamma)) + 
#    np.dot(kfl,((poros*kdry)/(1.+ np.complex(0,1)*omega*tal*gamma))) - np.complex(0,1)*np.dot(kfl,kuvperm)/(omega*visc))))
    
    Suuvv = Muv(S0)
    Kduuvv = Muv(kdry)
#    print kdry, Kduuvv
#    temp1 = 1. - np.dot(kfl,S0)
#    print kfl
    temp0 = np.dot((np.dot(ku,perm)),kv)
    delta = ((kfl*tal)/(poros*visc))*temp0
    
    temp1 = 1. - kfl*Suuvv
    temp2 = poros/(1.+np.complex(0,1)*omega*tal*gamma)
    temp3 = (poros*Kduuvv)/(1.+np.complex(0,1)*omega*tal*gamma)
#    temp4 = np.complex(0,1)*np.dot(ku, kv)
    temp4 = np.dot(ku, perm)
#    temp5 = np.dot(temp4, perm)
    temp5 = np.dot(temp4, kv)
#    temp6 = np.dot(temp5, kfl)
    temp6 = np.complex(0,1)*temp5*kfl
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
    return theta
    
def Zeta (tdry, S0, poros, omega, tal, gamma):
    I2 = np.identity(3)
    i2xi2 = np.tensordot(I2,I2,0)
    i2xi2m = TensortoKelvin(i2xi2)
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
    i2xi2m = TensortoKelvin(i2xi2)
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
    
def VoigttoKelvin (v):              #conversao de notacao voigt para kelvin
    i = np.identity(6)
    i[3][3] = i[4][4] = i[5][5] = np.sqrt(2)
    temp1 = np.dot(i,v)
    temp2 = np.dot(temp1,i)
    return temp2

def TensortoKelvin (t):             #notacao de kelvin para tensor de 4a ordem
    temp1 = TensortoVoigt(t)
    temp2 = VoigttoKelvin(temp1)
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
    
def DeltaC (Cr, C0):
    return Cr - C0

def T_r (deltaC, Gr):
    I4 = np.identity(6) #TensortoKelvin(Id4())
    temp = I4 - np.dot(Gr, deltaC)
#    print 't_r\n', temp,'\n\n'
    invtemp = np.linalg.inv(temp)
    tr = np.dot(deltaC,invtemp)
#    print C0, deltaC
    return tr

def T_dry(C0, Gr):
    I4 = np.identity(6)
    temp = I4 + np.dot(Gr, C0)
#    print 'tdry\n', temp,'\n\n'
    invtemp = np.linalg.inv(temp)
    tdry = -np.dot(C0,invtemp)
    return tdry
    
def T1 (vr, tr):
    T1 = vr*tr
    return T1
        
def T2 (T1, grs, vs, ts, I4):
    temp1 = T1
    temp2 = np.dot(temp1,grs)
    temp3 = np.dot(temp2, ts)
    temp4 = temp3*vs
    return temp4
   
def Ctotal (C0, tr, vr, gr, ts=0., vs=0.):
#    print vs
#    I4_ = Id4()    
#    I4 = Kelvin(I4_)
#    tr = 
    I4 = np.identity(6)
#    I4 = VoigttoKelvin(I4_)
#    print I4
#    print vs
    T1_ = T1(vr, tr)
    T2_ = T2(T1_, gr, vs, ts, I4)    
#    print T1_[0], vr, tr
#    print 'vr e tr', vr, tr
    invT1 = np.linalg.inv(T1_)
    temp1 = np.dot(invT1, T2_)
    temp2 = I4 + temp1
    invtemp2 = np.linalg.inv(temp2)
#    print '\ninv\n', T1_,'\n', T2_
    temp3 = np.dot(T1_,invtemp2)
    Ct = C0 + temp3
#    print C0[0][0]
#    print Ct[0][0], tr[0][0], gr[0][0], T1_[0][0], T2_[0][0]
    return Ct

def fun(cdem, vr, c, g):
    """Define the right-hand side of equation dy/dt = a*y""" 
#    f = a * y
    I=np.identity(6)
#    print c, I, g
    f = ((c-cdem)*((I-np.dot(g,(c-cdem)))**-1.0))/(1-vr)
    return f
    
def Tinclusion (C0, C1, G1, phi_, visc):  #ni = number of inclusions
#    C0 = mineral1[0] 
#    rhom = mineral1[1]
#    C1 = fluido[0]
#    rhoi = fluido[1]
#    poisi = fluido[2]
    kfl = C1[0][0]
#    S0 = np.linalg.inv(C0)
    omega = np.pi
    tal = 0.000001
#    visc = C0[3]          # 
    perm = 0.2*np.ones((3,3))*(10**(-12))
#    perm = 0.2*np.identity(3)*(10**(-12))# 0.2 Darcy * conversor m2
#    perm = np.identity(3)
#    print visco
#    visc = 0.1
#    phi_ = curphi
#    
#    
#    phi_ += step
#    dC = DeltaC(C1, C0)
#    G1 = TensorG(C0, alpha, poisi)
#    tdry = T_r(dC, G1)
#    kdry = Kdry (G1, C0, S0)
#    gamma = Gamma(kfl, kdry, S0)
#    theta = Theta(kfl, S0, omega, visc, phi_, gamma, tal, kdry, perm)
#    zeta = Zeta(tdry, S0, phi_, omega, tal, gamma)
#    qui = Qui(tdry, S0)
#    temp1 = np.dot(theta, zeta)
#    temp2 = np.complex(0,1)*omega*tal*np.dot(kfl, qui)
#    tsat = tdry + (temp1 + temp2)/(1.+np.complex(0,1)*omega*tal*gamma)
#    print tsat
    
#    dry = Dry()
#    dC_ = DeltaC(dry[0], C0)
    
    S0 = np.linalg.inv(C0)
#    G1_ = TensorG(C0, alpha, poisi)
    tdry = T_dry(C0, G1)
#    tdry = T_r(dC_, G1)
#    print tdry1, '\n\n', tdry,'\n\n\n\n\n\n'
#    print dC_ + C0
    
    kdry = Kdry (G1, C0, S0)
    gamma = Gamma(kfl, kdry, S0)
    theta = Theta(kfl, S0, omega, visc, phi_, gamma, tal, kdry, perm)
    zeta = Zeta(tdry, S0, phi_, omega, tal, gamma)
    qui = Qui(tdry, S0)
#    print qui
    temp1 = theta * zeta
    temp2 = np.complex(0,1)*omega*tal*kfl*qui
    
    tsat = tdry + (temp1 + temp2)/(1.+np.complex(0,1)*omega*tal*gamma)
    #, zeta[0][0]
#    print tdry[0][0]
    return tsat
   
#def IncludeCom (mineral1, fluido2, alpha, phimax, curphi, ni):  #ni = number of inclusions
##    fluido2 = Water()
##    fluido2 = Dry()
#    C0 = mineral1[0] 
#    rhom = mineral1[1]
#    poisi = mineral1[2]
#    C1 = fluido2[0]
#    rhoi = fluido2[1]
##    poisi = fluido2[2]
##    G1 = TensorG(C0, alpha, poisi)
##    print G1
#    step = phimax/ni
##    print step
#    phi_ = curphi
#    ts = np.zeros((6,6))
#    vs = np.zeros((6,6))
#    phi = np.zeros(ni)
#    rho = np.zeros(ni)
##    Vp = np.empty(ni)
##    Vs = np.empty(ni)
#    Ct = np.empty((ni,6,6))
#    for i in range (ni):
##        curphi = phi[i]
##        print i
#        phi_ = phi_ + step#*(1.-phi_)
##        print phi_
##        C0[2][0]/(C0[0][0] + C0[0][1])
#        G1 = TensorG(C0, alpha, poisi)
##        dC = DeltaC(C1, C0)
#        tsat = Tinclusion (C0, C1, G1, phi_)        
##        print C0[2][0]/(C0[0][0] + C0[0][1])
##        tdry = 
##        tma = tdry + tsat
#        Ct_ = Ctotal (C0, tsat, step, G1, ts, vs)
#        Ct[i] = Ct_
##        print Ct[i]
#        C0 = Ct_
#        ts = tsat
#        vs = phi_
##        print i, phi[0]
#        phi[i] = phi_
#        
#        rho[i] = (1.0 - phi[i])*rhom + phi[i]*rhoi
##        print rho[i]
##        Vp[i] = ((Ct[i][0][0]/rho)**-0.5)**-1
##        Vs[i] = ((Ct[i][3][3]/rho)**-0.5)**-1
#    Gct = Ct[-1][3][3]*0.5
#    Kct = Ct[-1][0][0] - 4.*Gct/3.
#    poisct = (3.*Kct-2.*Gct)/(2.*(3.*Kct+Gct))
#    print Ct[-1],'\n\n', rho[-1], '\n\n', mineral1[0] 
##    print Ct[0], Ct[500]
#    return [Ct, rho, phi, poisct]

def opa(cdem, vr, c,g):
    """Define the right-hand side of equation dy/dt = a*y""" 
    cdem = cdem.reshape(6,6)
    c = c.reshape(6,6)
    g = g.reshape(6,6)
    
    I=np.identity(np.shape(cdem)[0])
    dc = c-cdem
#    print cdem
    gdc = np.dot(g,dc)
    i_gdc = I-gdc
    invigdc = np.linalg.inv(i_gdc)
    tr = np.dot(dc,invigdc)
#    print tr
    f = tr/(1-vr)
    return f.reshape(-1)
    
def IncludeCom (mineral1, inclusionR, inclusionS, step):  #ni = number of inclusions
    C0 = mineral1[0] 
    rhom = mineral1[1]
    poisi = mineral1[2]
    alpha = inclusionR[0]
#    print step
    C1 = inclusionR[2][0]
    rhoi = inclusionR[2][1]
    visc = inclusionR[2][3]
    phi = step
#    step = phimax/ni
#    phi_ = phi
    Gi = TensorG(C0, alpha, poisi)
    if inclusionS[0] == 0:
#        phi = step
        ts = np.zeros((6,6))
        vs = np.zeros((6,6))
    else:
#        phi = step
#        Gs = TensorG(C0, inclusionS[0], inclusionS[2][2])
        vs = inclusionS[1]
        viscS = inclusionS[2][3]
        ts = Tinclusion (C0, inclusionS[2][0], Gi, vs, viscS)
#        print ts
#    print vs
    Gi = TensorG(C0, alpha, poisi)
#    dC = DeltaC(C1, C0)
    Tr = Tinclusion (C0, C1, Gi, phi, visc)      
#    print Tr
    Ct_ = Ctotal (C0, Tr, phi, Gi, ts, vs)
#    Cdem = C0
#    derCt_ = odeint(fun, Cdem, phi, args=(C0,Gi))#Ctotal (C0, Tr, phi, Gi, ts, vs)      cdem, vr, c, g):
#    Ct_ = derCt_[:,0]
#    Ct[i] = Ct_
    Gct = Ct_[3][3]*0.5
    Kct = Ct_[0][0] - 4.*Gct/3.
    
    rho = (1.0 - phi)*rhom + phi*rhoi
    Vp = ((Ct_[0][0]/rho)**0.5)
    Vs = ((Ct_[3][3]/(2.*rho))**0.5)
    cmin = Ct_
    poisct = (3.*Kct-2.*Gct)/(2.*(3.*Kct+Gct))
    mineral = [cmin, rho, poisct.real]
#    print Ct_,'\n\n', rho, '\n\n', Vp, Vs
#    return [Ct, rho, phimax, poisct, Vp, Vs]
    return Vp, Vs, rho, Kct, Gct, mineral

def MinQuad (x_data, y_data, grade):
#    print min(x_data), max(x_data)
    ni = (max(x_data)-min(x_data))/1000
    t = np.arange(min(x_data), max(x_data), ni)
    
    m =  len(x_data)
    x2 = np.array([x*x for x in x_data]).reshape((m, 1))
    x3 = np.array([x*x*x for x in x_data]).reshape((m, 1))
    x_data = np.reshape(x_data, (m, 1))
    y_data = np.reshape(y_data, (m, 1)) 
    
    ones = np.ones(m).reshape(m, 1)
    
    J1 = np.hstack([x_data, ones])
    J2 = np.hstack([x2, x_data, ones])
    J3 = np.hstack([x3, x2, x_data, ones])
    
    theta_g_1 = np.linalg.lstsq(J1, y_data )
    theta_g_2 = np.linalg.lstsq(J2, y_data )
    theta_g_3 = np.linalg.lstsq(J3, y_data )
#    print min(x_data)
    

    if grade == 1:
        lin = theta_g_1[0][1] + theta_g_1[0][0]*t
        return t, lin
    if grade == 2:
        quad = theta_g_2[0][2] + theta_g_2[0][1]*t + theta_g_2[0][0]*t*t
        return t, quad
    if grade == 3:
        thir = theta_g_3[0][3] + theta_g_3[0][2]*t + theta_g_3[0][1]*t*t + theta_g_3[0][0]*t*t*t
        return t, thir

def main3():
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
    aragonita = Aragonite()
    dolomita = Dolomite()
    clay = Clay()
    dry = Dry()
    water = Water()
    ni = 50
#    C0 = quartzo[0]
#    Ci = calcita[0]
#    alphai = 0.15
    #    poisi = calcita[2]
    
#    print 'ccal',Gi
    #fam = [alpha, poros, element]
    fam0 = [0.00, 0.00, calcita]
    fam1 = [0.20, 0.30, water]
    fam2 = [0.20, 0.10, water]
    fam3 = [0.20, 0.10, water]
#    fam4 = [0.05, 0.01, water]
#    fam5 = [0.01, 0.01, water]
#    fam = [fam0, fam1, fam2, fam3]
    fam = [fam0, fam1]
#    print fam, '\n\n'
#    fam = []
#    for i in np.linspace(0,0.4,20):
#        famx = [0.05, 0.02 , water]
#        fam.append(famx)        
#    print fam        
    porosi = 0.
#    porosf = 0.2
#    alpham = 0.25
#    porosm = 0.1
#    Vp = np.empty(ni)
#    Vs = np.empty(ni)
    Vp = []
    Vs = []
    Rho = []
    
#    k= np.empty(ni)
    K = []
    G = []
    
#    rhom = calcita[1]
##    rhoi = water[1]
#    elemento = element
#    elemento[0] = element[0][-1]
#    elemento[1] = element[1][-1]
#    elemento[2] = element[2][-1]
#    print water,elemento
#    [Ct, rho, phi, poisct] = IncludeIso (calcita, water, alphai, porosf, porosi, ni)
    
#    step = porosf/ni
    Phi=[]
    Mineral=[]
    
    
#    for v in family:
        
#    porosinc = 0.0
#    porosinc = np.linspace(porosi, porosf, ni)
#    porosminc = np.linspace(porosi, porosm, ni)
#    print porosinc
#    print water
    porosi = 0.
    matrix = calcita
#    print ((Ct_[0][0]/rho)**0.5)
    for inclusionR in fam[1:]:
#        print matrix[0][0]
#        if inclusionR != fam[-1]:
#            print 'diferente r e ultima'
        inclusionS = fam[fam.index(inclusionR)-1]
#        else:
#            inclusionS = 0.
#        print 'print dsps\n' , inclusionR, '\n',inclusionS
        porosf = porosi + inclusionR[1]#+inclusionS[1]
        porosinc = np.linspace(porosi, porosf, ni)
#        print porosi, porosf, inclusionR[1], inclusionS[1]
#        porosi = inclusionR[1]
#        print '\n', ma
        for step in porosinc[1:]:
#        for step in range(:
            print porosi, porosf, step, step-porosi
            [vp, vs, rho, k , g, mineral] = IncludeCom (matrix, inclusionR, inclusionS, step-porosi)# (step-porosi)/(1-porosf))# step-porosi)
#            print step, vp, rho, mineral [0][0][0]
#            print step, vp
#    for step in porosinc[1:]:
#        porosf = step+porosinc
#        [vp, vs, rho, k , g, mineral] = IncludeIso (calcita, water, alphai, step)
#        print porosi, step
#            print vp
            Vp.append(vp)
            Vs.append(vs)
            Rho.append(rho) 
            K.append(k)
            G.append(g)
            Mineral.append(mineral)
        #        porosinc = step+porosinc
        #        porosf = step+porosi
#            print porosi, porosf, step, Vp
            Phi.append(step)
#            print matrix, mineral, '\n\n'
        matrix = mineral
        porosi = porosf
        print 'rhoooo', Rho
        
    x, y = MinQuad (Phi, Vp, 2)
    

    import matplotlib.pyplot as plt
#    print Phi, K
    plt.plot(Phi, Vp,'red')
#    plt.show()
#    plt.plot(Phi, Vs)
    plt.grid()
#    plt.show()
#    plt.plot(Phi, Rho)
#    plt.xlim([0.0, 0.35])
#    plt.ylim([0.0, 7.0])
#    plt.show()
    
#    print theta_g_1
#    plt.plot(x, y)
#    plt.show()
    
def main4():
    Cmin = Calcite()
    Cinc = Water()
    C0 = Cmin[0]
    C1 = Cinc[0]
    alpha = 0.2
    v = Cmin[2]
    g = TensorG (C0, alpha, v)
    vr = np.linspace(0, 0.4, 50)
    rho = np.linspace(2.7,2.2,50)
    yy = odeint(opa, C0.reshape(-1), vr, args=(C1.reshape(-1),g.reshape(-1),))
    Ct = []
    vp=[]
    for i in range(len(yy)):
        Ct.append(yy[i].reshape(6,6))
        vp.append((yy[i][0]/rho[i])**0.5)
#        print '\n'*2,yy[i][0],'\n', Ct[i][0], vp[i]
#        print Ct[i][0][0]
#        print np.size(Ct)/36
#    print Ct[0][0]/rho
#    vp = (Ct[0][0]/rho)**0.5)
#    plt.plot(vr,yy[:,0])
    plt.plot(vr,vp)
#    print yy[-1].reshape(6,6)
#    plt.sho2w()
if __name__ == '__main__':
#    import matplotlib.pyplot as plt
#    import numpy as np
#    import Gassmann
#    rhof = 1.1 # g/cm3
#    phimax = 0.4
#    Gm= 76.8 #calcite
#    Km= 32.
#    rhom= 2.71
#    Kf = 2.2
#    alpha=0.1
#    K, G, phi = DEM(Km, Gm, np.array([0.0]), np.array([0.0]), np.array([alpha]), np.array([phimax]))
#    rho = (1.0 - phi)*rhom + phi*rhof
#    Ks = Gassmann.Ks(K, Km, Kf, phi)
#        
#    Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
#    Vs = np.sqrt(G/rho)
#    plt.plot(phi, Vp, 'g')
    
#    plt.plot(phi, Vs, 'g')    
#    plt.grid()
#    plt.show()
#    main3()
    main4()
    main1()