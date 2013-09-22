import numpy as np
from CoolProp.CoolProp import Props
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.optimize

def get_fluid_constants(Ref):
    if Ref == 'R407F':
        RefString = 'REFPROP-MIX:R32[0.47319469]&R125[0.2051091]&R134a[0.32169621]'
        # values from R410A
        N0 = np.array([0.0, 0.987252, -1.03017, 1.17666, -0.138991, 0.00302373, -2.53639, -1.96680, -0.830480, 0.172477, -0.261116, -0.0745473, 0.679757, -0.652431, 0.0553849, -0.0710970, -0.000875332, 0.0200760, -0.0139761, -0.0185110, 0.0171939, -0.00482049])
        T0 = np.array([0.0,0.44,1.2,2.97,2.95,0.2,1.93,1.78,3.0,0.2,0.74,3.0,2.1,4.3,0.25,7.0,4.7,13.0,16.0,25.0,17.0,7.4])
        D0 = np.array([0,1,1,1,2,5,1,2,3,5,5,5,1,1,4,4,9,2,2,4,5,6])
        L0 = np.array([0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3])
    elif Ref == 'R410A':
        RefString = 'REFPROP-MIX:R32[0.6976147]&R125[0.3023853]'
        # values from R410A
        N0 = np.array([0.0, 0.987252, -1.03017, 1.17666, -0.138991, 0.00302373, -2.53639, -1.96680, -0.830480, 0.172477, -0.261116, -0.0745473, 0.679757, -0.652431, 0.0553849, -0.0710970, -0.000875332, 0.0200760, -0.0139761, -0.0185110, 0.0171939, -0.00482049])
        T0 = np.array([0.0,0.44,1.2,2.97,2.95,0.2,1.93,1.78,3.0,0.2,0.74,3.0,2.1,4.3,0.25,7.0,4.7,13.0,16.0,25.0,17.0,7.4])
        D0 = np.array([0,1,1,1,2,5,1,2,3,5,5,5,1,1,4,4,9,2,2,4,5,6])
        L0 = np.array([0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3])
        
    return RefString, N0, T0, D0, L0

class IdealPartFitter():
    def __init__(self, Ref):
        self.RefString, N0, T0, D0, L0 = get_fluid_constants(Ref)
        self.T = np.linspace(100, 450, 200)
        self.C = Props('C', 'T', self.T, 'D', 1e-15, self.RefString)
        self.cp0_R = self.C/(8.314472/Props(self.RefString,'molemass'))
        
    def cp0_R_from_fit(self, a_e):
        a = a_e[0:len(a_e)//2]
        e = a_e[len(a_e)//2::]
        u1 = e[1]/self.T
        u2 = e[2]/self.T
        u3 = e[3]/self.T
        return a[0]*self.T**e[0]+a[1]*u1**2*np.exp(u1)/(np.exp(u1)-1)**2+a[2]*u2**2*np.exp(u2)/(np.exp(u2)-1)**2+a[3]*u3**2*np.exp(u3)/(np.exp(u3)-1)**2
        
    def OBJECTIVE_cp0_R(self, a_e):
        cp0_R_fit = self.cp0_R_from_fit(a_e)
        RMS = np.sqrt(np.mean(np.power((self.cp0_R-cp0_R_fit)/self.cp0_R, 2)))
        return RMS
        
    def fit(self):
        a_e = [2.8749, 2.0623, 5.9751, 1.5612, 0.1, 697.0, 1723.0, 3875.0]

        a_e = scipy.optimize.minimize(self.OBJECTIVE_cp0_R, a_e).x
        
        plt.plot(self.T, (self.cp0_R_from_fit(a_e)/self.cp0_R-1)*100, 'o-')
        plt.xlabel('Temperature [K]')
        plt.ylabel('($c_{p0}/R$ (fit) / $c_{p0}/R$ (REFPROP) -1)*100 [%]')
        plt.show()

def generate_1phase_data(Ref):
    RefString, N0, T0, D0, L0 = get_fluid_constants(Ref)
    
    Tc = Props(RefString, 'Tcrit')
    rhoc = Props(RefString, 'rhocrit')
    TTT, RHO, PPP = [], [], []
    
    for _T in np.linspace(220, 450, 300):
        print _T
        for _rho in np.logspace(np.log10(1e-10), np.log10(2.5*rhoc), 300):
            try:
                if _T > Tc:
                    p = Props('P', 'T', _T, 'D', _rho, RefString)
                else:
                    DL = Props('D', 'T', _T, 'Q', 0, RefString)
                    DV = Props('D', 'T', _T, 'Q', 1, RefString)
                    if _rho < DV or _rho > DL:
                        p = Props('P', 'T', _T, 'D', _rho, RefString)
                    else:
                        p = None
                if p is not None:
                    TTT.append(_T)
                    RHO.append(_rho)
                    PPP.append(p)
            except ValueError as VE:
                print VE
                pass
                
    
    import h5py
    h = h5py.File('T_rho_p.h5','w')
    grp = h.create_group(Ref)
    grp.create_dataset("T",data = np.array(TTT),compression = "gzip")
    grp.create_dataset("rho", data = np.array(RHO),compression = "gzip")
    grp.create_dataset("p", data = np.array(PPP),compression = "gzip")
    h.close()
    
def load_data(Ref):
    
    import h5py
    h = h5py.File('T_rho_p.h5','r')
    T = h.get(Ref + '/T').value
    rho = h.get(Ref + '/rho').value
    p = h.get(Ref + '/p').value
    return T,rho,p

def pressure_from_PPF(T,rho):
    RefString, N0, T0, D0, L0 = get_fluid_constants(Ref)
    Tc = Props(RefString, 'Tcrit')
    rhoc = Props(RefString, 'rhocrit')
    
    p = []
    R = 8.314472/Props(Ref,'molemass')
    for _T,_rho in zip(T,rho):
        delta = _rho/rhoc
        tau = Tc/_T
        p.append(Props("P",'T',_T,'D',_rho,Ref))
    
    return p
    
def pressure_from_EOS_wrap(N, Ref):
    return pressure_from_EOS(Ref, N=N)
    
def pressure_from_EOS(Ref, T, rho, N = None):
    
    RefString, N0, T0, D0, L0 = get_fluid_constants(Ref)
    
    if N is not None:
        PPF.n = N[0:len(N)//2]
        PPF.t = N[len(N)//2::]
    
    ddD1 = PPF.dphir_dDelta()

    Tc = Props(RefString, 'Tcrit')
    rhoc = Props(RefString, 'rhocrit')
    R = 8.314472/Props(RefString,'molemass')
    
    p = (rho*R*T)*(1+rho/rhoc*ddD1)
    return np.array(p,ndmin=1).T
    
def OBJECTIVE(N, Ref, pMixture):
    pPPF = pressure_from_EOS(Ref, T, rho, N = N)
    RMS = np.sqrt(np.mean(np.power((pPPF-pMixture)/pMixture, 2)))
    print RMS
    return RMS
    
def fit_residual_part(Ref):
    
    RefString, N0, T0, D0, L0 = get_fluid_constants(Ref)
    
    return scipy.optimize.minimize(OBJECTIVE, np.array(list(N0)+list(T0)), args = (Ref,p), options = dict(maxiter = 20)).x
    
if __name__=='__main__':
    Ref = 'R407F'
    
    IPF = IdealPartFitter(Ref)
    IPF.fit()
    
    quit()
    
    generate_1phase_data(Ref)
    
    T, rho, p = load_data(Ref)
    
##     print max(T), min(T)
    
##     # Generate a regular grid to interpolate the data.
##     xi = np.linspace(min(T), max(T), 100)
##     yi = np.linspace(min(rho), max(rho), 100)
##     xi, yi = np.meshgrid(xi, yi)
##     # Interpolate using delaunay triangularization 
##     zi = mlab.griddata(np.array(T),np.array(rho),np.array(p),xi,yi)
##     cont = plt.contourf(yi,xi,zi,30)
##     plt.colorbar()
##     plt.show()

    plt.plot(rho,T,'o')
    plt.show()
    
    RefString, N0, T0, D0, L0 = get_fluid_constants(Ref)
    
    from summer import PPF_summer
    Tc = Props(RefString, 'Tcrit')
    rhoc = Props(RefString, 'rhocrit')
    PPF = PPF_summer(Tc/T,rho/rhoc)
    PPF.set_constants(T0,D0,L0,6)
    
    N = fit_residual_part(Ref)
    
    print N
    
    pEOS = pressure_from_EOS(Ref, T, rho, N = N)
    
    plt.plot(p,pEOS,'o')
    plt.show()
    