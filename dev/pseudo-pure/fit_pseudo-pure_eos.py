import numpy as np
from CoolProp.CoolProp import Props
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.optimize
import h5py
from templates import *

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

class IdealPartFitter(object):
    def __init__(self, Ref):
        self.Ref = Ref
        self.RefString, N0, T0, D0, L0 = get_fluid_constants(Ref)
        self.molemass = Props(self.RefString,'molemass')
        self.Tc = Props(self.RefString, 'Tcrit')
        self.rhoc = Props(self.RefString, 'rhocrit')
        self.pc = Props(self.RefString, 'pcrit')
        self.T = np.linspace(100, 450, 200)
        self.C = Props('C', 'T', self.T, 'D', 1e-15, self.RefString)
        R = 8.314472/self.molemass
        self.cp0_R = self.C/R
        
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
        
        self.a = a_e[0:len(a_e)//2]
        self.e = a_e[len(a_e)//2::]
        
        plt.plot(self.T, (self.cp0_R_from_fit(a_e)/self.cp0_R-1)*100, 'o-')
        plt.xlabel('Temperature [K]')
        plt.ylabel('($c_{p0}/R$ (fit) / $c_{p0}/R$ (REFPROP) -1)*100 [%]')
        plt.close()

class ResidualPartFitter(object):
    
    def __init__(self, Ref):
        self.Ref = Ref
        self.RefString, self.N0, self.T0, self.D0, self.L0 = get_fluid_constants(Ref)
        
    def generate_1phase_data(self):
        
        Tc = Props(self.RefString, 'Tcrit')
        rhoc = Props(self.RefString, 'rhocrit')
        TTT, RHO, PPP = [], [], []
        
        for _T in np.linspace(220, 450, 300):
            print _T
            for _rho in np.logspace(np.log10(1e-10), np.log10(2.5*rhoc), 300):
                try:
                    if _T > Tc:
                        p = Props('P', 'T', _T, 'D', _rho, self.RefString)
                    else:
                        DL = Props('D', 'T', _T, 'Q', 0, self.RefString)
                        DV = Props('D', 'T', _T, 'Q', 1, self.RefString)
                        if _rho < DV or _rho > DL:
                            p = Props('P', 'T', _T, 'D', _rho, self.RefString)
                        else:
                            p = None
                    if p is not None:
                        TTT.append(_T)
                        RHO.append(_rho)
                        PPP.append(p)
                except ValueError as VE:
                    print VE
                    pass
                    
        h = h5py.File('T_rho_p.h5','w')
        grp = h.create_group(self.Ref)
        grp.create_dataset("T",data = np.array(TTT),compression = "gzip")
        grp.create_dataset("rho", data = np.array(RHO),compression = "gzip")
        grp.create_dataset("p", data = np.array(PPP),compression = "gzip")
        h.close()
    
    def load_data(self):
        h = h5py.File('T_rho_p.h5','r')
        self.T = h.get(self.Ref + '/T').value
        self.rho = h.get(self.Ref + '/rho').value
        self.p = h.get(self.Ref + '/p').value
    
    def pressure_from_EOS(self, N):
        
        self.PPF.n = N[0:len(N)//2]
        self.PPF.t = N[len(N)//2::]
        
        ddD1 = self.PPF.dphir_dDelta()

        Tc = Props(self.RefString, 'Tcrit')
        rhoc = Props(self.RefString, 'rhocrit')
        R = 8.314472/Props(self.RefString,'molemass')
        
        p = (self.rho*R*self.T)*(1+self.rho/rhoc*ddD1)
        return np.array(p,ndmin=1).T
    
    def OBJECTIVE(self, N):
        pPPF = self.pressure_from_EOS(N)
        RMS = np.sqrt(np.mean(np.power((pPPF-self.p)/self.p, 2)))
        print RMS
        return RMS
    
    def fit(self):
        
        from summer import PPF_summer
        Tc = Props(self.RefString, 'Tcrit')
        rhoc = Props(self.RefString, 'rhocrit')
        self.PPF = PPF_summer(Tc/self.T,self.rho/rhoc)
        self.PPF.set_constants(self.T0,self.D0,self.L0,6)
        
        self.N = scipy.optimize.minimize(self.OBJECTIVE, np.array(list(self.N0)+list(self.T0)), options = dict(maxiter = 20)).x
        
        h = h5py.File('fit_coeffs.h5','a')
        grp = h.create_group(self.Ref)
        grp.create_dataset("n", data = np.array(self.N[0:len(self.N)//2]), compression = "gzip")
        grp.create_dataset("t", data = np.array(self.N[len(self.N)//2::]), compression = "gzip")
        h.close()
    
class PPFFitterClass(object):
    
    def __init__(self, Ref, use_saved_data = True):
        
        self.Ref = Ref
        
        self.IPF = IdealPartFitter(Ref)
        self.IPF.fit()
    
        self.RPF = ResidualPartFitter(Ref)
        if not use_saved_data:
            self.RPF.generate_1phase_data()
            
        self.RPF.load_data()

        # Generate a regular grid to interpolate the data.
        xi = np.linspace(min(self.RPF.T), max(self.RPF.T), 100)
        yi = np.linspace(min(self.RPF.rho), max(self.RPF.rho), 100)
        xi, yi = np.meshgrid(xi, yi)
        # Interpolate using delaunay triangularization 
        zi = mlab.griddata(np.array(self.RPF.T),np.array(self.RPF.rho),np.array(self.RPF.p),xi,yi)
        cont = plt.contourf(yi,xi,zi,30)
        plt.colorbar()
        plt.show()
        
        if not use_saved_data:
            self.RPF.fit()
        self.output_files()
        
    def output_files(self):
        h = h5py.File('fit_coeffs.h5','r')
        n = h.get(self.Ref+'/n').value
        t = h.get(self.Ref+'/t').value

        # Output the header file
        header = PPF_h_template.format(Ref = self.Ref, RefUpper = self.Ref.upper())

                                    
        acoeffs = '0, '+', '.join(['{a:0.6f}'.format(a=_) for _ in self.IPF.a])
        # First one doesn't get divided by critical temperature, later ones do        
        bcoeffs = '0, '
        bcoeffs += str(self.IPF.e[0])+', '
        bcoeffs += ', '.join(['{b:0.4f}/{Tcrit:g}'.format(b=_,Tcrit = self.IPF.Tc) for _ in self.IPF.e[1::]])
            
        ncoeffs = ', '.join(['{a:0.6g}'.format(a=_) for _ in n])
        tcoeffs = ', '.join(['{a:0.6g}'.format(a=_) for _ in t])
        dcoeffs = ', '.join(['{a:0.6g}'.format(a=_) for _ in self.RPF.D0])
        lcoeffs = ', '.join(['{a:0.6g}'.format(a=_) for _ in self.RPF.L0])
            
        import sys
        sys.path.append('..')
        from fit_ancillary_ODRPACK import saturation_pressure, saturation_density
        pL = saturation_pressure(self.IPF.RefString, self.IPF.Ref, LV = 'L')
        pV = saturation_pressure(self.IPF.RefString, self.IPF.Ref, LV = 'V')
        rhoL = saturation_density(self.IPF.RefString, self.IPF.Ref, form='A', LV='L', add_critical = False)
        rhoV = saturation_density(self.IPF.RefString, self.IPF.Ref, form='B', LV='V', add_critical = False)
        
        code = PPF_cpp_template.format(Ref = self.Ref,
                                      RefUpper = self.Ref.upper(),
                                      acoeffs = acoeffs,
                                      bcoeffs = bcoeffs,
                                      Ncoeffs = ncoeffs,
                                      tcoeffs = tcoeffs,
                                      dcoeffs = dcoeffs,
                                      Lcoeffs = lcoeffs,
                                      N_phir = len(n),
                                      N_cp0 = len(self.IPF.a),
                                      molemass = self.IPF.molemass,
                                      Ttriple = 200,
                                      accentric = 0.7,
                                      pcrit = self.IPF.pc,
                                      Tcrit = self.IPF.Tc,
                                      rhocrit = self.IPF.rhoc,
                                      pL = pL,
                                      pV = pV,
                                      rhoL = rhoL,
                                      rhoV = rhoV
                                      )
                                      
        f = open(self.IPF.Ref+'.h','w')
        f.write(header)
        f.close()
        
        f = open(self.IPF.Ref+'.cpp','w')
        f.write(code)
        f.close()
        
    
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
    
if __name__=='__main__':
    Ref = 'R407F'
    
    PPFFitterClass(Ref)

    
    
    
    
    
    

    