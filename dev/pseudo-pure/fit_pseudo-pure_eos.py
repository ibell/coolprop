import numpy as np
from CoolProp.CoolProp import Props

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import scipy.optimize
import scipy.stats

import h5py
from templates import *

from Helmholtz import helmholtz

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    return r_value**2
                
def get_fluid_constants(Ref):
    if Ref == 'R407F':
        RefString = 'REFPROP-MIX:R32[0.47319469]&R125[0.2051091]&R134a[0.32169621]'
        
        #Values from R407C
##         N0 = np.array([0,1.0588, -1.12018, 0.629064, -0.351953, 0.00455978, -1.75725, -1.12009, 0.0277353, 0.898881, -1.17591, 0.0818591, -0.0794097, -0.000014047, 0.233779, -0.29179, 0.0154776, -0.0314579, -0.00442552, -0.0101254, 0.00915953, -0.00361575])
        
##         T = [0.6,1,5,1,7,0.3,0.7]
##         D = [1.0,1,1,2,2,4,6]
##         L = [0.0,0,0,0,0,0,0]
        
#         D = []
#         L = []
#         T = []
#         
#         for l in range(4):
#             for d in range(9):
#                 for t in range(4*l, 30, 5):
#                     L.append(float(l))
#                     D.append(float(d))
#                     T.append(float(t)/8.0)
#                
#         D0 = np.array(D)
#         L0 = np.array(L)
#         T0 = np.array(T)
#         N0 = 0.1*np.ones_like(D0)
        
        #Values from Span short(2003) (polar)
        D0 = np.array([0, 1.0, 1.0, 1.0, 3.0, 7.0, 1.0, 2.0, 5.0, 1.0, 1.0, 4.0, 2.0])
        L0 = np.array([0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0])
        T0 = np.array([0, 0.25,  1.25,  1.5, 0.25,  0.875, 2.375, 2.0,   2.125, 3.5,   6.5,   4.75,  12.5])
        N0 = 0.5*np.ones_like(D0)
        
##         D0 = np.array([0,1.0, 1, 1, 2, 5, 1, 2, 2, 3, 3, 5, 5, 5, 1, 1, 4, 4, 2, 4, 5, 6])
##         L0 = np.array([0,0.0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3])
        
##         # values from R410A
##         N0 = np.array([0.0, 0.987252, -1.03017, 1.17666, -0.138991, 0.00302373, -2.53639, -1.96680, -0.830480, 0.172477, -0.261116, -0.0745473, 0.679757, -0.652431, 0.0553849, -0.0710970, -0.000875332, 0.0200760, -0.0139761, -0.0185110, 0.0171939, -0.00482049])
##         T0 = np.array([0.0,0.44,1.2,2.97,2.95,0.2,1.93,1.78,3.0,0.2,0.74,3.0,2.1,4.3,0.25,7.0,4.7,13.0,16.0,25.0,17.0,7.4])
##         D0 = np.array([0,1.0,1,1,2,5,1,2,3,5,5,5,1,1,4,4,9,2,2,4,5,6])
##         L0 = np.array([0,0.0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3])
    elif Ref == 'R410A':
        RefString = 'REFPROP-MIX:R32[0.6976147]&R125[0.3023853]'
        # values from R410A
        N0 = np.array([0.0, 0.987252, -1.03017, 1.17666, -0.138991, 0.00302373, -2.53639, -1.96680, -0.830480, 0.172477, -0.261116, -0.0745473, 0.679757, -0.652431, 0.0553849, -0.0710970, -0.000875332, 0.0200760, -0.0139761, -0.0185110, 0.0171939, -0.00482049])
        T0 = np.array([0.0,0.44,1.2,2.97,2.95,0.2,1.93,1.78,3.0,0.2,0.74,3.0,2.1,4.3,0.25,7.0,4.7,13.0,16.0,25.0,17.0,7.4])
        D0 = np.array([0,1.0,1,1,2,5,1,2,3,5,5,5,1,1,4,4,9,2,2,4,5,6])
        L0 = np.array([0,0.0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3])
        
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
        self.tau = self.Tc/self.T
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
        
        cp0_over_R_check = 1-self.tau**2*self.d2phi0_dTau2(self.tau)
        
        plt.plot(self.T, (self.cp0_R_from_fit(a_e)/self.cp0_R-1)*100, '-', self.T, (cp0_over_R_check/self.cp0_R-1)*100, '^')
        plt.xlabel('Temperature [K]')
        plt.ylabel('($c_{p0}/R$ (fit) / $c_{p0}/R$ (REFPROP) -1)*100 [%]')
        plt.savefig('cp0.pdf')
        plt.close()
        
    def d2phi0_dTau2(self, tau):
        d = []
        for _tau in tau:
            #lead term is killed
            d.append(helmholtz.phi0_logtau(-1.0).dTau2(_tau, _tau)
                    + helmholtz.phi0_cp0_poly(self.a[0],self.e[0],self.Tc,298.15).dTau2(_tau, _tau) 
                    + helmholtz.phi0_Planck_Einstein(self.a,self.e/self.Tc,1,len(self.a)-1).dTau2(_tau, _tau)
                    )
        return np.array(d)

class ResidualPartFitter(object):
    
    def __init__(self, Ref, IPF):
        self.Ref = Ref
        self.IPF = IPF
        self.RefString, self.N0, self.T0, self.D0, self.L0 = get_fluid_constants(Ref)
        self.Tc = Props(self.RefString,'Tcrit')
        self.rhoc = Props(self.RefString,'rhocrit')
        molemass = Props(self.RefString,'molemass')
        self.R = 8.314472/ molemass
        
    def termwise_Rsquared(self):
        
        keepers = []
        values = []
        print len(self.N0), 'terms at start'
        for i in range(len(self.N0)):
            
            n = helmholtz.vectord([float(1)])
            d = helmholtz.vectord([self.D0[i]])
            t = helmholtz.vectord([self.T0[i]])
            l = helmholtz.vectord([self.L0[i]])
            
            self.phir = helmholtz.phir_power(n, d, t, l, 0, 0)
            
            PPF = self.evaluate_EOS(np.array(list(n)+list(t)))
            
            R2 = rsquared(PPF.p,self.phir.dDeltaV(self.tauV,self.deltaV))
            
            values.append((R2,i))
            if R2 > 0.9:
                keepers.append(i)
                
        values,indices = zip(*reversed(sorted(values)))

        keepers = list(indices[0:30])
                
        self.N0 = self.N0[keepers]
        self.T0 = self.T0[keepers]
        self.D0 = self.D0[keepers]
        self.L0 = self.L0[keepers]
        
        print len(self.N0), 'terms at end'
        
    def generate_1phase_data(self):
        
        Tc = Props(self.RefString, 'Tcrit')
        rhoc = Props(self.RefString, 'rhocrit')
        TTT, RHO, PPP, CPP, CVV, AAA = [], [], [], [], [], []
        
        for _T in np.linspace(220, 450, 300):
            print _T
            for _rho in np.logspace(np.log10(1e-2), np.log10(rhoc), 100):
                try:
                    if _T > Tc:
                        p = Props('P', 'T', _T, 'D', _rho, self.RefString)
                        cp = Props('C', 'T', _T, 'D', _rho, self.RefString)
                        cv = Props('O', 'T', _T, 'D', _rho, self.RefString)
                        a = Props('A', 'T', _T, 'D', _rho, self.RefString)
                    else:
                        DL = Props('D', 'T', _T, 'Q', 0, self.RefString)
                        DV = Props('D', 'T', _T, 'Q', 1, self.RefString)
                        if _rho < DV or _rho > DL:
                            p = Props('P', 'T', _T, 'D', _rho, self.RefString)
                            cp = Props('C', 'T', _T, 'D', _rho, self.RefString)
                            cv = Props('O', 'T', _T, 'D', _rho, self.RefString)
                            a = Props('A', 'T', _T, 'D', _rho, self.RefString)
                        else:
                            p = None
                    if p is not None:
                        TTT.append(_T)
                        RHO.append(_rho)
                        PPP.append(p)
                        CPP.append(cp)
                        CVV.append(cv)
                        AAA.append(a)
                        
                except ValueError as VE:
                    print VE
                    pass
            
            for _rho in np.linspace(rhoc, 3.36*rhoc, 50):
                try:
                    if _T > Tc:
                        p = Props('P', 'T', _T, 'D', _rho, self.RefString)
                        cp = Props('C', 'T', _T, 'D', _rho, self.RefString)
                        cv = Props('O', 'T', _T, 'D', _rho, self.RefString)
                        a = Props('A', 'T', _T, 'D', _rho, self.RefString)
                    else:
                        DL = Props('D', 'T', _T, 'Q', 0, self.RefString)
                        DV = Props('D', 'T', _T, 'Q', 1, self.RefString)
                        if _rho < DV or _rho > DL:
                            p = Props('P', 'T', _T, 'D', _rho, self.RefString)
                            cp = Props('C', 'T', _T, 'D', _rho, self.RefString)
                            cv = Props('O', 'T', _T, 'D', _rho, self.RefString)
                            a = Props('A', 'T', _T, 'D', _rho, self.RefString)
                        else:
                            p = None
                    if p is not None:
                        TTT.append(_T)
                        RHO.append(_rho)
                        PPP.append(p)
                        CPP.append(cp)
                        CVV.append(cv)
                        AAA.append(a)
                        
                except ValueError as VE:
                    print VE
                    pass
                    
        h = h5py.File('T_rho_p.h5','w')
        grp = h.create_group(self.Ref)
        grp.create_dataset("T",data = np.array(TTT),compression = "gzip")
        grp.create_dataset("rho", data = np.array(RHO),compression = "gzip")
        grp.create_dataset("p", data = np.array(PPP),compression = "gzip")
        grp.create_dataset("cp", data = np.array(CPP),compression = "gzip")
        grp.create_dataset("cv", data = np.array(CVV),compression = "gzip")
        grp.create_dataset("speed_sound", data = np.array(AAA),compression = "gzip")
        h.close()
    
    def load_data(self):
        h = h5py.File('T_rho_p.h5','r')
        self.T = h.get(self.Ref + '/T').value
        self.rho = h.get(self.Ref + '/rho').value
        self.p = h.get(self.Ref + '/p').value
        self.cp = h.get(self.Ref + '/cp').value
        self.cv = h.get(self.Ref + '/cv').value
        self.speed_sound = h.get(self.Ref + '/speed_sound').value
        
        self.tau = self.Tc/self.T
        self.delta = self.rho/self.rhoc
        self.tauV = helmholtz.vectord(self.tau)
        self.deltaV = helmholtz.vectord(self.delta)
        
        # Get the derivative d2phi0_dTau2 from the ideal part fitter
        self.d2phi0_dTau2 = self.IPF.d2phi0_dTau2(self.tau)
    
    def evaluate_EOS(self, N):
        
#         self.phir.n = helmholtz.vectord(N[0:len(N)//2])
#         self.phir.t = helmholtz.vectord(N[len(N)//2::])
        self.phir.n = helmholtz.vectord(N)
        
        dDelta = self.phir.dDeltaV(self.tauV,self.deltaV)
        dTau2 = self.phir.dTau2V(self.tauV,self.deltaV)
        dDelta2 = self.phir.dDelta2V(self.tauV,self.deltaV)
        dDelta_dTau = self.phir.dDelta_dTauV(self.tauV,self.deltaV)
        
        # Evaluate the pressure
        p = (self.rho*self.R*self.T)*(1 + self.delta*dDelta)
        # Evaluate the specific heat at constant volume
        cv_over_R = -self.tau**2*(self.d2phi0_dTau2 + dTau2)
        cv = cv_over_R*self.R
        # Evaluate the specific heat at constant pressure
        cp_over_R = cv_over_R+(1+self.delta*dDelta-self.delta*self.tau*dDelta_dTau)**2/(1+2*self.delta*dDelta+self.delta**2*dDelta2)
        cp = cp_over_R*self.R
        # Evaluate the speed of sound
        w = np.sqrt(1000*self.R*self.T*cp_over_R/cv_over_R*(1+2*self.delta*dDelta+self.delta**2*dDelta2))
        
        class stub: pass
        PPF = stub()
        PPF.p = np.array(p, ndmin = 1).T 
        PPF.cp = np.array(cp, ndmin = 1).T 
        PPF.cv = np.array(cv, ndmin = 1).T 
        PPF.w = np.array(w, ndmin = 1).T 
        
        return PPF
    
    def OBJECTIVE(self, N):
        PPF = self.evaluate_EOS(N)
            
##         plt.plot(PPF.p, self.p); plt.show()
##         plt.plot(PPF.cp, self.cp); plt.show()
##         plt.plot(PPF.cv, self.cv); plt.show()
##         plt.plot(PPF.w, self.speed_sound); plt.show()
        w_p = 1.0
        w_cv = 1.0
        w_w = 1.0
        w_cp = 1.0
        w_total = (w_p+w_cv+w_w+w_cp)/4
        
        w_p_norm = w_p/w_total
        w_cv_norm = w_cv/w_total
        w_cp_norm = w_cp/w_total
        w_w_norm = w_w/w_total
        
        residuals = np.r_[(PPF.p/self.p-1),w_cv_norm*(PPF.cv/self.cv-1),w_cp_norm*(PPF.cp/self.cp-1)]#,w_w_norm*(PPF.w**2/self.speed_sound**2-1)]
        RMS = np.sqrt(np.mean(np.power(residuals, 2)))
        
        print 'RMS:',RMS, '% Max',np.max(residuals),'%'
        return RMS
    
    def fit(self):
        
        # Kill off some not as good terms
        #self.termwise_Rsquared()
        
        # Load up the residual Helmholtz term with parameters
        n = helmholtz.vectord(self.N0)
        d = helmholtz.vectord(self.D0)
        t = helmholtz.vectord(self.T0)
        l = helmholtz.vectord(self.L0)
        self.phir = helmholtz.phir_power(n, d, t, l, 1, 12)
        
        # Solve for the coefficients
        Nbounds = [(-10,10) for _ in range(len(self.N0))]
        tbounds = [(-1,30) for _ in range(len(self.T0))]
        self.N = scipy.optimize.minimize(self.OBJECTIVE, np.array(list(self.N0)), bounds = Nbounds, options = dict(maxiter = 50)).x
        #self.N = scipy.optimize.minimize(self.OBJECTIVE, np.array(list(self.N0)+list(self.T0)), method = 'L-BFGS-B', bounds = Nbounds + tbounds, options = dict(maxiter = 100)).x

        # Write the coefficients to HDF5 file
        h = h5py.File('fit_coeffs.h5','w')
        grp = h.create_group(self.Ref)
        grp.create_dataset("n", data = np.array(self.N), compression = "gzip")
        #grp.create_dataset("t", data = np.array(self.N[len(self.N)//2::]), compression = "gzip")
        h.close()
    
    def evaluate_REFPROP(self, Ref, T, rho):
        
        p,cp,cv,w = [],[],[],[]
        R = 8.314472/Props(Ref,'molemass')
        for _T,_rho in zip(T, rho):
            p.append(Props("P",'T',_T,'D',_rho,Ref))
            cp.append(Props("C",'T',_T,'D',_rho,Ref))
            cv.append(Props("O",'T',_T,'D',_rho,Ref))
            w.append(Props("A",'T',_T,'D',_rho,Ref))
            
        class stub: pass
        PPF = stub()
        PPF.p = np.array(p, ndmin = 1).T 
        PPF.cp = np.array(cp, ndmin = 1).T 
        PPF.cv = np.array(cv, ndmin = 1).T 
        PPF.w = np.array(w, ndmin = 1).T 
        
        return PPF
    
    def check(self):
        # Load the coefficients from file
        h = h5py.File('fit_coeffs.h5','r')
        grp = h.get(self.Ref)
        n = grp.get('n').value
        h.close()
        
        import matplotlib.colors as colors
        cNorm  = colors.LogNorm(vmin=1e-3, vmax=50)
        PPF = self.evaluate_EOS(np.array(list(n)+list(self.T0)))
        
        SC1 = plt.scatter(self.rho, self.T, s = 8, c = np.abs(PPF.p/self.p-1)*100, edgecolors = 'none', cmap = plt.get_cmap('jet'), norm = cNorm)
        plt.gca().set_xscale('log')
        cb = plt.colorbar()
        cb.set_label('np.abs(PPF.p/self.p-1)*100')
        plt.savefig('pressure.png')
        plt.show()
        
        SC1 = plt.scatter(self.rho, self.T, s = 8, c = np.abs(PPF.cp/self.cp-1)*100, edgecolors = 'none', cmap = plt.get_cmap('jet'), norm = cNorm)
        plt.gca().set_xscale('log')
        cb  = plt.colorbar()
        cb.set_label('np.abs(PPF.cp/self.cp-1)*100')
        plt.savefig('cp.png')
        plt.show()

##         plt.plot(self.T,PPF.p/self.p,'.'); plt.show()
##         plt.plot(self.T,PPF.cp/self.cp,'.'); plt.show()
##         plt.plot(self.T,PPF.cv/self.cv,'.'); plt.show()
##         plt.plot(self.T,PPF.w/self.speed_sound,'.'); plt.show()
    
class PPFFitterClass(object):
    
    def __init__(self, Ref, regenerate_data = False, fit = True):
        
        self.Ref = Ref
        
        self.IPF = IdealPartFitter(Ref)
        self.IPF.fit()
    
        self.RPF = ResidualPartFitter(Ref, IPF = self.IPF)
        if regenerate_data:
            self.RPF.generate_1phase_data()
            
        self.RPF.load_data()
        
        if fit:
            self.RPF.fit()
        self.RPF.check()
        self.output_files()

    def contour_plot(values):
        """
        Parameters
        ----------
        values : iterable, same size as T and rho
        """
        
        plt.semilogx(self.RPF.rho,self.RPF.T,'o')
        plt.show()

        # Generate a regular grid to interpolate the data.
        xi = np.linspace(min(self.RPF.T), max(self.RPF.T), 100)
        yi = np.linspace(min(self.RPF.rho), max(self.RPF.rho), 100)
        xi, yi = np.meshgrid(xi, yi)
        # Interpolate using delaunay triangularization 
        zi = mlab.griddata(np.array(self.RPF.T),np.array(self.RPF.rho),np.array(values),xi,yi)
        cont = plt.contourf(yi,xi,zi,30)
        plt.colorbar()
        plt.show()
        
    def output_files(self):
        h = h5py.File('fit_coeffs.h5','r')
        n = h.get(self.Ref+'/n').value
        #t = h.get(self.Ref+'/t').value

        # Output the header file
        header = PPF_h_template.format(Ref = self.Ref, RefUpper = self.Ref.upper())
                                    
        acoeffs = '0, '+', '.join(['{a:0.6f}'.format(a=_) for _ in self.IPF.a])
        # First one doesn't get divided by critical temperature, later ones do        
        bcoeffs = '0, '
        bcoeffs += str(self.IPF.e[0])+', '
        bcoeffs += ', '.join(['{b:0.4f}/{Tcrit:g}'.format(b=_,Tcrit = self.IPF.Tc) for _ in self.IPF.e[1::]])
            
        ncoeffs = ', '.join(['{a:0.6g}'.format(a=_) for _ in n])
        tcoeffs = ', '.join(['{a:0.6g}'.format(a=_) for _ in self.RPF.T0])
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
        
    

    
if __name__=='__main__':
    Ref = 'R407F'
    
    PPFFitterClass(Ref)

    
    
    
    
    
    

    