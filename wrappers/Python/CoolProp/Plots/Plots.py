import numpy as np
import CoolProp.CoolProp as cp
from scipy.interpolate import interp1d
from types import NoneType
from matplotlib import pylab
import math

def InlineLabel(xv,yv,x = None, y= None, axis = None, fig = None):
    """
    This will give the coordinates and rotation required to align a label with
    a line on a plot
    """
    
    def ToPixelCoords(xv,yv,axis,fig):
        [Axmin,Axmax]=axis.get_xlim()
        [Aymin,Aymax]=axis.get_ylim()
        DELTAX_axis=Axmax-Axmin
        DELTAY_axis=Aymax-Aymin
        
        width=fig.get_figwidth()
        height=fig.get_figheight()
        pos=axis.get_position().get_points()
        [[Fxmin,Fymin],[Fxmax,Fymax]]=pos
        DELTAX_fig=width*(Fxmax-Fxmin)
        DELTAY_fig=height*(Fymax-Fymin)
        
        #Convert coords to pixels
        x=(xv-Axmin)/DELTAX_axis*DELTAX_fig+Fxmin
        y=(yv-Aymin)/DELTAY_axis*DELTAY_fig+Fymin
        
        return x,y
    
    def ToDataCoords(xv,yv,axis,fig):
        [Axmin,Axmax]=axis.get_xlim()
        [Aymin,Aymax]=axis.get_ylim()
        DELTAX_axis=Axmax-Axmin
        DELTAY_axis=Aymax-Aymin
        
        width=fig.get_figwidth()
        height=fig.get_figheight()
        pos=axis.get_position().get_points()
        [[Fxmin,Fymin],[Fxmax,Fymax]]=pos
        DELTAX_fig=(Fxmax-Fxmin)*width
        DELTAY_fig=(Fymax-Fymin)*height
        
        #Convert back to measurements
        x=(xv-Fxmin)/DELTAX_fig*DELTAX_axis+Axmin
        y=(yv-Fymin)/DELTAY_fig*DELTAY_axis+Aymin
        
        return x,y
    
    if axis is None:
        axis=pylab.gca()
    
    if fig is None:
        fig=pylab.gcf()
    
    
    
    if y is None and x is not None:
        trash=0
        (xv,yv)=ToPixelCoords(xv,yv,axis,fig)
        #x is provided but y isn't
        (x,trash)=ToPixelCoords(x,trash,axis,fig)
    
        #Get the rotation angle
        f = interp1d(xv, yv)
        y = f(x)
        h = 0.001*x
        dy_dx = (f(x+h)-f(x-h))/(2*h)
        rot = np.arctan(dy_dx)/np.pi*180.
        
    elif x is None and y is not None:
        #y is provided, but x isn't
        
        _xv = xv[::-1]
        _yv = yv[::-1]
        #Find x by interpolation
        x = interp1d(yv, xv)(y)
        trash=0
        (xv,yv)=ToPixelCoords(xv,yv,axis,fig)
        (x,trash)=ToPixelCoords(x,trash,axis,fig)
        
        f = interp1d(xv, yv)
        y = f(x)
        h = 0.001*x
        dy_dx = (f(x+h)-f(x-h))/(2*h)
        rot = np.arctan(dy_dx)/np.pi*180.
        
    (x,y)=ToDataCoords(x,y,axis,fig)
    return (x,y,rot)

def show():
    """
    A convenience function to call pylab.show()
    """
    pylab.show()
    

def drawIsoLines(Ref, iMin, iMax, which='', plot='', axis=None, fig=None):
    
    if axis is None:
        axis=pylab.gca()
        
    if fig is None:
        fig=pylab.gcf()
    
    # define which isolines are allowed for the different plots
    plots = {
      'ts'   : ['p','h','d'],
      'ph'   : ['s','t','d'],
      'hs'   : ['p','t','d'],
      'ps'   : ['t','h','d'],
      'prho' : ['t','s','h'],
      'trho' : ['s','p','h'],
      'pt'   : ['s','p','d']
    }
    
    plot = 'ts'
    which='p'
    
    if type(plot)!=NoneType:
        plot = str(plot).strip().lower()
        if plot in plots:
            if type(which)!=NoneType:
                which = str(which).strip().lower()
                if which in plots[plot]:
                    lines = getLines(Ref, iMin, iMax, which, plot, axis)
                    for line in lines:
                        axis.plot(line['x'],line['y'],color=line['colour'])
                elif which=='all':
                    for l in plots[plot]:
                        drawIsoLines(Ref, which=l, plot=plot, axis=axis, fig=fig)
                else:
                    which=False
            else:
                which=False
        else:
            plot = False
    else: 
        plot = False
    
    if not plot:
        raise ValueError('You have to specify the kind of plot, use Ts, Ph, hs, Ps, Prho, Trho or PT.')
    
    if not which:
        raise ValueError('This kind of line is not supported for your plot. Please choose another one.')
    
def getLines(Ref, iMin, iMax, which, plot, axis):
    
    # Get strings to feed to Props function
    yName = str(plot[0] ).upper()
    xName = str(plot[1:]).upper()
    iName = str(which[0]).upper()
    
    # Get current axis limits, be sure to set those before drawing isolines
    [Axmin,Axmax]=axis.get_xlim()
    [Aymin,Aymax]=axis.get_ylim()
    
    # Determine x range for plotting
    samples = 100
    x0 = np.linspace(Axmin,Axmax,samples)
    
#    patterns = {
#      'P' : [1.,2.5,5],
#      'D' : [1.,2.5,5],
#      'H' : [1.,2.5,5],
#      'T' : [1.],
#      'S' : [1.,2.5,5]
#    }
    
    numberOfLines = 10
    patterns = {
      'P' : (lambda x: np.around(np.logspace(math.log10(x[0]), math.log10(x[1]), num=numberOfLines),decimals=-2)),
      'D' : [1.,2.5,5],
      'H' : [1.,2.5,5],
      'T' : (lambda x: np.around(np.linspace(x[0], x[1], num=numberOfLines))+273.15),
      'S' : [1.,2.5,5]
    }
    
    # Use the y range to determine spacing of isolines    
#    iVal = [cp.Props(iName,yName,Aymin,xName,Axmin,Ref),
#            cp.Props(iName,yName,Aymax,xName,Axmin,Ref),
#            cp.Props(iName,yName,Aymin,xName,Axmax,Ref),
#            cp.Props(iName,yName,Aymax,xName,Axmax,Ref) ]
#    iVal = patterns[iName]([np.min(iVal),np.max(iVal)]) 
    iVal = patterns[iName]([iMin,iMax]) 
    
    # TODO: Determine saturation state if two phase region present
    xVal = [x0 for i in iVal]
    
    (Y,X) = _getI_YX(Ref,iName,xName,yName,iVal,xVal)
    
    lines = []
    for j,a in enumerate(X):
        line = {
          'x' : X[j],
          'y' : Y[j],
          'label' : getIsoLineLabel(iName,iVal[j]),
          'colour': getIsoLineColour(iName)
          }
        lines.append(line)
        
    return lines
        
#    line1 = {
#      'x' : [2,4],
#      'y' : [9,14],
#      'label' : '1st label',
#      'colour': getIsoLineColour('p')
#    }
#    line2 = {
#      'x' : [3,5],
#      'y' : [12,15],
#      'label' : '2nd label',
#      'colour': getIsoLineColour('t')
#    }
#    return [line1,line2]

def _getI_YX(Ref,iName,xName,yName,iVal,xVal):
    
    if (len(iVal)!=len(xVal)):
        raise ValueError('We need the same number of x value arrays as iso quantities.')
    
    y  = []
    x  = []
    for j in range(len(iVal)):
        jiVal = iVal[j] #constant quantity
        jxVal = xVal[j] #x-axis array
        Y = np.array([cp.Props(yName,iName,jiVal,xName,ijxVal,Ref) for ijxVal in jxVal])
        y.append(Y)
        x.append(jxVal)
        
    return x,y


#def getP_TS(Ref, pVals, x0):
#    """ Get isobars for a Ts-plot """
#    y  = []
#    x  = []   
#    for iVal in pVals:
#        xi = x0
#        if (iVal < cp.Props(Ref,'pcrit')):
#            xi = np.append(xi, [cp.Props('S','P',iVal,'Q',0,Ref), cp.Props('S','P',iVal,'Q',1,Ref)])
#            xi = np.sort(xi)
#        (Y,X) = _getI_YX(Ref,'P','S','T',[iVal],[xi])
#        y.append(Y)
#        x.append(X)
#    return x,y
#
#def getH_TS():
#    return [[1],[2]]
#
#def getD_TS():
#    return [[1],[2]]
#
#def getS_PH():
#    return [[1],[2]]
#
#def getT_PH():
#    return [[1],[2]]
#
#def getD_PH():
#    return [[1],[2]]
#
#def getP_HS():
#    return [[1],[2]]
#
#def getT_HS():
#    return [[1],[2]]
#
#def getD_HS():
#    return [[1],[2]]
#
#def getT_PS():
#    return [[1],[2]]
#
#def getH_PS():
#    return [[1],[2]]
#
#def getD_PS():
#    return [[1],[2]]
#
#def getT_PRHO():
#    return [[1],[2]]
#
#def getS_PRHO():
#    return [[1],[2]]
#
#def getH_PRHO():
#    return [[1],[2]]
#
#def getS_TRHO():
#    return [[1],[2]]
#
#def getP_TRHO():
#    return [[1],[2]]
#
#def getH_TRHO():
#    return [[1],[2]]
#
#def getS_PT():
#    return [[1],[2]]
#
#def getP_PT():
#    return [[1],[2]]
#
#def getD_PT():
#    return [[1],[2]]
   
      

def getIsoLineColour(which):
    colourMap = {                 
      't' : 'red',
      'p' : 'grey',
      'h' : 'green',
      'd' : 'blue',
      's' : 'yellow'
    }
    return colourMap[str(which).lower()]

def getIsoLineLabel(which,num):
    labelMap = {                 
      't' : [r'$T = ','$ K'],
      'p' : [r'$p = ','$ kPa'],
      'h' : [r'$h = ','$ kJ/kg'],
      'd' : [r'$\rho = ','$ kg/m$^3$'],
      's' : [r'$s = ','$ kJ/kg-K']
    }
    l = labelMap[str(which).lower()]
    return l[0]+str(num)+l[1]


##    if len(pressures)==0:
##        #Calculate pressures
##        [Axmin,Axmax]=axis.get_xlim()
##        [Aymin,Aymax]=axis.get_ylim()
#
##def drawIsobars(pressures=[], plot='', axis=None, fig=None):
    

    
    
def Ts(Ref,Tmin = None, Tmax = None, show=False, axis=None, **kwargs):
    """
    Make a temperature-entropy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """

    ax = axis if axis is not None else pylab.gca()
    if Tmin is None:
        Tmin = cp.Props(Ref,'Tmin')
    if Tmax is None:
        Tmax = cp.Props(Ref,'Tcrit')-1e-5
        
    if Tmin > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmin cannot be greater than fluid critical temperature')
    if Tmax > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmax cannot be greater than fluid critical temperature')
    Tmin = max(Tmin, cp.Props(Ref,'Tmin')+0.01)
    Tmax = min(Tmax, cp.Props(Ref,'Tcrit')-1e-5)
    Tsat = np.linspace(Tmin,Tmax,1000)
    (ssatL,ssatV)=(0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        ssatL[i] = cp.Props('S','T',Tsat[i],'Q',0,Ref)
        ssatV[i] = cp.Props('S','T',Tsat[i],'Q',1,Ref)
        
    ax.plot(ssatL,Tsat,'k')
    ax.plot(ssatV,Tsat,'k')
    
    if Tmax>cp.Props(Ref,'Tcrit')-2e-5:
        ax.plot(np.r_[ssatL[-1],ssatV[-1]],np.r_[Tsat[-1],Tsat[-1]],'k')
        ax.plot((ssatL[-1]+ssatV[-1])/2.,Tsat[-1],'o')

    ax.set_xlabel('Entropy [kJ/kg$\cdot$K]')
    ax.set_ylabel('Temperature [K]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()

def Ph(Ref, Tmin=None, Tmax = None, show = False, axis=None, **kwargs):
    
    """
    Make a pressure-enthalpy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    
    ax = axis if axis is not None else pylab.gca()
    if Tmin is None:
        Tmin = cp.Props(Ref,'Tmin')
    if Tmax is None:
        Tmax = cp.Props(Ref,'Tcrit')-1e-5
        
    if Tmin > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmin cannot be greater than fluid critical temperature')
    if Tmax > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmax cannot be greater than fluid critical temperature')
    Tmin = max(Tmin, cp.Props(Ref,'Tmin')+0.01)
    Tmax = min(Tmax, cp.Props(Ref,'Tcrit')-1e-5)
    
    Tsat = np.linspace(Tmin,Tmax,1000)
    (hsatL,psatL,hsatV,psatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        hsatL[i] = cp.Props('H','T',Tsat[i],'Q',0,Ref)
        hsatV[i] = cp.Props('H','T',Tsat[i],'Q',1,Ref)
        psatL[i] = cp.Props('P','T',Tsat[i],'Q',0,Ref)
        psatV[i] = cp.Props('P','T',Tsat[i],'Q',1,Ref)

    ax.plot(hsatL,psatL,'k')
    ax.plot(hsatV,psatV,'k')

    if Tmax>cp.Props(Ref,'Tcrit')-2e-5:
        ax.plot(np.r_[hsatL[-1],hsatV[-1]],np.r_[psatL[-1],psatV[-1]],'k')
        ax.plot((hsatL[-1]+hsatV[-1])/2.,(psatL[-1]+psatV[-1])/2.,'o')
        #ax.plot(np.r_[hsatL[-1],hsatV[-1]],np.r_[psatL[-1],psatV[-1]],'k')
    
    ax.set_xlabel('Enthalpy [kJ/kg]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()
    
def hs(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    
    """
    Make a enthalpy-entropy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    ax = axis if axis is not None else pylab.gca()
    if Tmin is None:
        Tmin = cp.Props(Ref,'Tmin')
    if Tmax is None:
        Tmax = cp.Props(Ref,'Tcrit')-1e-5
        
    if Tmin > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmin cannot be greater than fluid critical temperature')
    if Tmax > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmax cannot be greater than fluid critical temperature')
        
    Tmin = max(Tmin, cp.Props(Ref,'Tmin')+0.01)
    Tmax = min(Tmax, cp.Props(Ref,'Tcrit')-1e-5)
    
    Tsat = np.linspace(Tmin,Tmax,1000)
    (ssatL,hsatL,ssatV,hsatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        ssatL[i] = cp.Props('S','T',Tsat[i],'Q',0,Ref)
        ssatV[i] = cp.Props('S','T',Tsat[i],'Q',1,Ref)
        hsatL[i] = cp.Props('H','T',Tsat[i],'Q',0,Ref)
        hsatV[i] = cp.Props('H','T',Tsat[i],'Q',1,Ref)

    ax.plot(ssatL,hsatL,'k')
    ax.plot(ssatV,hsatV,'k')

    if Tmax>cp.Props(Ref,'Tcrit')-2e-5:
        ax.plot(np.r_[ssatL[-1],ssatV[-1]],np.r_[hsatL[-1],hsatV[-1]],'k')
        ax.plot((ssatL[-1]+ssatV[-1])/2.,(hsatL[-1]+hsatV[-1])/2.,'o')
        #ax.plot(ssatL[-1],hsatL[-1],'o')
        #ax.plot(np.r_[ssatL[-1],ssatV[-1]],np.r_[hsatL[-1],hsatV[-1]],'k')
    
    ax.set_xlabel('Entropy [kJ/kg/K]')
    ax.set_ylabel('Enthalpy [kJ/kg]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()
        
def Ps(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    
    """
    Make a pressure-entropy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    ax = axis if axis is not None else pylab.gca()
    if Tmin is None:
        Tmin = cp.Props(Ref,'Tmin')
    if Tmax is None:
        Tmax = cp.Props(Ref,'Tcrit')-1e-5
        
    if Tmin > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmin cannot be greater than fluid critical temperature')
    if Tmax > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmax cannot be greater than fluid critical temperature')
        
    Tmin = max(Tmin, cp.Props(Ref,'Tmin')+0.01)
    Tmax = min(Tmax, cp.Props(Ref,'Tcrit')-1e-5)
    
    Tsat = np.linspace(Tmin,Tmax,1000)
    (ssatL,psatL,ssatV,psatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        ssatL[i] = cp.Props('S','T',Tsat[i],'Q',0,Ref)
        ssatV[i] = cp.Props('S','T',Tsat[i],'Q',1,Ref)
        psatL[i] = cp.Props('P','T',Tsat[i],'Q',0,Ref)
        psatV[i] = cp.Props('P','T',Tsat[i],'Q',1,Ref)

    ax.plot(ssatL,psatL,'k')
    ax.plot(ssatV,psatV,'k')

    if Tmax>cp.Props(Ref,'Tcrit')-2e-5:
        ax.plot(np.r_[ssatL[-1],ssatV[-1]],np.r_[psatL[-1],psatV[-1]],'k')
        ax.plot((ssatL[-1]+ssatV[-1])/2.,(psatL[-1]+psatV[-1])/2.,'o')
        #ax.plot(np.r_[ssatL[-1],ssatV[-1]],np.r_[psatL[-1],psatV[-1]],'k')
    
    ax.set_xlabel('Entropy [kJ/kg/K]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()
        
def Prho(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    
    """
    Make a pressure-density plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    ax = axis if axis is not None else pylab.gca()
    if Tmin is None:
        Tmin = cp.Props(Ref,'Tmin')
    if Tmax is None:
        Tmax = cp.Props(Ref,'Tcrit')-1e-5
        
    if Tmin > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmin cannot be greater than fluid critical temperature')
    if Tmax > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmax cannot be greater than fluid critical temperature')
    Tmin = max(Tmin, cp.Props(Ref,'Tmin')+0.01)
    Tmax = min(Tmax, cp.Props(Ref,'Tcrit')-1e-5)
    
    Tsat = np.linspace(Tmin,Tmax,1000)
    (rhosatL,psatL,rhosatV,psatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        rhosatL[i] = cp.Props('D','T',Tsat[i],'Q',0,Ref)
        rhosatV[i] = cp.Props('D','T',Tsat[i],'Q',1,Ref)
        psatL[i] = cp.Props('P','T',Tsat[i],'Q',0,Ref)
        psatV[i] = cp.Props('P','T',Tsat[i],'Q',1,Ref)

    ax.plot(rhosatL,psatL,'k')
    ax.plot(rhosatV,psatV,'k')

    if Tmax>cp.Props(Ref,'Tcrit')-2e-5:
        ax.plot(np.r_[rhosatL[-1],rhosatV[-1]],np.r_[psatL[-1],psatV[-1]],'k')
        ax.plot((rhosatL[-1]+rhosatV[-1])/2.,(psatL[-1]+psatV[-1])/2.,'o')
        #ax.plot(np.r_[rhosatL[-1],rhosatV[-1]],np.r_[psatL[-1],psatV[-1]],'k')
    
    ax.set_xlabel('Density [kg/m$^3$]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()
        
def Trho(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    
    """
    Make a temperature-density plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    ax = axis if axis is not None else pylab.gca()
    if Tmin is None:
        Tmin = cp.Props(Ref,'Tmin')
    if Tmax is None:
        Tmax = cp.Props(Ref,'Tcrit')-1e-5
        
    if Tmin > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmin cannot be greater than fluid critical temperature')
    if Tmax > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmax cannot be greater than fluid critical temperature')
    Tmin = max(Tmin, cp.Props(Ref,'Tmin')+0.01)
    Tmax = min(Tmax, cp.Props(Ref,'Tcrit')-1e-5)
    
    Tsat = np.linspace(Tmin,Tmax,1000)
    (rhosatL,rhosatV)=(0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        rhosatL[i] = cp.Props('D','T',Tsat[i],'Q',0,Ref)
        rhosatV[i] = cp.Props('D','T',Tsat[i],'Q',1,Ref)

    ax.plot(rhosatL,Tsat,'k')
    ax.plot(rhosatV,Tsat,'k')

    if Tmax>cp.Props(Ref,'Tcrit')-2e-5:
        ax.plot(np.r_[rhosatL[-1],rhosatV[-1]],np.r_[Tsat[-1],Tsat[-1]],'k')
        ax.plot((rhosatL[-1]+rhosatV[-1])/2.,Tsat[-1],'o')
        #ax.plot(np.r_[rhosatL[-1],rhosatV[-1]],np.r_[Tsat[-1],Tsat[-1]],'k')
    
    ax.set_xlabel('Density [kg/m$^3$]')
    ax.set_ylabel('Temperature [K]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()

def PT(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    
    """
    Make a pressure-temperature plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    ax = axis if axis is not None else pylab.gca()
    if Tmin is None:
        Tmin = cp.Props(Ref,'Tmin')
    if Tmax is None:
        Tmax = cp.Props(Ref,'Tcrit')-1e-5
        
    if Tmin > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmin cannot be greater than fluid critical temperature')
    if Tmax > cp.Props(Ref,'Tcrit'):
        raise ValueError('Tmax cannot be greater than fluid critical temperature')
    Tmin = max(Tmin, cp.Props(Ref,'Tmin')+0.01)
    Tmax = min(Tmax, cp.Props(Ref,'Tcrit')-1e-5)
    
    Tsat = np.linspace(Tmin,Tmax,1000)
    (psatL,psatV)=(0.0*Tsat,0.0*Tsat)
    for i in range(len(Tsat)):
        psatL[i] = cp.Props('P','T',Tsat[i],'Q',0,Ref)
        psatV[i] = cp.Props('P','T',Tsat[i],'Q',1,Ref)

    ax.plot(Tsat,psatL,'k')
    ax.plot(Tsat,psatV,'k')

    if Tmax>cp.Props(Ref,'Tcrit')-2e-5:
        ax.plot(np.r_[Tsat[-1],Tsat[-1]],np.r_[psatL[-1],psatV[-1]],'k')
        ax.plot(Tsat[-1],(psatL[-1]+psatV[-1])/2.,'o')
        #ax.plot(np.r_[Tsat[-1],Tsat[-1]],np.r_[psatL[-1],psatV[-1]],'k')
    
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()

        
if __name__=='__main__':
    hs('R245fa', show = True)
    PT('R245fa', show = True)
    Ph('Helium', show = True)
    Trho('R245fa', show = True)
    Prho('R245fa', show = True)
    Ps('R290', show = True)
    Ps('R290', show = True)
    Ph('R290', show = True)
    Ts('R290', show = True)
    
