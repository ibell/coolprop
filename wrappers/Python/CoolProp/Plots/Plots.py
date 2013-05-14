import numpy as np
import CoolProp.CoolProp as cp
from scipy.interpolate import interp1d
from types import NoneType
from matplotlib import pylab
import math
import re

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
    

def getI_YX(Ref,iName,xName,yName,iVal,xVal):
    """
    Calculates lines for constant iName (iVal) over an interval of xName (xVal). 
    
    Returns (x[],y[]) - a tuple of arrays containing the values in x and y dimensions.    
    """
    
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


def drawIsoLines(Ref, iMin, iMax, numberOfLines=5, which=None, plot=None, axis=None, fig=None):
    """
    Draw lines with constant values of type 'which' in terms of x and y as 
    defined by 'plot'. 'iMin' and 'iMax' are minimum and maximum value between
    which 'numberOfLines' get drawn. 
    
    There should also be helpful error messages...
    """
    
    if axis is None:
        axis=pylab.gca()
        
    if fig is None:
        fig=pylab.gcf()
    
    if not plot is None:
        if not which is None:
            if not which=='all':
                lines = getIsoLines(Ref, plot, which, iMin, iMax, numberOfLines, axis)
                _drawLines(Ref,lines,axis)
            else:
                # TODO: assign limits to values automatically
                raise ValueError('Plotting all lines automatically is not supported, yet..')
                for l in getIsoLineIds(plot):
                    drawIsoLines(Ref, iMin, iMax, which=l, plot=plot, axis=axis, fig=fig) 
    
    
def satBounds(Ref,kind,xmin=None,xmax=None):
    """
    Generates limits for the saturation line in either T or p determined by 'kind'. 
    If xmin or xmax are provided, values will be checked against the allowable 
    range for the EOS and an error might be generated. 
    
    Returns a tuple containing (xmin,xmax)
    """
    
    kind = str(kind).upper()
    name = ''
    minKey = ''
    if (str(kind)=='P'):
        kind = 'p'
        name = 'pressure'
        minKey = 'ptriple'
    elif (str(kind)=='T'):
        name = 'temperature'
        minKey = 'Tmin'
    else:
        raise ValueError('Saturation curves can only be computed for T or p.')
    
    if xmin is None:
        xmin = cp.Props(Ref,str(minKey)     ) + 1e-5
    if xmax is None:
        xmax = cp.Props(Ref,str(kind)+'crit') - 1e-5
        
    if xmin > cp.Props(Ref,str(kind)+'crit'):
        raise ValueError('Minimum '+str(name)+' cannot be greater than fluid critical '+str(name)+'.')
    if xmax > cp.Props(Ref,str(kind)+'crit'):
        raise ValueError('Maximum '+str(name)+' cannot be greater than fluid critical '+str(name)+'.')
    
    xmin = max(xmin, cp.Props(Ref,str(minKey)    )  + 1e-5)
    xmax = min(xmax, cp.Props(Ref,str(kind)+'crit') - 1e-5)
    
    return (xmin,xmax)


def getSatLines(Ref, plot, kind=None, kmin=None, kmax=None, x=[0.,1.]):
    """
    Calculates bubble and dew line in the quantities for your plot. 
    
    You can specify if you need evenly spaced entries in either 
    pressure or temperature by supplying kind='p' and kind='T' 
    (default), respectively.
    
    Limits can be set with kmin (default: minimum from EOS) and 
    kmax (default: critical value). 
    
    Returns lines[] - a 2D array of dicts containing 'x' and 'y' coordinates 
    for bubble and dew line. Additionally, the dict holds the keys 
    'kmax', 'label' and 'opts', those can be used for plotting as well.
    """
    
    xName,yName,plot = _plotXY(plot)
    
    if (str(kind).lower()=='p'):
        kind = 'P'
    else:
        kind = 'T' 
    
    (kmin,kmax) = satBounds(Ref, kind, xmin=kmin, xmax=kmax)   
    k0          = np.linspace(kmin,kmax,1000)
    
    iName       = 'Q'
    iVal        = x
    kVal        = [k0 for i in iVal]
    
    if xName!=kind:
        (Xx,Yx) = getI_YX(Ref,iName,kind,xName,iVal,kVal)
    else:
        (Xx,Yx) = (kVal,kVal)
    
    if yName!=kind:
        (Xy,Yy) = getI_YX(Ref,iName,kind,yName,iVal,kVal)
    else:
        (Xy,Yy) = (kVal,kVal)
    
    # Merge the two lines, capital Y holds important information. We merge on X values
    # Every entry, eg. Xy, contains two arrays of values.   
    lines = []
    for j in range(len(Yx)): # two dimensions: i = {0,1} 
        line = {
          'x' : Yx[j],
          'y' : Yy[j], 
          'kmax' : kmax
          }
        if iVal[j]==0.:
            line['label'] = 'bubble line'
            line['opts'] = { 'color':getIsoLineColour(iName), 'lw':1 }
        elif iVal[j]==1.:
            line['label'] = 'dew line'
            line['opts'] = { 'color':getIsoLineColour(iName), 'lw':1 }
        else:
            line['label'] = getIsoLineLabel(iName,iVal[j]),
            line['opts'] = { 'color':getIsoLineColour(iName), 'lw':1, 'alpha':0.25}
        
        lines.append(line)
        
    return lines

  
def getIsoLines(Ref, plot, iName, iMin, iMax, numberOfLines, axis):
    """
    Calculates lines of constant 'iName' between the values for 'iMin'
    and 'iMax'. The input 'plot' determines what will be returned 
    as x and y coordinates.
    
    Returns lines[numberOfLines] - an array of dicts containing 
    'x' and 'y' coordinates for bubble and dew line. Additionally, 
    the dict holds the keys 'label' and 'opts', those can be used 
    for plotting as well.
    """
    which = False
    if not plot is None:
        xName,yName,plot = _plotXY(plot)
        if not iName is None:
            if (iName in getIsoLineIds(plot)) or (iName=='Q'):
                which = True            
            else:
                which=False
        else:
            which=False
    else: 
        plot = False
        
    if not plot:
        raise ValueError('You have to specify the kind of plot.')
    
    #xName,yName,plot = _plotXY(plot)
    if not which:
        raise ValueError('This kind of isoline is not supported for a '+str(plot)+'-plot. Please choose from '+str(getIsoLineIds(plot))+' or Q.')
     
    # Problematic inputs that cannot be handled by the 
    # internal CoolProp solver have to be avoided. 
    # Converting everything to TD values:
    switchXY = False
    if plot=='TS':
        if iName=='D':
            switchXY = True # TD is defined, SD is not
    elif plot=='PH':
        if iName=='S':
            switchXY = True # PS is more stable than HS
        if iName=='T':
            switchXY = True # PT is defined, HT is not
        if iName=='D':
            switchXY = True # PD is defined, HD is not
    elif plot=='HS':
        if iName=='T':
            raise ValueError('You should not reach this point!')
        if iName=='D':
            raise ValueError('You should not reach this point!')
    elif plot=='PS':
        if iName=='T':
            switchXY = True # PT is defined, ST is not
        if iName=='D':
            switchXY = True # PD is defined, SD is not
    elif plot=='PD':
        if iName=='S':
            switchXY = True # PS is defined, DS is not
        if iName=='H':
            switchXY = True # PH is defined, DH is not
    elif plot=='TD':
        if iName=='S':
            raise ValueError('You should not reach this point!')
        if iName=='H':
            raise ValueError('You should not reach this point!')
    elif plot=='PT':
        if iName=='S':
            switchXY = True # PS is defined, TS is not

    # Get current axis limits, be sure to set those before drawing isolines
    if switchXY:
        [Axmin,Axmax]=axis.get_ylim()
        #[Aymin,Aymax]=axis.get_xlim()
        tmpName = yName
        yName = xName
        xName = tmpName
    else:
        [Axmin,Axmax]=axis.get_xlim()
        #[Aymin,Aymax]=axis.get_ylim()    
    
    # Determine x range for plotting
    if xName=='T': #Sacrifice steps 
        Axmin = max(cp.Props(Ref,'Tmin'), Axmin)
        #Axmax = Axmax + 273.15 
    x0 = np.linspace(Axmin,Axmax,1000)

    patterns = {
      'P' : (lambda x: np.logspace(math.log(x[0],2), math.log(x[1],2), num=numberOfLines, base=2)),
      'D' : (lambda x: np.logspace(math.log(x[0],2), math.log(x[1],2), num=numberOfLines, base=2)),
      'H' : (lambda x: np.linspace(x[0], x[1], num=numberOfLines)),
      'T' : (lambda x: np.linspace(x[0], x[1], num=numberOfLines)),
      'S' : (lambda x: np.linspace(x[0], x[1], num=numberOfLines)),
      'Q' : (lambda x: np.linspace(x[0], x[1], num=numberOfLines))
    }
    
    # TODO: Use the y range to determine spacing of isolines    
#    iVal = [cp.Props(iName,yName,Aymin,xName,Axmin,Ref),
#            cp.Props(iName,yName,Aymax,xName,Axmin,Ref),
#            cp.Props(iName,yName,Aymin,xName,Axmax,Ref),
#            cp.Props(iName,yName,Aymax,xName,Axmax,Ref) ]
#    iVal = patterns[iName]([np.min(iVal),np.max(iVal)]) 
    iVal = patterns[iName]([iMin,iMax])       
                
    if iName=='Q':
        lines = getSatLines(Ref, plot, x=iVal)
        return lines
    
    # TODO: Determine saturation state if two phase region present
    xVal = [x0 for i in iVal]
    
    (X,Y) = getI_YX(Ref,iName,xName,yName,iVal,xVal)
    
    if switchXY:
        tmpY = Y
        Y    = X
        X    = tmpY
        
    lines = []
    for j in range(len(X)):
        line = {
          'x' : X[j],
          'y' : Y[j],
          'label' : getIsoLineLabel(iName,iVal[j]),
          'opts': { 'color':getIsoLineColour(iName), 'lw':1, 'alpha':0.25 }
          }
        lines.append(line)
        
    return lines


def _plotXY(plot):
    """
    Creates strings for the x and y-axis values from a given plot 
    name like Ts, Ph, hs, Ps, Prho, Trho and PT.
    """
    
    # check if plot has the right format    
    plots = ['TS','PH','HS','PS','PD','TD','PT']
    
    rhoPat = re.compile("rho", re.IGNORECASE)
    plot = rhoPat.sub("D", plot).upper()
    
    if not plot in plots:
        raise ValueError('You have to specify the kind of plot, use Ts, Ph, hs, Ps, Prho, Trho or PT.')
    
    # Get strings to feed to Props function
    yName = str(plot[0] )
    xName = str(plot[1:])
    
    return xName,yName,yName+xName


def getIsoLineIds(plot):
    # define which isolines are allowed for the different plots
#    plots = {
#      'TS' : ['P','H'],#,'D'],
#      'PH' : ['S'],#,'T','D'],
#      'HS' : ['P'],#,'T','D'],
#      'PS' : ['H'],#,'T','D'],
#      'PD' : ['T'],#,'S','H'],
#      'TD' : ['P'],#,'S','H'],
#      'PT' : ['D','P']#,'S']
#    }
    # Changing the XY coordinates allows for more
    # combinations of isolines.
    plots = {
      'TS' : ['P','D'],#,'H'],
      'PH' : ['S','T','D'],
      'HS' : ['P'],#,'T','D'],
      'PS' : ['H','T','D'],
      'PD' : ['T','S','H'],
      'TD' : ['P'],#,'S','H'],
      'PT' : ['D','P','S']
    }
    return plots[str(plot)]


def getIsoLineColour(which):
    colourMap = {                 
      'T' : 'red',
      'P' : 'cyan',
      'H' : 'green',
      'D' : 'blue',
      'S' : 'yellow',
      'Q' : 'black'
    }
    return colourMap[str(which)]


def getIsoLineLabel(which,num):
    labelMap = {                 
      'T' : [r'$T = ','$ K'],
      'P' : [r'$p = ','$ kPa'],
      'H' : [r'$h = ','$ kJ/kg'],
      'D' : [r'$\rho = ','$ kg/m$^3$'],
      'S' : [r'$s = ','$ kJ/kg-K'],
      'Q' : [r'$x = ','$']
    }
    l = labelMap[str(which)]
    val = l[0]+str(num)+l[1]
    return val


def _drawLines(Ref,lines,axis):
    """
    Just an internal method to systematically plot values from
    the generated 'line' dicts, is able to cover the whole
    saturation curve. Closes the gap at the critical point and
    adds a marker between the two last points of bubble and 
    dew line if they reach up to critical point.
    """
       
    if len(lines)==2 and (
      'bubble' in str(lines[0]['label']).lower() 
      and 'dew' in str(lines[1]['label']).lower()): 
        # We plot the saturation curve
        bubble = lines[0]
        dew    = lines[1]
        axis.plot(bubble['x'],bubble['y'],**bubble['opts'])
        axis.plot(dew['x'],   dew['y'],   **dew['opts'])
        # Do we need to test if this is T or p?
        Tmax = min(bubble['kmax'],dew['kmax'])
        if Tmax>cp.Props(Ref,'Tcrit')-2e-5:
            axis.plot(np.r_[bubble['x'][-1],dew['x'][-1]],np.r_[bubble['y'][-1],dew['y'][-1]],**bubble['opts'])
            axis.plot((bubble['x'][-1]+dew['x'][-1])/2.,(bubble['y'][-1]+dew['y'][-1])/2.,'o')
    else:
        for line in lines:
            axis.plot(line['x'],line['y'],**line['opts'])
        
    
def Ts(Ref,Tmin=None, Tmax=None, show=False, axis=None, **kwargs):
    """
    Make a temperature-entropy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """

    ax = axis if axis is not None else pylab.gca()
    lines = getSatLines(Ref, 'Ts', kind='T', kmin=Tmin, kmax=Tmax)  
    _drawLines(Ref,lines,ax)
    # or alternatively:
    #drawIsoLines(Ref, 0, 1, which='Q', plot='Ts', axis=ax)
    
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
    lines  = getSatLines(Ref, 'ph', kind='T', kmin=Tmin, kmax=Tmax)  
    _drawLines(Ref,lines,ax)
    
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
    lines  = getSatLines(Ref, 'hs', kind='T', kmin=Tmin, kmax=Tmax)  
    _drawLines(Ref,lines,ax)
    
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
    lines  = getSatLines(Ref, 'ps', kind='T', kmin=Tmin, kmax=Tmax)  
    _drawLines(Ref,lines,ax)
    
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
    lines  = getSatLines(Ref, 'pd', kind='T', kmin=Tmin, kmax=Tmax)  
    _drawLines(Ref,lines,ax)
    
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
    lines  = getSatLines(Ref, 'Td', kind='T', kmin=Tmin, kmax=Tmax)  
    _drawLines(Ref,lines,ax)
    
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
    lines  = getSatLines(Ref, 'pT', kind='T', kmin=Tmin, kmax=Tmax)  
    _drawLines(Ref,lines,ax)
    
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        pylab.show()

        
if __name__=='__main__':
#    hs('R245fa', show = True)
#    raw_input("Press Enter to continue...")
#    PT('R245fa', show = True)
#    raw_input("Press Enter to continue...")
#    Ph('Helium', show = True)
#    raw_input("Press Enter to continue...")
#    Trho('R245fa', show = True)
#    raw_input("Press Enter to continue...")
#    Prho('R245fa', show = True)
#    raw_input("Press Enter to continue...")
#    Ps('R290', show = True)
#    raw_input("Press Enter to continue...")
#    Ph('R290', show = True)
#    raw_input("Press Enter to continue...")
#    Ts('R290', show = True)
#    raw_input("Press Enter to continue...")

    ax = pylab.gca()
    Ref = 'n-Pentane'
    # Set limits, be sure to use the 
    # internal CoolProp Units!
    ax.set_xlim([-0.5,1.5])
    ax.set_ylim([300,530])
    
    lines = []
    lines.extend(getSatLines(Ref, 'Ts'))
    lines.extend(getIsoLines(Ref, 'Ts', 'Q', 0.6, 0.8, 3, ax))
    lines.extend(getIsoLines(Ref, 'Ts', 'P', 100, 2000, 5, ax))
    lines.extend(getIsoLines(Ref, 'Ts', 'D', 2,    600, 5, ax))
#    lines.extend(getIsoLines(Ref, 'Ts', 'H', 100,  300, 5, ax))
    for line in lines:
        ax.plot(line['x'],np.array(line['y'])-273.15,**line['opts'])
    # Adjust the T limits to Celsius
    ax.set_ylim([25,250])
    pylab.show()

#    # Example for getSatLines
#    ax = pylab.gca()
#    (x,y) = getSatLines('n-Pentane', 'Ts')
#    ax.plot(x[0],y[0],'k') # bubble line
#    ax.plot(x[1],y[1],'k') # dew line
##    ax.set_xlim([-0.5,1.5])
##    ax.set_ylim([300,450])
#    drawIsoLines('n-Pentane', 100,  100, numberOfLines=1, which='P', plot='Ts', axis=ax)
#    x = getIsoLines('n-Pentane', 'Ts', 'P', 200, 200, 1, ax)
##    (x,y) = getI_YX('n-Pentane','P','S','T',[200],[np.linspace(-0.5, 1.5)])
#    ax.plot(x[0],x[0],'blue') 
#    ax.plot(line['x'],numpy.array(line['y'])-273.15,color=line['colour'])
#    pylab.show()
    
    
