
import numpy, matplotlib, matplotlib.pyplot, math, re
from scipy.interpolate import interp1d

import CoolProp.CoolProp as CP


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
        axis=matplotlib.pyplot.gca()
    
    if fig is None:
        fig=matplotlib.pyplot.gcf()
    
    
    
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
        rot = numpy.arctan(dy_dx)/numpy.pi*180.
        
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
        rot = numpy.arctan(dy_dx)/numpy.pi*180.
        
    (x,y)=ToDataCoords(x,y,axis,fig)
    return (x,y,rot)

def show():
    """
    A convenience function to call pylab.show()
    """
    matplotlib.pyplot.show()
    

def drawIsoLines(Ref, plot, which, iValues=[], num=0, axis=None, fig=None):
    """
    Draw lines with constant values of type 'which' in terms of x and y as 
    defined by 'plot'. 'iMin' and 'iMax' are minimum and maximum value between
    which 'num' get drawn. 
    
    There should also be helpful error messages...
    """
    
    if axis is None:
        axis=matplotlib.pyplot.gca()
        
    if fig is None:
        fig=matplotlib.pyplot.gcf()
    
    if not plot is None:
        if not which is None:
            if not which=='all':
                lines = getIsoLines(Ref, plot, which, iValues=iValues, num=num, axis=axis)
                return drawLines(Ref,lines,axis)
            else:
                # TODO: assign limits to values automatically
                raise ValueError('Plotting all lines automatically is not supported, yet..')
            
                ll = _getIsoLineIds(plot)    
                if not len(ll)==len(iValues):
                    raise ValueError('Please provide a properly sized array of bounds.')
                for c,l in enumerate(ll):
                    drawIsoLines(Ref, plot, l, iValues=iValues[c], num=num, axis=axis, fig=fig)
    
    
def drawLines(Ref,lines,axis,plt_kwargs=None):
    """
    Just an internal method to systematically plot values from
    the generated 'line' dicts, method is able to cover the whole
    saturation curve. Closes the gap at the critical point and
    adds a marker between the two last points of bubble and 
    dew line if they reach up to critical point.
    
    Returns the an array of line objects that can be used to change 
    the colour or style afterwards.  
    """
    if not plt_kwargs is None:
        for line in lines:
            line['opts'] = plt_kwargs
    plottedLines = []
    if len(lines)==2 and (
      'bubble' in str(lines[0]['label']).lower() 
      and 'dew' in str(lines[1]['label']).lower()): 
        # We plot the saturation curve
        bubble = lines[0]
        dew    = lines[1]
        line, = axis.plot(bubble['x'],bubble['y'],**bubble['opts'])
        plottedLines.extend([line])
        line, = axis.plot(dew['x'],   dew['y'],   **dew['opts'])
        plottedLines.extend([line])
        # Do we need to test if this is T or p?
        Tmax = min(bubble['kmax'],dew['kmax'])
        if Tmax>CP.Props(Ref,'Tcrit')-2e-5:
            axis.plot(numpy.r_[bubble['x'][-1],dew['x'][-1]],numpy.r_[bubble['y'][-1],dew['y'][-1]],**bubble['opts'])
            #axis.plot((bubble['x'][-1]+dew['x'][-1])/2.,(bubble['y'][-1]+dew['y'][-1])/2.,'o',color='Tomato')
    else:
        for line in lines:
            line, = axis.plot(line['x'],line['y'],**line['opts'])
            plottedLines.extend([line])

    return plottedLines 


def getIsoLines(Ref, plot, iName, iValues=[], num=0, axis=None):
    """
    This is the core method to obtain lines in the dimensions defined
    by 'plot' that describe the behaviour of fluid 'Ref'. The constant 
    value is determined by 'iName' and has the values of 'iValues'. 
    
    'iValues' is an array-like object holding at least one element. Lines 
    are calculated for every entry in 'iValues'. If the input 'num' is 
    larger than the amount of entries in 'iValues', an internally defined 
    pattern is used to calculate an appropriate line spacing between the maximum
    and minimum values provided in 'iValues'. 
    
    Returns lines[num] - an array of dicts containing 'x' and 'y' 
    coordinates for bubble and dew line. Additionally, the dict holds 
    the keys 'label' and 'opts', those can be used for plotting as well.
    """
    if axis is None:
        axis=matplotlib.pyplot.gca()
        
    which = False
    if not plot is None:
        xName,yName,plot = _plotXY(plot)
        if not iName is None:
            if (iName in _getIsoLineIds(plot)) or (iName=='Q'):
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
        raise ValueError('This kind of isoline is not supported for a '+str(plot)+'-plot. Please choose from '+str(_getIsoLineIds(plot))+' or Q.')
     
    # Enforce the bounds!
    ((Axmin,Axmax), (Aymin,Aymax)) = _setBounds(Ref, plot, axis=axis)
    switchXY = False
#    # Problematic inputs that cannot be handled by the 
#    # internal CoolProp solver have to be avoided. 
#    # Converting everything to TD values:
#    rho_0 =  CP.Props(Ref,'rhocrit')
#    T_0   = (CP.Props(Ref,'Tcrit')+CP.Props(Ref,'Tmin'))/2. 
#    # Promote to lists for the 4 corners of the plot
#    rho_c = [rho_0, rho_0, rho_0, rho_0]
#    T_c   = [T_0  , T_0  , T_0  , T_0]
#    # find s from T and D inputs, solve for D
#    bounds_rho = (0.0                 , 1e10)
#    bounds_T   = (CP.Props(Ref,'Tmin'), 1e10)
    if plot=='TS':
#        # SciPy for scalar functions
#        from scipy.optimize._minimize import minimize_scalar        
#        # residual from known x and y to find i 
#        def f(i,x,y): return numpy.power((CP.Props('S','T',y,'D',i,Ref)-x),2)
#        # first case, upper right corner
#        res = minimize_scalar(f, rho_0, bounds=bounds_rho, args=(Axmax, Aymax), method='bounded')
#        T_c[0]   = Aymax
#        rho_c[0] = res.x
#        # second case, upper left corner
#        res = minimize_scalar(f, rho_0, bounds=bounds_rho, args=(Axmin, Aymax), method='bounded')
#        T_c[1]   = Aymax
#        rho_c[1] = res.x
#        # third case, lower left corner
#        res = minimize_scalar(f, rho_0, bounds=bounds_rho, args=(Axmin, Aymin), method='bounded')
#        T_c[2]   = Aymin
#        rho_c[2] = res.x
#        # fourth case, lower right corner
#        res = minimize_scalar(f, rho_0, bounds=bounds_rho, args=(Axmax, Aymin), method='bounded')
#        T_c[3]   = Aymin
#        rho_c[3] = res.x
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
    
    if switchXY:
        xmin = Aymin
        xmax = Aymax
        tmpName = yName
        yName = xName
        xName = tmpName
    else:
        xmin = Axmin
        xmax = Axmax

    # Calculate the points
    x0 = numpy.linspace(xmin,xmax,1000.)

    patterns = {
      'P' : (lambda x: numpy.logspace(math.log(x[0],2.), math.log(x[1],2.), num=x[2], base=2.)),
      'D' : (lambda x: numpy.logspace(math.log(x[0],2.), math.log(x[1],2.), num=x[2], base=2.)),
      'H' : (lambda x: numpy.linspace(x[0], x[1], num=x[2])),
      'T' : (lambda x: numpy.linspace(x[0], x[1], num=x[2])),
      'S' : (lambda x: numpy.linspace(x[0], x[1], num=x[2])),
      'Q' : (lambda x: numpy.linspace(x[0], x[1], num=x[2]))
    }
    
    # Get isoline spacing while honouring the inputs
    if len(iValues)<1: # No values given, use plot boundaries to determine limits
        raise ValueError('Automatic interval detection for isoline boundaries is not supported yet, use the iValues=[min, max] parameter.')
        iVal = [CP.Props(iName,'T',T_c[i],'D',rho_c[i],Ref) for i in range(len(T_c))]
        iVal = patterns[iName]([numpy.min(iVal),numpy.max(iVal),num]) 
    elif 1<len(iValues)<num: # We need more numbers
        iNum = num - len(iValues) 
        iVal = patterns[iName]([numpy.min(iValues),numpy.max(iValues),iNum+2]) # include min and max
        iVal = numpy.append(iVal[1:-1],iValues) # exclude min and max again
    else: # Either len(iValues)==1 or len(iValues)>=num
        iVal = numpy.array(iValues)
         
    iVal = numpy.sort(iVal) # sort again
        
    
    if iName=='Q':
        lines = _getSatLines(Ref, plot, x=iVal)
        return lines
    
    # TODO: Determine saturation state if two phase region present
    xVal = [x0 for i in iVal]
    
    (X,Y) = _getI_YX(Ref,iName,xName,yName,iVal,xVal)
    
    if switchXY:
        tmpY = Y
        Y    = X
        X    = tmpY
        
    lines = []
    for j in range(len(X)):
        line = {
          'x' : X[j],
          'y' : Y[j],
          'label' : _getIsoLineLabel(iName,iVal[j]),
          'opts': { 'color':getIsoLineColour(iName), 'lw':0.75, 'alpha':0.5 }
          }
        lines.append(line)
        
    return lines
    
def _satBounds(Ref,kind,xmin=None,xmax=None):
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
        xmin = CP.Props(Ref,str(minKey)     ) + 1e-5
    if xmax is None:
        xmax = CP.Props(Ref,str(kind)+'crit') - 1e-5
        
    if xmin > CP.Props(Ref,str(kind)+'crit'):
        raise ValueError('Minimum '+str(name)+' cannot be greater than fluid critical '+str(name)+'.')
    if xmax > CP.Props(Ref,str(kind)+'crit'):
        raise ValueError('Maximum '+str(name)+' cannot be greater than fluid critical '+str(name)+'.')
    
    xmin = max(xmin, CP.Props(Ref,str(minKey)    )  + 1e-5)
    xmax = min(xmax, CP.Props(Ref,str(kind)+'crit') - 1e-5)
    
    return (xmin,xmax)


def _setBounds(Ref, plot, axis=None):
    """
    Generates limits for the axes in terms of x,y defined by 'plot'
    based on temperature and pressure.  
    
    Returns a tuple containing ((xmin,xmax), (ymin,ymax))    
    """
    
    xName,yName,plot = _plotXY(plot)
    
    # Get current axis limits, be sure to set those before drawing isolines
    # if no limits are set, use triple point and critical conditions
    X = []
    X.append(CP.Props(xName,'T',1.5*CP.Props(Ref,'Tcrit'), 'P',   CP.Props(Ref,'ptriple'),Ref))
    X.append(CP.Props(xName,'T',1.1*CP.Props(Ref,'Tmin') , 'P', 1.5*CP.Props(Ref,'pcrit'),Ref))
    X.append(CP.Props(xName,'T',1.5*CP.Props(Ref,'Tcrit'), 'P', 1.5*CP.Props(Ref,'pcrit'),Ref))
    X.append(CP.Props(xName,'T',1.1*CP.Props(Ref,'Tmin') , 'P',   CP.Props(Ref,'ptriple'),Ref))
    Y = []
    Y.append(CP.Props(yName,'T',1.5*CP.Props(Ref,'Tcrit'), 'P',   CP.Props(Ref,'ptriple'),Ref))
    Y.append(CP.Props(yName,'T',1.1*CP.Props(Ref,'Tmin') , 'P', 1.5*CP.Props(Ref,'pcrit'),Ref))
    Y.append(CP.Props(yName,'T',1.1*CP.Props(Ref,'Tcrit'), 'P', 1.5*CP.Props(Ref,'pcrit'),Ref))
    Y.append(CP.Props(yName,'T',1.5*CP.Props(Ref,'Tmin') , 'P',   CP.Props(Ref,'ptriple'),Ref))
    
    minX = numpy.min(X)
    maxX = numpy.max(X)
    minY = numpy.min(Y)
    maxY = numpy.max(Y)
    
    if axis.get_autoscalex_on():
        axis.set_xlim([minX,maxX])
    else:
        [cuiX,cuaX] = axis.get_xlim()
        axis.set_xlim(left=numpy.max([minX,cuiX]))
        #axis.set_xlim(right=numpy.min(maxX,cuaX))
        
    if axis.get_autoscaley_on():
        axis.set_ylim([minY,maxY])
    else:
        [cuiY,cuaY] = axis.get_ylim()
        axis.set_ylim(bottom=numpy.max([minY,cuiY]))
        #axis.set_ylim(right=numpy.min(maxY,cuaY))
    
    [cuiX,cuaX] = axis.get_xlim()
    [cuiY,cuaY] = axis.get_ylim()
    
    return ((cuiX,cuaX),(cuiY,cuaY))


def _getSatLines(Ref, plot, kind=None, kmin=None, kmax=None, x=[0.,1.]):
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
    
    (kmin,kmax) = _satBounds(Ref, kind, xmin=kmin, xmax=kmax)   
    k0          = numpy.linspace(kmin,kmax,1000)
    
    iName       = 'Q'
    iVal        = x
    kVal        = [k0 for i in iVal]
    
    if xName!=kind:
        (Xx,Yx) = _getI_YX(Ref,iName,kind,xName,iVal,kVal)
    else:
        (Xx,Yx) = (kVal,kVal)
    
    if yName!=kind:
        (Xy,Yy) = _getI_YX(Ref,iName,kind,yName,iVal,kVal)
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
            line['opts'] = { 'color':getIsoLineColour(iName), 'lw':1.00 }
        elif iVal[j]==1.:
            line['label'] = 'dew line'
            line['opts'] = { 'color':getIsoLineColour(iName), 'lw':1.00 }
        else:
            line['label'] = _getIsoLineLabel(iName,iVal[j]),
            line['opts'] = { 'color':getIsoLineColour(iName), 'lw':0.75, 'alpha':0.5}
        
        lines.append(line)
        
    return lines


def _getI_YX(Ref,iName,xName,yName,iVal,xVal):
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
        Y = numpy.array([CP.Props(yName,iName,[jiVal],xName,ijxVal,Ref) for ijxVal in jxVal])
        y.append(Y)
        x.append(jxVal)
        
    return x,y


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


def _getIsoLineIds(plot):
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
#    colourMap = {                 
#      'T' : 'red',
#      'P' : 'cyan',
#      'H' : 'green',
#      'D' : 'blue',
#      'S' : 'yellow',
#      'Q' : 'black'
#    }
    colourMap = {                 
      'T' : 'Darkred',
      'P' : 'DarkCyan',
      'H' : 'DarkGreen',
      'D' : 'DarkBlue',
      'S' : 'DarkOrange',
      'Q' : 'black'
    }
    return colourMap[str(which)]


def _getIsoLineLabel(which,num):
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

    
def Ts(Ref,Tmin=None, Tmax=None, show=False, axis=None, **kwargs):
    """
    Make a temperature-entropy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """

    ax = axis if axis is not None else matplotlib.pyplot.gca()
    lines = _getSatLines(Ref, 'Ts', kind='T', kmin=Tmin, kmax=Tmax)  
    drawLines(Ref,lines,ax)
    # or alternatively:
    #drawIsoLines(Ref, 0, 1, which='Q', plot='Ts', axis=ax)
    
    ax.set_xlabel('Entropy [kJ/kg$\cdot$K]')
    ax.set_ylabel('Temperature [K]')
    ax.autoscale(enable=True)
    if show:
        matplotlib.pyplot.show()
    return ax


def Ph(Ref, Tmin=None, Tmax = None, show = False, axis=None, **kwargs):
    """
    Make a pressure-enthalpy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    
    ax = axis if axis is not None else matplotlib.pyplot.gca()
    lines  = _getSatLines(Ref, 'ph', kind='T', kmin=Tmin, kmax=Tmax)  
    drawLines(Ref,lines,ax)
    
    ax.set_xlabel('Enthalpy [kJ/kg]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        matplotlib.pyplot.show()
    return ax
    
    
def hs(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    """
    Make a enthalpy-entropy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    
    ax = axis if axis is not None else matplotlib.pyplot.gca()
    lines  = _getSatLines(Ref, 'hs', kind='T', kmin=Tmin, kmax=Tmax)  
    drawLines(Ref,lines,ax)
    
    ax.set_xlabel('Entropy [kJ/kg/K]')
    ax.set_ylabel('Enthalpy [kJ/kg]')
    ax.autoscale(enable=True)
    if show:
        matplotlib.pyplot.show()
    return ax
        
        
def Ps(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    """
    Make a pressure-entropy plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    
    ax = axis if axis is not None else matplotlib.pyplot.gca()
    lines  = _getSatLines(Ref, 'ps', kind='T', kmin=Tmin, kmax=Tmax)  
    drawLines(Ref,lines,ax)
    
    ax.set_xlabel('Entropy [kJ/kg/K]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        matplotlib.pyplot.show()
    return ax
        
        
def Prho(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    """
    Make a pressure-density plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    
    ax = axis if axis is not None else matplotlib.pyplot.gca()
    lines  = _getSatLines(Ref, 'pd', kind='T', kmin=Tmin, kmax=Tmax)  
    drawLines(Ref,lines,ax)
    
    ax.set_xlabel('Density [kg/m$^3$]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        matplotlib.pyplot.show()
    return ax
        
        
def Trho(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    """
    Make a temperature-density plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    
    ax = axis if axis is not None else matplotlib.pyplot.gca()
    lines  = _getSatLines(Ref, 'Td', kind='T', kmin=Tmin, kmax=Tmax)  
    drawLines(Ref,lines,ax)
    
    ax.set_xlabel('Density [kg/m$^3$]')
    ax.set_ylabel('Temperature [K]')
    ax.autoscale(enable=True)
    if show:
        matplotlib.pyplot.show()
    return ax


def PT(Ref, Tmin=None, Tmax = None, show = False, axis = None, **kwargs):
    """
    Make a pressure-temperature plot for the given fluid
    
    Will plot in the current axis unless the optional parameter *axis* gives the name for the axis to use
    """
    
    ax = axis if axis is not None else matplotlib.pyplot.gca()
    lines  = _getSatLines(Ref, 'pT', kind='T', kmin=Tmin, kmax=Tmax)  
    drawLines(Ref,lines,ax)
    
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Pressure [kPa]')
    ax.autoscale(enable=True)
    if show:
        matplotlib.pyplot.show()
    return ax

        
if __name__=='__main__':
#    import matplotlib
#    matplotlib.use('Qt4Agg')
#    
#    import matplotlib.pyplot
#    
#    fluid = "n-Pentane"
#    # define custom styles
#    plt_kwargs_1 = {"color": "green","linewidth": 1.5}
#    plt_kwargs_2 = {"color": "red","linewidth": 0.75}
#    
#    raw_input("Press Enter to draw quality lines...")
#    fig, ((ax1, ax2)) = matplotlib.pyplot.subplots(1, 2, sharey='row')
#    drawIsoLines(fluid, 'Ts', 'Q', [0.0, 1.0], axis=ax1) # for predefined styles
#    drawIsoLines(fluid, 'Ts', 'Q', [0.3, 0.7], axis=ax1) # for predefined styles
#    # Get the data points
#    saturation = getIsoLines(fluid, 'Ts', 'Q', [0.0,  1.0 ], axis=ax2)
#    quality    = getIsoLines(fluid, 'Ts', 'Q', [0.3,  0.7 ], axis=ax2)
#    # and draw the lines
#    drawLines(fluid,saturation[:],ax2,plt_kwargs=plt_kwargs_1)
#    drawLines(fluid,quality[:],ax2,plt_kwargs=plt_kwargs_2)
#    matplotlib.pyplot.show()
#    raw_input("Press Enter to draw more lines...")
#    fig, (ax1) = matplotlib.pyplot.subplots(1, 1)
#    drawLines(fluid,saturation,ax1)
#    drawIsoLines(fluid, 'Ts', 'Q', [0.0, 1.0], axis=ax1) # for predefined styles
#    drawIsoLines(fluid, 'Ts', 'Q', [0.25, 0.75], num=3, axis=ax1) # for predefined styles
#    drawIsoLines(fluid, 'Ts', 'H', [50,  100], num=10, axis=ax1) # for predefined styles
#    drawIsoLines(fluid, 'Ts', 'P', [10,  5000], num=10, axis=ax1) # for predefined styles
#    matplotlib.pyplot.show()
#    hs('n-Pentane', show = True)
#    raw_input("Press Enter to continue...")
    hs('R245fa', show = True)
    raw_input("Press Enter to continue...")
    PT('R245fa', show = True)
    raw_input("Press Enter to continue...")
    Ph('Helium', show = True)
    raw_input("Press Enter to continue...")
    Trho('R245fa', show = True)
    raw_input("Press Enter to continue...")
    Prho('R245fa', show = True)
    raw_input("Press Enter to continue...")
    Ps('R290', show = True)
    raw_input("Press Enter to continue...")
    Ph('R290', show = True)
    raw_input("Press Enter to continue...")
    Ts('R290', show = True)
    raw_input("Press Enter to continue...")
