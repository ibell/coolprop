import numpy as np
import pylab
#from CoolProp.Plots.Plots import getIsoLines, drawIsoLines, Ts

#def test_Trho():
#    for Fluid in CoolProp.__fluids__:
#        for T in np.linspace(Props(Fluid,'Tmin')+1e-5,Props(Fluid,'Tcrit')+100,20):
#            for rho in np.linspace(1e-10,Props(Fluid,'rhocrit')*3,20):
#                yield check_rho,Fluid,T,rho
#
#def check_rho(Fluid,T,rho):
#    p = Props('P','T',T,'D',rho,Fluid)
#
#def test_namedPlots():
#    return True
#
#def test_functions():
#
#    Ref = 'n-Pentane'
#
#    # Version A, get lines and do the plotting here
#    fig, ax = pylab.subplots(1, 1)
#    # Set limits, be sure to use internal CoolProp Units!
#    ax.set_xlim([-0.5,1.5])
#    ax.set_ylim([300,530])
#    lines = []
#    lines.extend(getIsoLines(Ref, 'Ts', 'Q', [0.0, 0.6,  0.75, 0.775, 0.8, 1.0], axis=ax))
#    lines.extend(getIsoLines(Ref, 'Ts', 'P', [100, 2000], num=5, axis=ax))
#    lines.extend(getIsoLines(Ref, 'Ts', 'D', [2,    600], num=5, axis=ax))
##    lines.extend(getIsoLines(Ref, 'Ts', 'H', 100,  300, 5, ax))
#    for line in lines:
#        ax.plot(line['x'],np.array(line['y'])-273.15,**line['opts'])
#    # Adjust the T limits to Celsius
#    ax.set_ylim([25,250])
#
#    # Version B, use built-in drawing functions and receiver line objects
#    fig, ax = pylab.subplots(1, 1)
#    ax = Ts(Ref)
#    ax.set_xlim([-0.5,1.5])
#    ax.set_ylim([300,530])
#    quality    = drawIsoLines(Ref, 'Ts', 'Q', [0.3,  0.75, 0.775, 0.8], axis=ax)
#    isobars    = drawIsoLines(Ref, 'Ts', 'P', [100, 2000]             , num=5, axis=ax)
#    isochores  = drawIsoLines(Ref, 'Ts', 'D', [2,    600]             , num=7, axis=ax)
##    isenthalps = drawIsoLines(Ref, 'Ts', 'H', 100,  300, 5, ax)
#
#
#    #pylab.show()


def test_back_compatibility():
    def Ts_plot_tests():
        from CoolProp.Plots import Ts
        Ts('R290', show=True)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Ts('R290', show=True, axis=ax)

        Ts('R290', show=True, Tmin=200, Tmax=300)

    def Ph_plot_tests():
        from CoolProp.Plots import Ph
        Ph('R290', show=True)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Ph('R290', show=True, axis=ax)

        Ph('R290', show=True, Tmin=200, Tmax=300)

    def Isolines_plot_tests():
        from matplotlib import pyplot
        from CoolProp.Plots import Ts, drawIsoLines
        Ref = 'n-Pentane'
        ax = Ts(Ref)
        ax.set_xlim([-0.5, 1.5])
        ax.set_ylim([300, 530])
        quality = drawIsoLines(Ref, 'Ts', 'Q', [0.3, 0.5, 0.7, 0.8], axis=ax)
        isobars = drawIsoLines(Ref, 'Ts', 'P', [100, 2000], num=5, axis=ax)
        isochores = drawIsoLines(Ref, 'Ts', 'D', [2, 600], num=7, axis=ax)
        pyplot.show()

    #Ts_plot_tests()
    #Ph_plot_tests()
    Isolines_plot_tests()


def test_new_code():
    def Ts_plot_tests():
        from CoolProp.Plots import PropsPlot
        plt = PropsPlot('R290', 'Ts')
        plt.show()

    def Ph_plot_tests():
        from CoolProp.Plots import PropsPlot
        plt = PropsPlot('R290', 'Ph')
        plt.show()

    def Isolines_plot_tests():
        from CoolProp.Plots import PropsPlot
        plt = PropsPlot('n-Pentane', 'Ts')
        plt.set_axis_limits([-0.5, 1.5, 300, 530])
        plt.draw_isolines('Q', [0.3, 0.5, 0.7, 0.8])
        plt.draw_isolines('P', [100, 2000], num=5)
        plt.draw_isolines('D', [2, 600], num=7)
        plt.show()

    #Ts_plot_tests()
    #Ph_plot_tests()
    Isolines_plot_tests()


if __name__=='__main__':
    import nose
    nose.runmodule()