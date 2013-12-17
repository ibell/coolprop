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
    fluid_ref = 'R290'

    def Ts_plot_tests():
        from CoolProp.Plots import Ts
        Ts(fluid_ref, show=True)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Ts(fluid_ref, show=True, axis=ax)

        Ts(fluid_ref, show=True, Tmin=200, Tmax=300)

    def Ph_plot_tests():
        from CoolProp.Plots import Ph
        Ph(fluid_ref, show=True)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Ph(fluid_ref, show=True, axis=ax)

        Ph(fluid_ref, show=True, Tmin=200, Tmax=300)

    def PT_plot_tests():
        from CoolProp.Plots import PT
        PT(fluid_ref, show=True)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        PT(fluid_ref, show=True, axis=ax)

        PT(fluid_ref, show=True, Tmin=200, Tmax=300)

    def Ps_plot_tests():
        from CoolProp.Plots import Ps
        Ps(fluid_ref, show=True)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Ps(fluid_ref, show=True, axis=ax)

        Ps(fluid_ref, show=True, Tmin=200, Tmax=300)

    def Prho_plot_tests():
        from CoolProp.Plots import Prho
        Prho(fluid_ref, show=True)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Prho(fluid_ref, show=True, axis=ax)

        Prho(fluid_ref, show=True, Tmin=200, Tmax=300)

    def Trho_plot_tests():
        from CoolProp.Plots import Trho
        Trho(fluid_ref, show=True)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        Trho(fluid_ref, show=True, axis=ax)

        Trho(fluid_ref, show=True, Tmin=200, Tmax=300)

    def hs_plot_tests():
        from CoolProp.Plots import hs
        hs(fluid_ref, show=True)

        from matplotlib import pyplot
        fig = pyplot.figure(2)
        ax = fig.gca()
        hs(fluid_ref, show=True, axis=ax)

        hs(fluid_ref, show=True, Tmin=200, Tmax=300)

    def Isolines_plot_tests():
        from matplotlib import pyplot
        from CoolProp.Plots import Ts, drawIsoLines
        ax = Ts(fluid_ref)
        #ax.set_xlim([-0.5, 1.5])
        #ax.set_ylim([300, 530])
        quality = drawIsoLines(fluid_ref, 'Ts', 'Q', [0.3, 0.5, 0.7, 0.8], axis=ax)
        isobars = drawIsoLines(fluid_ref, 'Ts', 'P', [100, 2000], num=5, axis=ax)
        isochores = drawIsoLines(fluid_ref, 'Ts', 'D', [2, 600], num=7, axis=ax)
        pyplot.show()

    Ts_plot_tests()
    Ph_plot_tests()
    Ps_plot_tests()
    PT_plot_tests()
    Prho_plot_tests()
    Trho_plot_tests()
    hs_plot_tests()
    Isolines_plot_tests()


def test_new_code():
    fluid_ref = 'Water'

    def Ts_plot_tests():
        from CoolProp.Plots import PropsPlot
        plt = PropsPlot(fluid_ref, 'Ts')
        plt.show()

    def Ph_plot_tests():
        from CoolProp.Plots import PropsPlot
        plt = PropsPlot(fluid_ref, 'Ph')
        plt.show()

    def Isolines_plot_tests():
        from CoolProp.Plots import PropsPlot
        plt = PropsPlot(fluid_ref, 'Ts')
        #plt.set_axis_limits([-0.5, 1.5, 300, 530])
        plt.draw_isolines('Q', [0.3, 0.5, 0.7, 0.8])
        plt.draw_isolines('P', [100, 2000], num=5)
        plt.draw_isolines('D', [2, 600], num=7)
        plt.show()

    def Graph_annotations():
        from CoolProp.Plots import PropsPlot, IsoLines
        plt = PropsPlot(fluid_ref, 'Ts')
        plt.draw_isolines('Q', [0.3, 0.5, 0.7, 0.8])
        plt.draw_isolines('P', [100, 2000], num=5)
        plt.draw_isolines('D', [2, 600], num=7)
        plt.title('New Title')
        plt.xlabel('New x label')
        plt.ylabel('New y label')
        plt.show()
        plt = IsoLines(fluid_ref, 'Ts', 'P')
        plt.draw_isolines([100, 2000], num=5)
        plt.show()

    def Mixture():
        from CoolProp.Plots import PropsPlot
        plt = PropsPlot('REFPROP-MIX:R32[0.47319469]&R125[0.2051091]&R134a[0.32169621]', 'TD')
        plt._plot_default_annotations()
        plt.show()

    Ts_plot_tests()
    Ph_plot_tests()
    Isolines_plot_tests()
    Graph_annotations()
    Mixture()


if __name__=='__main__':
    import nose
    nose.runmodule()
