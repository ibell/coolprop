Python Plotting
===============

.. note::
    The plotting API has been changed as of v4.0
    Examples using the 

The following example can be used to create a Temperature-Entropy plot for
propane (R290):

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    plt = PropsPlot('R290', 'Ts')
    plt.show()


The following example can be used to create a Pressure-Enthalpy plot for R410A:

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    plt = PropsPlot('R410A', 'Ph')
    plt.show()

The available plots are:

== ====================
PT Pressure-Temperature
PD Pressure-Density
PH Pressure-Enthalpy
PS Pressure-Entropy
TD Temperature-Density
TS Temperatre-Entropy
HS Enthalpy-Entropy
== ====================


The following, more advanced example, can be used to draw lines of constant
properties for n-Pentane. Note the different ways to invoke the
:py:func:`CoolProp.Plots.Plots.PropsPlot.draw_isolines` function draw:

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    ref_fluid = 'n-Pentane'
    plt = PropsPlot(ref_fluid, 'Ts')
    plt.draw_isolines('Q', [0.3, 0.5, 0.7, 0.8])
    plt.draw_isolines('P', [100, 2000], num=5)
    plt.draw_isolines('D', [2, 600], num=7)
    plt.set_axis_limits([-2, 1.5, 200, 500])
    plt.show()

Some of the commonly used `Matplotlib <http://www.matplotlib.org>`_ functions,
such as :func:`title`, :func:`xlabel` and :func:`ylabel` have been wrapped in
the :py:class:`CoolProp.Plots.Plots.PropsPlot` class to make the plotting of
graphs a little simpler, for example:

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    plt = PropsPlot('Water', 'Ts')
    plt.title('Ts Graph for Water')
    plt.xlabel(r's $[{kJ}/{kg K}]$')
    plt.ylabel(r'T $[K]$')
    plt.grid()
    plt.show()

The following two examples show how the :class:`matplotlib.pyplot` functions
and :class:`matplotlib.pyplot.axes` functions can also be used along side
the :py:class:`CoolProp.Plots.Plots.PropsPlot` class

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    plt = PropsPlot('Water', 'Ph')
    ax = plt.axis
    ax.set_yscale('log')
    ax.text(400, 5500, 'Saturated Liquid', fontsize=15, rotation=40)
    ax.text(2700, 3500, 'Saturated Vapour', fontsize=15, rotation=-100)
    plt.show()

.. plot::
    :include-source:

    from matplotlib import pyplot
    from CoolProp.Plots import PropsPlot

    ref_fluid = 'R600a'
    fig = pyplot.figure(1, figsize=(10, 10), dpi=100)
    for i, gtype in enumerate(['PT', 'PD', 'PS', 'PH', 'TD', 'TS', 'HS']):
        ax = pyplot.subplot(4, 2, i)
        if gtype.startswith('P'):
            ax.set_yscale('log')
        plt = PropsPlot(ref_fluid, gtype, axis=ax)
        plt.title(gtype)
        plt._draw_graph()
    pyplot.tight_layout()
    pyplot.show()

