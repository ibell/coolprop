Python Plotting
===============

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    plt = PropsPlot('R290', 'Ts')
    plt.show()

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    plt = PropsPlot('R410A', 'Ph')
    plt.show()

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

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    plt = PropsPlot('Water', 'Ts')
    plt.title('Ts Graph for Water')
    plt.xlabel(r's $[{kJ}/{kg K}]$')
    plt.ylabel(r'T $[K]$')
    plt.grid()
    plt.show()

.. plot::
    :include-source:

    from CoolProp.Plots import PropsPlot

    plt = PropsPlot('Water', 'Ph')
    ax = plt.axis
    ax.set_yscale('log')
    ax.text(400, 5500, 'Saturated Liquid', fontsize=15, rotation=40)
    ax.text(2700, 3500, 'Saturated Vapour', fontsize=15, rotation=-100)
    plt.show()
