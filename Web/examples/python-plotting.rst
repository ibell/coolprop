.. _python-plotting:

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
