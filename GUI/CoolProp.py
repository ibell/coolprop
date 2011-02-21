
import matplotlib   
matplotlib.use('GTKAgg')

from matplotlib.figure import Figure   
from matplotlib.axes import Subplot   
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg, NavigationToolbar2GTKAgg 
from matplotlib.mlab import csv2rec
import numpy as np
from scipy.optimize import fmin_bfgs,fmin_cg
import sys,copy,os

try:
    import pygtk
except:
    print 'No PyGTK'
    sys.exit(1)
    
try:
    import gtk
    import gobject
except:
    print 'No gtk'
    sys.exit(1)

def set_model_from_list (cb, items):
    """Setup a ComboBox or ComboBoxEntry based on a list of strings."""           
    model = gtk.ListStore(str)
    for i in items:
        model.append([i])
    cb.set_model(model)
    if type(cb) == gtk.ComboBoxEntry:
        cb.set_text_column(0)
    elif type(cb) == gtk.ComboBox:
        cell = gtk.CellRendererText()
        cb.pack_start(cell, True)
        cb.add_attribute(cell, 'text', 0)

class CoolPropGUIGTK:
    
    """This is an Hello World GTK application"""
    def __init__(self):
        self.wTree=gtk.Builder()
        self.wTree.add_from_file('CoolPropGUI.glade')
        
        cb = self.wTree.get_object("cmbRef")
        set_model_from_list(cb,['R410A','R404A','R744 (CO2)','R717 (Ammonia)','R134a','R290','R32'])
        
        cb = self.wTree.get_object("cmbOutput")
        set_model_from_list(cb,['Pressure','Temperature','Enthalpy','Density','Entropy'])
        
        cb = self.wTree.get_object("cmbParam1")
        set_model_from_list(cb,['Temperature','Pressure'])
   
        cb = self.wTree.get_object("cmbParam2")
        set_model_from_list(cb,['Pressure','Temperature','Enthalpy','Density','Entropy'])

        dic = {"on_window1_destroy" : gtk.main_quit}
        self.wTree.connect_signals(dic)
        
if __name__ == "__main__":
    hwg = CoolPropGUIGTK()
    gtk.main()