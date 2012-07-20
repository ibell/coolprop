import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as WXCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as WXToolbar
import matplotlib as mpl
import CoolProp as CP
from CoolProp.Plots.Plots import Ph

class PlotPanel(wx.Panel):
    def __init__(self, parent, **kwargs):
        wx.Panel.__init__(self, parent, **kwargs)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.figure = mpl.figure.Figure(dpi=100)
        self.canvas = WXCanvas(self, -1, self.figure)
        self.ax = self.figure.add_axes((0.15,0.15,0.8,0.8))
        #self.toolbar = WXToolbar(self.canvas)
        #self.toolbar.Realize()
        sizer.Add(self.canvas,1,wx.EXPAND)
        #sizer.Add(self.toolbar)
        self.SetSizer(sizer)
        sizer.Layout()
        
class PHPlotFrame(wx.Frame):
    def __init__(self, Fluid):
        wx.Frame.__init__(self, None,title='p-h plot: '+Fluid)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.PP = PlotPanel(self, size = (-1,-1))
        
        sizer.Add(self.PP, 1, wx.EXPAND)
        self.SetSizer(sizer)
        Ph(Fluid, axis = self.PP.ax)
        sizer.Layout()
        
        self.add_menu()
    
    def add_menu(self):
        # Menu Bar
        self.MenuBar = wx.MenuBar()
        self.File = wx.Menu()
        
        mnuItem  = wx.MenuItem(self.File, -1, "Edit...", "", wx.ITEM_NORMAL)
        
        self.File.AppendItem(mnuItem)
        self.MenuBar.Append(self.File, "File")
        
        self.SetMenuBar(self.MenuBar)
        
class MainFrame(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self,None)
        
        self.build()
    
    def build(self):
        # Menu Bar
        self.MenuBar = wx.MenuBar()
        
        self.plots = wx.Menu()
        self.PHPlot  = wx.Menu()#Item(self.plots, -1, "p-h plot", "", wx.ITEM_NORMAL)
        
        for Fluid in sorted(CP.__fluids__):
            mnuItem  = wx.MenuItem(self.PHPlot, -1, Fluid, Fluid, wx.ITEM_NORMAL)
            self.PHPlot.AppendItem(mnuItem)
            self.Bind(wx.EVT_MENU, lambda(event): self.OnPHPlot(event, mnuItem), mnuItem)
        
        self.MenuBar.Append(self.plots, "Plots")
        self.plots.AppendMenu(-1,'p-h plot',self.PHPlot)
        
        self.SetMenuBar(self.MenuBar)
    
    def OnPHPlot(self, event, mnuItem):
        #Make a p-h plot instance in a new frame
        #Get the label (Fluid name)
        Fluid = self.PHPlot.FindItemById(event.Id).Label
        PH = PHPlotFrame(Fluid)
        PH.Show()
        
if __name__=='__main__':
    app = wx.App(False)
    wx.InitAllImageHandlers()
    
    frame = MainFrame() 
    frame.Show(True) 
    app.MainLoop()