import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as WXCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as WXToolbar
import matplotlib as mpl
import CoolProp as CP
from CoolProp.Plots.Plots import Ph, Ts

# Munge the system path if necessary to add the lib folder (only really needed
# for packaging using cx_Freeze)
#if os.path.exists('lib') and os.path.abspath(os.path.join(os.curdir,'lib')) not in os.:

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

class TSPlotFrame(wx.Frame):
    def __init__(self, Fluid):
        wx.Frame.__init__(self, None,title='T-s plot: '+Fluid)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.PP = PlotPanel(self, size = (-1,-1))
        
        sizer.Add(self.PP, 1, wx.EXPAND)
        self.SetSizer(sizer)
        Ts(str(Fluid), axis = self.PP.ax)
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
                
class PHPlotFrame(wx.Frame):
    def __init__(self, Fluid):
        wx.Frame.__init__(self, None,title='p-h plot: '+Fluid)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.PP = PlotPanel(self, size = (-1,-1))
        
        sizer.Add(self.PP, 1, wx.EXPAND)
        self.SetSizer(sizer)
        Ph(str(Fluid), axis = self.PP.ax)
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
        self.PHPlot = wx.Menu()
        self.TSPlot = wx.Menu()
        
        for Fluid in sorted(CP.__fluids__):
            mnuItem  = wx.MenuItem(self.PHPlot, -1, Fluid, "", wx.ITEM_NORMAL)
            self.PHPlot.AppendItem(mnuItem)
            self.Bind(wx.EVT_MENU, lambda(event): self.OnPHPlot(event, mnuItem), mnuItem)
            
            mnuItem  = wx.MenuItem(self.TSPlot, -1, Fluid, "", wx.ITEM_NORMAL)
            self.TSPlot.AppendItem(mnuItem)
            self.Bind(wx.EVT_MENU, lambda(event): self.OnTSPlot(event, mnuItem), mnuItem)
        
        self.MenuBar.Append(self.plots, "Plots")
        self.plots.AppendMenu(-1,'p-h plot',self.PHPlot)
        self.plots.AppendMenu(-1,'T-s plot',self.TSPlot)
        
        self.SetMenuBar(self.MenuBar)
    
    def OnPHPlot(self, event, mnuItem):
        #Make a p-h plot instance in a new frame
        #Get the label (Fluid name)
        Fluid = self.PHPlot.FindItemById(event.Id).Label
        PH = PHPlotFrame(Fluid)
        PH.Show()
        
    def OnTSPlot(self, event, mnuItem):
        #Make a p-h plot instance in a new frame
        #Get the label (Fluid name)
        Fluid = self.TSPlot.FindItemById(event.Id).Label
        TS = TSPlotFrame(Fluid)
        TS.Show()
        
if __name__=='__main__':
    app = wx.App(False)
    wx.InitAllImageHandlers()
    
    frame = MainFrame() 
    frame.Show(True) 
    app.MainLoop()