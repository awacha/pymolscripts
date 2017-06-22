import pymol
import sys
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends import backend_tkagg

def _new_figure_manager(num, *args, **kwargs):
    if pymol._ext_gui is None:
        return new_figure_manager(num, *args, **kwargs)
    backend_tkagg.show._needmain = False
    if sys.version().major >= 3:
        import tkinter as Tk
    else:
        import Tkinter as Tk
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, FigureManagerTkAgg
    FigureClass = kwargs.pop('FigureClass', Figure)
    figure = FigureClass(*args, **kwargs)
    window = Tk.Toplevel(master=pymol._ext_gui.root)
    canvas = FigureCanvasTkAgg(figure, master=window)
    figManager = FigureManagerTkAgg(canvas, num, window)
    if matplotlib.is_interactive():
        figManager.show()
    return figManager

new_figure_manager = backend_tkagg.new_figure_manager
backend_tkagg.new_figure_manager=_new_figure_manager