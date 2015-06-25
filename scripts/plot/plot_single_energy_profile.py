'''Script to plot energy profile'''
import threading

from simple_plot import *
from plot_input import *


class SinglePlotThread(threading.Thread):
    "Sub thread to plot single energy profile."
    def __init__(self, func, args):
        threading.Thread.__init__(self)
        self.func = func
        self.args = args

    def run(self):
        apply(self.func, self.args)


