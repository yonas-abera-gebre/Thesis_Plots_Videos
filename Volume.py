import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use( 'tkagg' )
from numpy import sin, cos, pi, exp, sqrt, dot
from scipy import integrate
import scipy.special as special
import mpmath as mp
from scipy.signal import find_peaks
from math import floor, ceil
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys

def PEM_Data_Plotter(data):
    
    X, Y, Z = np.mgrid[-10:10:60j, -10:10:60j,-10:10:60j]


    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=PEM_data.flatten(),
        isomin=0.0,
        isomax=1.0,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=17, # needs to be a large number for good volume rendering
        ))
    fig.show()


if __name__=="__main__":
    PEM_Data_Plotter(data)