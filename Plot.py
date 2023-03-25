
import sys
import json 
import cmath
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# matplotlib.use('TkAgg')
from numpy import sin, log, pi, angle, sqrt
import numpy as np
import mpmath as mp
from scipy.special import sph_harm, lpmv
from scipy import special
from os.path import expanduser
from scipy.signal import find_peaks
import h5py
import plotly.graph_objects as go
from plotly.subplots import make_subplots
sys.path.append('/mpdata/becker/yoge8051/Research/TDSE/General_Functions/')
import Module as Mod
import PModule as PMod
import plotly.graph_objects as go
from math import floor, atan2, acos
from plotly import offline
import k3d
# import plotly.io as pio
# pio.renderers.default = 'vscode'

def closest(lst, k): 
    r = lst[min(range(len(lst)), key = lambda i: abs(float(lst[i])-k))] 
    return r

def Wave_Function_Value(input_par, psi, grid, r, theta, phi, m):
    r = closest(grid, r)
    grid_idx = np.where(grid == r)[0][0]

    l_max_bs = input_par["l_max_bs"]
    return_val = 0.0j

    for l in range(0, l_max_bs):
        psi_idx = grid.size*l + grid_idx
        return_val += 1/r*psi[psi_idx]*sph_harm(m, l, phi, theta)
    
    return return_val

def PAD_Momentum(input_par, psi, grid, m):

    x_axis = np.linspace(-10. , 10., 60)
    y_axis = np.linspace(-10. , 10., 60)
    z_axis = np.linspace(-10. , 10., 60)

    pad_value = np.zeros((z_axis.size,y_axis.size))
    pad_value_magnitude = np.zeros((z_axis.size,y_axis.size, x_axis.size))
    pad_value_phase = np.zeros((z_axis.size,y_axis.size, x_axis.size))
    pad_value_total = np.zeros((z_axis.size,y_axis.size, x_axis.size), dtype=complex)

    z_print = z_axis[0:-1:int(len(z_axis)/10)]
    for i, z in enumerate(z_axis):

        if z in z_print:
            print(round(z,3))

        for j, y in enumerate(y_axis):
            pad_value_temp = 0.0j

            for k, x in enumerate(x_axis):

                r = np.sqrt(x*x + y*y + z*z)
                if r == 0:
                    r = 0.01

                phi = atan2(y, x)
                if phi < 0:
                    phi = 2*pi + phi

                theta = np.arccos(z/r)
                
                val = Wave_Function_Value(input_par, psi, grid, r, theta, phi, m)

                
                # pad_value_temp +=  np.power(np.abs(val),2)
                pad_value_magnitude[i, j, k] = np.power(np.abs(val),2)
                pad_value_phase[i, j, k] = cmath.phase(val)

                pad_value_total[i, j, k] = val
            # pad_value[j, i] = pad_value_temp.real

    # pad_value_save = pad_value_save / pad_value_save.max()

    # pad_value_magnitude_reshaped = pad_value_magnitude.reshape(pad_value_magnitude.shape[0], -1)
    # pad_value_phase_reshaped = pad_value_phase.reshape(pad_value_phase.shape[0], -1)

    # np.savetxt("O2_M1_13_Mag.txt",  pad_value_magnitude_reshaped)
    # np.savetxt("O2_M1_13_Phase.txt",  pad_value_phase_reshaped)

    pad_value_magnitude /= pad_value_magnitude.max()

    np.save("H2_M1_1_Mag", pad_value_magnitude)
    np.save("H2_M1_1_Phase", pad_value_phase)
    np.save("H2_M1_1_Total", pad_value_total)


def Psi_Plotter(input_par, psi, grid):
    block_to_qn, qn_to_block = Mod.Index_Map(input_par)
    n_max = input_par["n_max"]
    n_values = np.arange(1, n_max + 1)

    psi_grid = np.zeros(len(grid))
    for idx in qn_to_block:
        psi_grid += np.abs(psi[idx])

    return psi_grid

def Target_File_Reader_WO_Parity(input_par):
    file = h5py.File(file_location + input_par["Target_File"], 'r')
    energy = {}
    wave_function = {}
    for m in range(input_par["m_max_bs"] + 1):
        n_quantum_number = 1
        for i in range(input_par["n_max"]):
            
            energy_temp = file["Energy_" + str(abs(m)) + "_" + str(i)]
            energy_temp = np.array(energy_temp[:,0] + 1.0j*energy_temp[:,1])
            wave_function_temp = file["Psi_" + str(abs(m)) + "_" + str(i)]
            wave_function_temp = np.array(wave_function_temp[:,0] + 1.0j*wave_function_temp[:,1])
            
            energy[(n_quantum_number, m)] = energy_temp
            wave_function[(n_quantum_number, m)] = wave_function_temp
            n_quantum_number += 1
    
    return energy, wave_function 

def Plot_Wave_Function():
    


    X, Y, Z = np.mgrid[-2.5:2.5:80j, -2.5:2.5:80j, -3.5:3.5:80j]

    data = np.load("H2.npy")

    data /= data.max()

 
    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=data.flatten(),
        isomin=0.1,
        isomax=0.8,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=17, # needs to be a large number for good volume rendering
        ))
    offline.plot(fig)


def Plot_Wave_Function2():
    


    X, Y, Z = np.mgrid[-2.5:2.5:80j, -2.5:2.5:80j, -3.5:3.5:80j]

    data = np.load("H2.npy")

    data /= data.max()

 
    plt_volume = k3d.volume(data)

    plot = k3d.plot()
    plot += plt_volume
    plot.display()




if __name__=="__main__":
    file_location = "/mpdata/becker/yoge8051/Research/Data/Diatomic/Hydrogen/"
    input_par = Mod.Input_File_Reader(file_location + "input.json")
    grid = Mod.Grid(input_par["grid_spacing"], input_par["grid_size"])
    grid= grid.grid
    energy, bound_states = Target_File_Reader_WO_Parity(input_par)

    m = 1
    n = 1
    psi = bound_states[(n, m)]

    print(energy[(n, m)])
    PAD_Momentum(input_par, psi, grid, m)

    
    # Plot_Wave_Function2()
