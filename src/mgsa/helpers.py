"""Helper functions

"""
import numpy as np
import matplotlib.pyplot as plt

from mgsa.io import pH_soil
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap, ScalarMappable

def say_hello():
    print("Hello world")



def plot_old(data, title, cmap = 'Blues', vmin = None, vmax = None, fontname='Helvetica'):
    '''
    plot_native_perturbed
    inputs: an array with 10 rows (native pHs) and 11 columns (perturbed pHs)
    produces plot
    '''
    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    native = [4.987, 5.324, 5.405, 5.822, 6.186, 6.255, 6.545, 6.789, 6.86,  7.052]
 
    plt.imshow(data, aspect='auto', cmap=cmap, origin='lower', vmin = vmin, vmax = vmax)
    plt.colorbar()

    plt.xlabel('Perturbed pH (approximate)', fontname='Helvetica', fontsize = 'x-large')
    plt.ylabel('Native pH', fontname='Helvetica', fontsize = 'x-large')

    x = np.linspace(3.8, 8.4, 11)
    plt.xticks(ticks=np.linspace(0, 10, 11), labels=[f"{val:.1f}" for val in x], rotation=45, fontname='Helvetica', fontsize = 'x-large')
    plt.yticks(ticks=np.linspace(0, 9, 10), labels=[f"{val:.1f}" for val in native], rotation=45, fontname='Helvetica', fontsize = 'x-large')


    plt.title(label=title, fontname='Helvetica', fontsize = 'x-large')


    plt.tight_layout()
    plt.show()
    
def plot(data, title, cmap = 'Blues', vmin = None, vmax = None, fontname = 'Helvetica'):
    
    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    #native = [4.987, 5.324, 5.405, 5.822, 6.186, 6.255, 6.545, 6.789, 6.86,  7.052]
    native = [5.0, 5.3, 5.41, 5.8, 6.15, 6.3, 6.5, 6.75, 6.9,  7.1] #some rounding done intentionally so ticklabels can be large and not overlap


        
    x = []
    y = []
    for i in range(10):
        pert = pH_soil(soils[i])
        nat = native[i]*np.ones(11)
        x.append(pert)
        y.append(nat)
        
    x = np.concatenate(x)
    y = np.concatenate(y)
    data = data.flatten()
        

    cmap = get_cmap(cmap)

    vmin = vmin if vmin is not None else np.min(data)
    vmax = vmax if vmax is not None else np.max(data)
    norm = Normalize(vmin=vmin, vmax=vmax)


    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(x, y, c=data, cmap=cmap, norm=norm, 
                        edgecolors='black', linewidths=1, s = 200)


    cbar = plt.colorbar(scatter, ticks=[vmin, vmax])
    cbar.ax.set_yticklabels([f'{vmin:.1f}', f'{vmax:.1f}'])

    ax.set_xlabel('Perturbed pH', fontname='Helvetica', fontsize = 'x-large')
    ax.set_ylabel('Native pH', fontname='Helvetica', fontsize = 'x-large')

    x = np.linspace(3.8, 8.4, 11)
    ax.set_xticks(ticks=x, labels=[f"{val:.1f}" for val in x], rotation=45, fontname=fontname, fontsize = 'x-large')
    ax.set_yticks(ticks=native, labels=[f"{val:.1f}" for val in native], rotation=45, fontname=fontname, fontsize = 'x-large')


    ax.set_title(title, fontname='Helvetica', fontsize = 'x-large')


    plt.show()
        
        
        
        