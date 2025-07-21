"""Helper functions

"""
import numpy as np
import matplotlib.pyplot as plt

def say_hello():
    print("Hello world")


'''
plot_native_perturbed
inputs: an array with 10 rows (native pHs) and 11 columns (perturbed pHs)
produces plot
'''

def plot(data, title, cmap = 'Blues', vmin = None, vmax = None, fontname='Helvetica'):
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
    