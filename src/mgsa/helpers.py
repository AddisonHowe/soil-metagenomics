"""Helper functions

"""
import numpy as np
import matplotlib.pyplot as plt

from mgsa.io import pH_soil
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.cm import get_cmap, ScalarMappable

def say_hello():
    print("Hello world")



def plot_old(data, title, cmap = "Blues", vmin = None, vmax = None, fontname="Helvetica"):
    """
    plot_native_perturbed
    inputs: an array with 10 rows (native pHs) and 11 columns (perturbed pHs)
    produces plot
    """
    soils = ["Soil3", "Soil5", "Soil6", "Soil9", "Soil11", "Soil12", "Soil14", "Soil15", "Soil16", "Soil17"]
    native = [4.987, 5.324, 5.405, 5.822, 6.186, 6.255, 6.545, 6.789, 6.86,  7.052]
 
    plt.imshow(data, aspect="auto", cmap=cmap, origin="lower", vmin = vmin, vmax = vmax)
    plt.colorbar()

    plt.xlabel("Perturbed pH (approximate)", fontname="Helvetica", fontsize = "x-large")
    plt.ylabel("Native pH", fontname="Helvetica", fontsize = "x-large")

    x = np.linspace(3.8, 8.4, 11)
    plt.xticks(ticks=np.linspace(0, 10, 11), labels=[f"{val:.1f}" for val in x], rotation=45, fontname='Helvetica', fontsize = 'x-large')
    plt.yticks(ticks=np.linspace(0, 9, 10), labels=[f"{val:.1f}" for val in native], rotation=45, fontname='Helvetica', fontsize = 'x-large')


    plt.title(label=title, fontname='Helvetica', fontsize = 'x-large')


    plt.tight_layout()
    plt.show()
    
def plot(
        data, 
        title, 
        cmap="Blues", 
        norm=None,
        vmin=None, 
        vmax=None, 
        fontname="Helvetica",
        fontsize="x-large", 
        cbarfontsize=None,
        show=True,
        figsize=(8, 6),
        circ_size=200,
        edge_colors="black",
        linewidths=1,
        include_line=False,
        datdir="../data",
        c_bar=True
):
    
    soils = [
        "Soil3", "Soil5", "Soil6", "Soil9", "Soil11", 
        "Soil12", "Soil14", "Soil15", "Soil16", "Soil17"
    ]
    #native = [4.987, 5.324, 5.405, 5.822, 6.186, 6.255, 6.545, 6.789, 6.86,  7.052]
    #some rounding done intentionally so ticklabels can be large and not overlap
    native = [5.0, 5.3, 5.41, 5.8, 6.15, 6.3, 6.5, 6.75, 6.9,  7.1] 

    x = []
    y = []
    for i in range(10):
        pert = pH_soil(soils[i], DATDIR=datdir)
        nat = native[i]*np.ones(11)
        x.append(pert)
        y.append(nat)
        
    x = np.concatenate(x)
    y = np.concatenate(y)
    data = data.flatten()
    
    cmap = get_cmap(cmap)

    vmin = vmin if vmin is not None else np.nanmin(data)
    vmax = vmax if vmax is not None else np.nanmax(data)

    if norm is None:
        norm = Normalize(vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(figsize=figsize)
    scatter = ax.scatter(
        x, y, c=data, cmap=cmap, norm=norm, 
        linewidths=linewidths, 
        edgecolors=edge_colors, 
        s=circ_size,
    )

    if c_bar:
        cbar = plt.colorbar(scatter, ticks=[vmin, vmax])
        cbar.ax.set_yticklabels([f"{vmin:.2f}", f"{vmax:.2f}"], fontsize=cbarfontsize)

    ax.set_xlabel("Perturbed pH", fontname=fontname, fontsize=fontsize)
    ax.set_ylabel("Native pH", fontname=fontname, fontsize=fontsize)

    x = np.linspace(3.8, 8.4, 3)
    y = np.linspace(5, 7, 3)
    ax.set_xticks(
        ticks=x, 
        labels=[f"{val:.0f}" for val in x], 
        rotation=45, fontname=fontname, fontsize=fontsize,
    )
    ax.set_yticks(
        ticks=y, 
        labels=[f"{val:.0f}" for val in y], 
        rotation=45, fontname=fontname, fontsize=fontsize,
    )

    # Add y=x line
    if include_line:
        xlims, ylims = ax.get_xlim(), ax.get_ylim()
        xs = np.linspace(max(xlims[0], ylims[0]), min(xlims[1], ylims[1]), 10)
        ax.plot(xs, xs, "k--", zorder=0)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

    ax.set_title(title, fontname=fontname, fontsize=fontsize)

    if show:
        plt.show()

    return ax


def plot_lfc(
        data_y,
        data_x, 
        title, 
        cmap="Blues", 
        colposinf="blue",
        colneginf="red",
        norm=None,
        vmin=None, 
        vmax=None, 
        fontname="Helvetica",
        fontsize="x-large", 
        cbarfontsize=None,
        show=True,
        figsize=(8, 6),
        circ_size=200,
        edge_colors="black",
        linewidths=1,
        include_line=False,
        datdir="../data",
):
    """Plot log2(y/x) with nan and inf handling"""
    
    soils = [
        "Soil3", "Soil5", "Soil6", "Soil9", "Soil11", 
        "Soil12", "Soil14", "Soil15", "Soil16", "Soil17"
    ]
    #native = [4.987, 5.324, 5.405, 5.822, 6.186, 6.255, 6.545, 6.789, 6.86,  7.052]
    #some rounding done intentionally so ticklabels can be large and not overlap
    native = [5.0, 5.3, 5.41, 5.8, 6.15, 6.3, 6.5, 6.75, 6.9,  7.1] 

    x = []
    y = []
    for i in range(10):
        pert = pH_soil(soils[i], DATDIR=datdir)
        nat = native[i]*np.ones(11)
        x.append(pert)
        y.append(nat)
        
    x = np.concatenate(x)
    y = np.concatenate(y)
    data_y = data_y.flatten()
    data_x = data_x.flatten()

    mask_finite = ((data_y != 0) & (data_x != 0))
    mask_posinf = (data_x == 0) & (data_y != 0)
    mask_neginf = (data_y == 0) & (data_x != 0)

    lfc = np.zeros_like(data_y, dtype=float)
    lfc[mask_finite] = np.log2(data_y[mask_finite] / data_x[mask_finite])
    lfc[mask_posinf] = np.inf
    lfc[mask_neginf] = -np.inf

    # Recompute finite mask to include case where both are 0.
    mask_finite = ((data_y != 0) & (data_x != 0)) | ((data_y == 0) & (data_x == 0))
    
    cmap = get_cmap(cmap)

    vmin = vmin if vmin is not None else np.nanmin(lfc[mask_finite])
    vmax = vmax if vmax is not None else np.nanmax(lfc[mask_finite])
    if vmin >= 0:
        vmin = -1.0
    if vmax <= 0:
        vmax = 1.0

    if norm is None:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    fig, ax = plt.subplots(figsize=figsize)
    scatter = ax.scatter(
        x[mask_finite], y[mask_finite], c=lfc[mask_finite], 
        cmap=cmap, norm=norm, 
        linewidths=linewidths, 
        edgecolors=edge_colors, 
        s=circ_size, 
    )

    ax.scatter(
        x[mask_posinf], y[mask_posinf], c=colposinf,
        linewidths=linewidths, 
        edgecolors=edge_colors, 
        s=circ_size * 0.80, 
    )

    ax.scatter(
        x[mask_neginf], y[mask_neginf], c=colneginf,
        linewidths=linewidths, 
        edgecolors=edge_colors, 
        s=circ_size * 0.80, 
    )

    cbar = plt.colorbar(scatter, ticks=[vmin, vmax])

    cbar.ax.set_yticks(
        [vmin, 0, vmax], 
        [f"{vmin:.2f}", "0", f"{vmax:.2f}"], fontsize=cbarfontsize
    )

    ax.set_xlabel("Perturbed pH", fontname=fontname, fontsize=fontsize)
    ax.set_ylabel("Native pH", fontname=fontname, fontsize=fontsize)

    x = np.linspace(3.8, 8.4, 3)
    y = np.linspace(5, 7, 3)
    ax.set_xticks(
        ticks=x, 
        labels=[f"{val:.0f}" for val in x], 
        rotation=45, fontname=fontname, fontsize=fontsize,
    )
    ax.set_yticks(
        ticks=y, 
        labels=[f"{val:.0f}" for val in y], 
        rotation=45, fontname=fontname, fontsize=fontsize,
    )

    # Add y=x line
    if include_line:
        xlims, ylims = ax.get_xlim(), ax.get_ylim()
        xs = np.linspace(max(xlims[0], ylims[0]), min(xlims[1], ylims[1]), 10)
        ax.plot(xs, xs, "k--", zorder=0)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

    ax.set_title(title, fontname=fontname, fontsize=fontsize)

    if show:
        plt.show()

    return ax

    
def plot10(
        data, title, 
        cmap="Blues", 
        vmin=None, 
        vmax=None, 
        fontname="Helvetica", 
        fontsize="x-large", 
        show=True,
        figsize=(8, 6),
        circ_size=200,
        edge_colors="black",
        linewidths=1,
        include_line=False,
        datdir="../data"
):
    
    soils = [
        "Soil3", "Soil5", "Soil6", "Soil9", "Soil11", 
        "Soil12", "Soil14", "Soil15", "Soil16", "Soil17"
    ]
    
    #native = [4.987, 5.324, 5.405, 5.822, 6.186, 6.255, 6.545, 6.789, 6.86,  7.052]
    # Some rounding done intentionally so ticklabels can be large and not overlap
    native = [5.0, 5.3, 5.41, 5.8, 6.15, 6.3, 6.5, 6.75, 6.9,  7.1] 


        
    x = []
    y = []
    for i in range(10):
        pert = pH_soil(soils[i], DATDIR=datdir)[:10]
        nat = native[i]*np.ones(10)
        x.append(pert)
        y.append(nat)
        
    x = np.concatenate(x)
    y = np.concatenate(y)
    data = data.flatten()
        

    cmap = get_cmap(cmap)

    vmin = vmin if vmin is not None else np.min(data)
    vmax = vmax if vmax is not None else np.max(data)
    norm = Normalize(vmin=vmin, vmax=vmax)


    fig, ax = plt.subplots(figsize=figsize)
    scatter = ax.scatter(
        x, y, c=data, cmap=cmap, norm=norm, 
        linewidths=linewidths, 
        edgecolors=edge_colors, 
        s=circ_size,
    )

    cbar = plt.colorbar(scatter, ticks=[vmin, vmax])
    cbar.ax.set_yticklabels([f"{vmin:.2f}", f"{vmax:.2f}"])

    ax.set_xlabel("Perturbed pH", fontname=fontname, fontsize=fontsize)
    ax.set_ylabel("Native pH", fontname=fontname, fontsize=fontsize)

    # x = np.linspace(3.8, 8, 10)
    x = np.linspace(3.8, 8, 3)
    y = np.linspace(5, 7, 3)
    ax.set_xticks(ticks=x, labels=[f"{val:.0f}" for val in x], rotation=45, fontname=fontname, fontsize=fontsize)
    ax.set_yticks(ticks=y, labels=[f"{val:.0f}" for val in y], rotation=45, fontname=fontname, fontsize=fontsize)
    
    # Add y=x line
    if include_line:
        xlims, ylims = ax.get_xlim(), ax.get_ylim()
        xs = np.linspace(max(xlims[0], ylims[0]), min(xlims[1], ylims[1]), 10)
        ax.plot(xs, xs, "k--", zorder=0)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

    ax.set_title(title, fontname=fontname, fontsize=fontsize)

    if show:
        plt.show()

    return ax
        
        
        
        