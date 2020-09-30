import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import ndimage as ndi
from skimage import feature


emin = 'default'
emax = 'default'
sigma = 'default'
output_path = 'default'
parameters_string = "#Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font= 'helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nimage\npolygon"




def sigma(seqid,mod,emin,emax,sigma):
    '''Opens event file, performs energy cut, performs and plots canny edge detection for chosen sigma and returns plots'''
    evt_file = 'nu'+str(seqid)+str(mod)+'01_cl.evt'
    hdul = fits.open(evt_file)
    data = hdul[1].data
    DET1X = data["DET1X"]
    DET1Y = data["DET1Y"]
    PI = data["PI"]
    if emin == 'default':
        emin = 36
    else:
        emin = emin
    if emax == 'default':
        emax == 209
    else:
        emax = emax
    del_PI_indices = [i for i, x in enumerate(PI) if x<=emin or x>=emax]
    cut_DET1X = np.delete(DET1X, del_PI_indices)
    cut_DET1Y = np.delete(DET1Y, del_PI_indices)
    cut_counts_arrays = np.histogram2d(cut_DET1X, cut_DET1Y, [360,360], range=[[0,360],[0,360]])
    cut_counts = np.hstack(cut_counts_arrays[0])
    cut_counts_binned = np.split(cut_counts,360)
    im = np.column_stack(cut_counts_binned)
    if sigma == 'default':
        sigma = 10*np.std(cut_counts)
    else:
        sigma = sigma
    edges = feature.canny(im, sigma=sigma)
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 5), sharex=True, sharey=True)

    ax1.imshow(im, cmap=plt.cm.gray)
    ax1.axis('off')
    ax1.set_title('Original', fontsize=15)

    ax2.imshow(edges, cmap=plt.cm.gray)
    ax2.axis('off')
    ax2.set_title(r'Canny Filter, $\sigma='+str(sigma)+'$', fontsize=15)
    fig.tight_layout()
    plot = plt.show()
    return plot



def sigmarange(seqid,mod,emin,emax):
    '''Opens event file, performs energy cut, performs and plots canny edge detection for range of sigmas and returns plots'''
    evt_file = 'nu'+seqid+mod+'01_cl.evt'
    hdul = fits.open(evt_file)
    data = hdul[1].data
    DET1X = data["DET1X"]
    DET1Y = data["DET1Y"]
    PI = data["PI"]
    if emin == 'default':
        emin = 36
    else:
        emin = emin
    if emax == 'default':
        emax == 209
    else:
        emax = emax
    del_PI_indices = [i for i, x in enumerate(PI) if x<=emin or x>=emax]
    cut_DET1X = np.delete(DET1X, del_PI_indices)
    cut_DET1Y = np.delete(DET1Y, del_PI_indices)
    cut_counts_arrays = np.histogram2d(cut_DET1X, cut_DET1Y, [360,360], range=[[0,360],[0,360]])
    cut_counts = np.hstack(cut_counts_arrays[0])
    cut_counts_binned = np.split(cut_counts,360)
    im = np.column_stack(cut_counts_binned)
    sigma_values = np.concatenate(([0], np.arange(3,12.5,0.5)))
    edges_array = [im]
    for s in sigma_values:
        edges = feature.canny(im,sigma=s)
        edges_array.append(edges)
    nr, nc = 5, 4
    assert len(sigma_values) == nr * nc, "number of sigma values should match grid cells"
    fig, axs = plt.subplots(nrows=nr, ncols=nc, figsize=(12, 15), sharex=True, sharey=True)
    for i, (s, edge) in enumerate(zip(sigma_values, edges_array)):
        if i == 0:
            title = 'Original'
        else:
            title = r'Canny Filter, $\sigma={}$'.format(s)
        ax = axs[i // nc][i % nc]
        ax.imshow(edge, cmap=plt.cm.gray)
        ax.axis('off')
        ax.set_title(title, fontsize=12)
    #fig.tight_layout()
    plot = plt.show()
    return plot

def write_reg_file(seqid,mod,emin,emax,sigma,output_path,name_of_region_file):
    '''Opens event file, performs energy cut, performs and plots canny edge detection for given sigma,identifies corners,writes/overwrites region file'''
    evt_file = 'nu'+str(seqid)+str(mod)+'01_cl.evt'
    hdul = fits.open(evt_file)
    data = hdul[1].data
    DET1X = data["DET1X"]
    DET1Y = data["DET1Y"]
    PI = data["PI"]
    if emin == 'default':
        emin = 36
    else:
        emin = emin
    if emax == 'default':
        emax == 209
    else:
        emax = emax
    del_PI_indices = [i for i, x in enumerate(PI) if x<=emin or x>=emax]
    cut_DET1X = np.delete(DET1X, del_PI_indices)
    cut_DET1Y = np.delete(DET1Y, del_PI_indices)
    cut_counts_arrays = np.histogram2d(cut_DET1X, cut_DET1Y, [360,360], range=[[0,360],[0,360]])
    cut_counts = np.hstack(cut_counts_arrays[0])
    cut_counts_binned = np.split(cut_counts,360)
    im = np.column_stack(cut_counts_binned)
    if sigma == 'default':
        sigma = 10*np.std(cut_counts)
    else:
        sigma = sigma
    edges = feature.canny(im, sigma=sigma)
    indices = np.where(edges != [0])
    x_coords = list(indices[1])
    y_coords = list(indices[0])
    region_cols = list(set(x_coords))
    curve_perimeter_coords_1 = []
    for i in region_cols:
        if i in curve_perimeter_coords_1:
            pass
        else:
            curve_perimeter_coords_1.append(i)
            y_coord_indices =[j for j, x in enumerate(x_coords) if x == i]
            y_coords_for_i = []
            for i in y_coord_indices:
                y_coords_for_i.append(y_coords[i])
            curve_perimeter_coords_1.append(max(y_coords_for_i))
    corners = []
    for i in curve_perimeter_coords_1:
        corners.append(i)
    region_cols_2 = region_cols[::-1]
    curve_perimeter_coords_2 = []
    for i in region_cols_2:
        if i in curve_perimeter_coords_2:
            pass
        else:
            curve_perimeter_coords_2.append(i)
            y_coord_indices =[j for j, x in enumerate(x_coords) if x == i]
            y_coords_for_i = []
            for i in y_coord_indices:
                y_coords_for_i.append(y_coords[i])
            curve_perimeter_coords_2.append(min(y_coords_for_i))
    for i in curve_perimeter_coords_2:
        corners.append(i)
    if output_path == 'default':
        output_path = os.getcwd()
    else:
        output_path = os.chdir(output_path)
    region_file = open(name_of_region_file, "w")
    region_file.write(parameters_string + str(tuple(corners)))
    region_file.close()
