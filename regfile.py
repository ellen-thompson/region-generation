import os
import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
from astropy.io import fits
import astropy.units as u
from nustar_gen import info, utils
ns = info.NuSTAR()
from scipy import ndimage as ndi
from skimage import feature


parameters_string = "#Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font= 'helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nimage\npolygon"


def image(seqid,mod,scale):
    '''Opens event file, returns plotted image of event file w/ lin or log scale'''
    evt_file = 'nu'+str(seqid)+str(mod)+'01_cl.evt'
    hdul = fits.open(evt_file)
    data = hdul[1].data
    DET1X = data["DET1X"]
    DET1Y = data["DET1Y"]
    counts_arrays = np.histogram2d(DET1X, DET1Y, [360,360], range=[[0,360],[0,360]])
    counts = np.hstack(counts_arrays[0])
    counts_binned = np.split(counts,360)
    im = np.column_stack(counts_binned)
    fig,ax = plt.subplots(1,figsize=(7,7))
    ax.axis('off')
    assert scale == 'log' or scale == 'lin',"Scale must be 'log' or 'lin'"
    if scale == 'log':
        my_cmap = copy.copy(plt.cm.get_cmap('viridis'))
        my_cmap.set_bad((0,0,0))
        ax.imshow(im,norm=colors.LogNorm(),interpolation='nearest',cmap=my_cmap, origin = 'lower')
    if scale == 'lin':
        ax.imshow(im,interpolation='nearest',cmap='viridis', origin = 'lower')

    plot = plt.show()
    return plot

def filter_source(seqid,mod,limit,scale):
    '''Opens source and event files, interpolates source positions onto event times, filters event file and writes filtered file, plots image with point source removed'''
    sky2det = 'nu'+str(seqid)+str(mod)+'_sky2det.fits'
    hdul = fits.open(sky2det)
    src_data = hdul[1].data
    src_time = src_data["TIME"]
    src_DET1X = src_data["DET1X"]
    src_DET1Y = src_data["DET1Y"]
    evt_file = 'nu'+str(seqid)+str(mod)+'01_cl.evt'
    hdul = fits.open(evt_file)
    evt_data = hdul[1].data
    evt_time = evt_data["TIME"]
    evt_DET1X = evt_data["DET1X"]
    evt_DET1Y = evt_data["DET1Y"]
    xr = np.interp(evt_time,src_time,src_DET1X)
    yr = np.interp(evt_time,src_time,src_DET1Y)
    dr = np.sqrt((evt_DET1X-xr)**2 + (evt_DET1Y - yr)**2)
    limit = limit*u.arcmin
    limit_pix = (limit/ns.pixel).cgs
    filt = np.where((dr>limit_pix))[0]
    keep_DET1X = evt_DET1X[filt]
    keep_DET1Y = evt_DET1Y[filt]
    hdul = fits.open(evt_file)
    hdul[1].data = evt_data[filt]
    hdul.writeto('nu'+str(seqid)+'filt'+str(mod)+'01_cl.evt',overwrite=True)
    with fits.open('nu'+str(seqid)+str(mod)+'_cl.evt', mode='update') as hdul:
        hdul.flush()
    keep_counts_arrays = np.histogram2d(keep_DET1X, keep_DET1Y, [360,360], range=[[0,360],[0,360]])
    keep_counts = np.hstack(keep_counts_arrays[0])
    keep_counts_binned = np.split(keep_counts,360)
    im = np.column_stack(keep_counts_binned)
    fig,ax = plt.subplots(1,figsize=(7,7))
    ax.axis('off')
    assert scale == 'log' or scale == 'lin',"Scale must be 'log' or 'lin'"
    if scale == 'log':
        my_cmap = copy.copy(plt.cm.get_cmap('viridis'))
        my_cmap.set_bad((0,0,0))
        ax.imshow(im,norm=colors.LogNorm(),interpolation='nearest',cmap=my_cmap, origin = 'lower')
    if scale == 'lin':
        ax.imshow(im,interpolation='nearest',cmap='viridis', origin = 'lower')

    plot = plt.show()
    return plot


def sigma(seqid,mod,emin,emax,sigma):
    '''Opens event file, performs energy cut, performs and plots canny edge detection for chosen sigma and returns plots'''
    evt_file = 'nu'+str(seqid)+str(mod)+'01_cl.evt'
    hdul = fits.open(evt_file)
    data = hdul[1].data
    DET1X = data["DET1X"]
    DET1Y = data["DET1Y"]
    PI = data["PI"]
    del_PI_indices = [i for i, x in enumerate(PI) if x<=emin or x>=emax]
    cut_DET1X = np.delete(DET1X, del_PI_indices)
    cut_DET1Y = np.delete(DET1Y, del_PI_indices)
    cut_counts_arrays = np.histogram2d(cut_DET1X, cut_DET1Y, [360,360], range=[[0,360],[0,360]])
    cut_counts = np.hstack(cut_counts_arrays[0])
    cut_counts_binned = np.split(cut_counts,360)
    im = np.column_stack(cut_counts_binned)
    edges = feature.canny(im, sigma=sigma)
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 5), sharex=True, sharey=True)

    ax1.imshow(im, cmap=plt.cm.gray, origin = 'lower')
    ax1.axis('off')
    ax1.set_title('Original', fontsize=15)

    ax2.imshow(edges, cmap=plt.cm.gray, origin = 'lower')
    ax2.axis('off')
    ax2.set_title(r'Canny Filter, $\sigma='+str(sigma)+'$', fontsize=15)
    fig.tight_layout()
    plot = plt.show()
    return plot



def sigma_range(seqid,mod,emin,emax):
    '''Opens event file, performs energy cut, performs and plots canny edge detection for range of sigmas and returns plots'''
    evt_file = 'nu'+seqid+mod+'01_cl.evt'
    hdul = fits.open(evt_file)
    data = hdul[1].data
    DET1X = data["DET1X"]
    DET1Y = data["DET1Y"]
    PI = data["PI"]
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
        ax.imshow(edge, cmap=plt.cm.gray, origin = 'lower')
        ax.axis('off')
        ax.set_title(title, fontsize=12)
    #fig.tight_layout()
    plot = plt.show()
    return plot

def check_region(seqid,mod,emin,emax,sigma,scale):
    '''Opens event file, performs energy cut, identifies polygon region corners from canny edge detection,plots region over original image in linear scale'''
    evt_file = 'nu'+str(seqid)+str(mod)+'01_cl.evt'
    hdul = fits.open(evt_file)
    data = hdul[1].data
    DET1X = data["DET1X"]
    DET1Y = data["DET1Y"]
    PI = data["PI"]
    del_PI_indices = [i for i, x in enumerate(PI) if x<=emin or x>=emax]
    cut_DET1X = np.delete(DET1X, del_PI_indices)
    cut_DET1Y = np.delete(DET1Y, del_PI_indices)
    cut_counts_arrays = np.histogram2d(cut_DET1X, cut_DET1Y, [360,360], range=[[0,360],[0,360]])
    cut_counts = np.hstack(cut_counts_arrays[0])
    cut_counts_binned = np.split(cut_counts,360)
    im = np.column_stack(cut_counts_binned)
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
    polygon_corners = []

    for i in range(0,len(corners)-1,2):
        corner = [corners[i],corners[i+1]]
        polygon_corners.append(corner)
    fig,ax = plt.subplots(1,figsize=(7,7))
    ax.axis('off')
    assert scale == 'log' or scale == 'lin',"Scale must be 'log' or 'lin'"
    if scale == 'log':
        my_cmap = copy.copy(plt.cm.get_cmap('viridis'))
        my_cmap.set_bad((0,0,0))
        ax.imshow(im,norm=colors.LogNorm(),interpolation='nearest',cmap=my_cmap, origin = 'lower')
    if scale == 'lin':
        ax.imshow(im,interpolation='nearest',cmap='viridis', origin = 'lower')
    polygon_region = patches.Polygon(polygon_corners, edgecolor = 'white', facecolor ='none',linewidth = 2)
    ax.add_patch(polygon_region)
    plot = plt.show()
    return plot



def write_regfile(seqid,mod,emin,emax,sigma,output_path,name_of_region_file):
    '''Opens event file, performs energy cut,identifies polygon region corners from canny edge detection, writes/overwrites region file'''
    evt_file = 'nu'+str(seqid)+str(mod)+'01_cl.evt'
    hdul = fits.open(evt_file)
    data = hdul[1].data
    DET1X = data["DET1X"]
    DET1Y = data["DET1Y"]
    PI = data["PI"]
    del_PI_indices = [i for i, x in enumerate(PI) if x<=emin or x>=emax]
    cut_DET1X = np.delete(DET1X, del_PI_indices)
    cut_DET1Y = np.delete(DET1Y, del_PI_indices)
    cut_counts_arrays = np.histogram2d(cut_DET1X, cut_DET1Y, [360,360], range=[[0,360],[0,360]])
    cut_counts = np.hstack(cut_counts_arrays[0])
    cut_counts_binned = np.split(cut_counts,360)
    im = np.column_stack(cut_counts_binned)
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
    os.chdir(output_path)
    region_file = open(name_of_region_file, "w")
    region_file.write(parameters_string + str(tuple(corners)))
    region_file.close()

def save_image(seqid,mod,emin,emax,sigma,scale,output_path,name_of_image):
    '''Opens event file, performs energy cut, identifies polygon region corners from canny edge detection,plots region over original image in linear scale'''
    evt_file = 'nu'+str(seqid)+str(mod)+'01_cl.evt'
    hdul = fits.open(evt_file)
    data = hdul[1].data
    DET1X = data["DET1X"]
    DET1Y = data["DET1Y"]
    PI = data["PI"]
    del_PI_indices = [i for i, x in enumerate(PI) if x<=emin or x>=emax]
    cut_DET1X = np.delete(DET1X, del_PI_indices)
    cut_DET1Y = np.delete(DET1Y, del_PI_indices)
    cut_counts_arrays = np.histogram2d(cut_DET1X, cut_DET1Y, [360,360], range=[[0,360],[0,360]])
    cut_counts = np.hstack(cut_counts_arrays[0])
    cut_counts_binned = np.split(cut_counts,360)
    im = np.column_stack(cut_counts_binned)
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
    polygon_corners = []

    for i in range(0,len(corners)-1,2):
        corner = [corners[i],corners[i+1]]
        polygon_corners.append(corner)
    fig,ax = plt.subplots(1,figsize=(7,7))
    ax.axis('off')
    assert scale == 'log' or scale == 'lin',"Scale must be 'log' or 'lin'"
    if scale == 'log':
        my_cmap = copy.copy(plt.cm.get_cmap('viridis'))
        my_cmap.set_bad((0,0,0))
        ax.imshow(im,norm=colors.LogNorm(),interpolation='nearest',cmap=my_cmap, origin = 'lower')
    if scale == 'lin':
        ax.imshow(im,interpolation='nearest',cmap='viridis', origin = 'lower')
    polygon_region = patches.Polygon(polygon_corners, edgecolor = 'white', facecolor ='none',linewidth = 2)
    ax.add_patch(polygon_region)
    os.chdir(output_path)
    plt.savefig(name_of_image)
    return
