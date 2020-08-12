import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.io import fits
from astropy.table import Table
from itertools import compress

#get in directory w/ data (this observation has "Easy" stray light)
cd Downloads/40111002002/event_cl/

#open data w/ FITS
hdul = fits.open("nu40111002002A01_cl.evt")

#move into 2nd hdu
data = hdul[1].data

#make table of columns
tab = Table(data).columns

#define variables for det1x, det1y, and detid columns
DET1X = data["DET1X"]
DET1Y = data["DET1Y"]
DET_ID = data["DET_ID"]

#make 2D histogram of det1x and det1y for visualization
NBINS=(360,360)
fig,ax = plt.subplots(1)
ax.hist2d(data['DET1X'], data['DET1Y'], NBINS, range=[[1, 360], [1,360]], cmap='viridis')
ax.set_xlabel('DET1X')
ax.set_ylabel('DET1Y')
#plt.show()
plt.savefig('2dhist.png')
#save file for reference with threshhold histogram & highlighted region image

#define function that returns #counts for a given bin
def counts(xlower, xupper, ylower, yupper):
    x = tab['DET1X']
    y = tab['DET1Y']
    events = np.where((x > xlower) & (x < xupper) & (y > ylower) & (y < yupper))
    return (np.size(events))

#loop through 10x10 bins and create array of counts
xlowervalues = list(range(0, 360, 10))
xuppervalues = list(range(10, 370, 10))
ylowervalues = list(range(0, 360, 10))
yuppervalues = list(range(10, 370, 10))
def counts_array():
    array =[]
    for i in xlowervalues:
        for j in xuppervalues:
                if j - i == 10:
                    for k in ylowervalues:
                        for l in yuppervalues:
                            if l - k == 10:
                                a = counts(i, j, k, l)
                                array.append(a)
    return(array)

#make counts_array function a variable
array = counts_array()

#

#plot histogram
plt.hist(array, range=[1, 360], bins = 360)
plt.xlabel("Per Bin Count Rate")
plt.ylabel("Frequency of Counts")
#plt.show()
plt.savefig('threshold.png')

#identify neighbors(up, down, left, and right bins of a given bin)
def neighbors(i):
    neighbors_array = []
    for neighbor in [i+1, i-1, i-36, i+36]:
        neighbors_array.append(array[neighbor % len(array)])
    return neighbors_array

#call sl threshold ~250 from histogram, if bin has at least 1 neighbor at or above SL threshold, count as in SL region
#boolean region test
def region_test(i):
    if neighbors(i)[0] >= 250 or neighbors(i)[1] >= 250 or neighbors(i)[2] >= 250 or neighbors(i)[3] >= 250:
        i = True
    else: i = False
    return i

#make array of passing bin indices
region_bins = [i for i, x in enumerate(array) if region_test(i)]

#visualizing the region
#plot patches over SL bins if region test True
indices = range(1, 1295, 1)
NBINS=(360,360)
fig,ax = plt.subplots(1)
ax.hist2d(data['DET1X'], data['DET1Y'], NBINS, range=[[1, 360], [1,360]], cmap='gray')
ax.set_xlabel('DET1X')
ax.set_ylabel('DET1Y')
for i in region_bins:
    cornerx = (indices[i]// 36)*10
    cornery = (indices[i] % 36)*10
    SL = patches.Rectangle((cornerx,cornery), 10, 10, edgecolor='none', facecolor='yellow', alpha = 0.4)
    ax.add_patch(SL)
#plt.savefig('region.png')
plt.show()

#define initial Figure Of Merit rectangle sides
mincol = min(region_bins)//36
maxcol = max(region_bins)//36
minrow = min(region_bins)%36
maxrow = max(region_bins)%36

#test FOM, check with bin count to optimize rectangle.
#for this file, the min and max cols, minrow-1, and max row give best FOM value
all_bins = range(1,1297,1)
def rectangle():
    rectangle = []
    for i in all_bins:
        if i//36 >= (mincol) and i//36 <= (maxcol) and i%36 >= (minrow-1) and i%36 <= (maxrow):
            rectangle.append(i)
    return rectangle
rectangle = rectangle()

def bincount():
    bincount = 0
    for i in rectangle:
        if i in region_bins:
            bincount += 1
    return bincount

def FOM():
    FOM = 0
    for i in rectangle:
        if i in region_bins:
            FOM += 1
        if i not in region_bins:
            FOM += -1
    return FOM

#visualization: plot binned image w/rectangle (did some tweaking for later work, can disregard here)
indices = range(1, 1295, 1)
NBINS=(360,360)
fig,ax = plt.subplots(1)
ax.hist2d(data['DET1X'], data['DET1Y'], NBINS, range=[[1, 360], [1,360]], cmap='gray')
ax.set_xlabel('DET1X')
ax.set_ylabel('DET1Y')

x = mincol*10
y = minrow*10 - 10
width = (maxcol*10) - x #returns 250
height =  (maxrow*10) - y #returns 140

rectangle_region = patches.Rectangle((x-10, y+10), width+20, height, edgecolor = 'green', facecolor ='none')
ax.add_patch(rectangle_region)
plt.show()
#try saving this rectangle as a region file
#use standard DS9 header and style options, rectangle format is (centerx, centery, width, height,angle)
centerx = (width / 2) + mincol*10 #returns 215
centery = (height / 2) + minrow*10 #returns 90

regionfile = open("40111002002A01_SL_Region.reg", "w")
regionfile.write("# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font= 'helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\ndetector\nbox(215,90,250,140,0)")
regionfile.close()

#to read in python:
cat 40111002002A01_SL_Region.reg

#polygon region from rectangles

def corners(centerx, centery, width, height):
    upperleft = (centerx - (width/2) , centery + (height/2))
    upperright = (centerx + (width/2) , centery + (height/2))
    lowerleft = (centerx - (width/2) , centery - (height/2))
    lowerright = (centerx + (width/2), centery - (height/2))
    print(upperleft, upperright, lowerleft, lowerright)

regionfile = open("polygon.reg", "w")
regionfile.write("# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font= 'helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\ndetector\npolygon(90,40,140,40,140,110,190,110,190,150,240,150,240,170,290,170,290,180,340,180,340,120,90,20))")
regionfile.close()

#NEXT STEPS: once we have the region files, we need to know the detector(s) that the region covers to plug into XSPEC
#define function that returns the detector id for all events
#events_indices is a list of arrays. convert to a list of lists and then to one list
def det_ids_for_bin(bin):
    xlower = (bin//36)*10
    xupper = xlower +10
    ylower = (bin%36)*10
    yupper = ylower + 10
    x = tab['DET1X']
    y = tab['DET1Y']
    events_indices = np.argwhere((x > xlower) & (x < xupper) & (y > ylower) & (y < yupper))
    events_indices_list_of_lists = [arr.tolist() for arr in events_indices]
    events_indices_flat_list = []
    for sublist in events_indices_list_of_lists:
        for item in sublist:
            events_indices_flat_list.append(item)
    det_ids = []
    for i in events_indices_flat_list:
        det_ids.append(DET_ID[i])
    return det_ids

#run det_ids_for_bin for all region bins
#all_det_ids a list of arrays, convert to list of lists and then to one list
def all_det_ids():
    all_det_ids = []
    for bin in region_bins:
        det_ids = det_ids_for_bin(bin)
        all_det_ids.append(det_ids)
    all_det_ids_list = []
    for sublist in all_det_ids:
        for item in sublist:
            all_det_ids_list.append(item)
    return all_det_ids_list

all_det_ids_list = all_det_ids()

#running all_det_ids_list.count(i) for i = 0,1,2,3: 85, 0, 21052, 66350

#check that bins not split on detectors
def check_bins():
    all_det_ids = []
    for bin in region_bins:
        det_ids = det_ids_for_bin(bin)
        all_det_ids.append(det_ids)
    all_det_ids_list = []
    check_bins = []
    for sublist in all_det_ids:
        d = sublist[0]
        if sublist.count(d) == len(sublist):
            check_bins.append(True)
        else:
            check_bins.append(False)
    return check_bins

check_bins = check_bins()

bad_bins = [i for i, val in enumerate(check_bins) if not val]
#[84, 85, 87, 88, 90, 91, 92, 93, 94, 95, 96, 97, 278]

def det_ids_for_bad_bins():
    bad_det_ids = []
    for i in bad_bins:
        bin = region_bins[i]
        bad_det_ids.append(det_ids_for_bin(bin))
    return bad_det_ids

det_ids_for_bad_bins = det_ids_for_bad_bins()

from statistics import mode

def mode_bad_bins():
    modes = []
    for list in det_ids_for_bad_bins:
        modes.append(mode(list))
    return modes




#make separate region files for each detector with events on it
def modes():
    all_det_ids = []
    for bin in region_bins:
        det_ids = det_ids_for_bin(bin)
        all_det_ids.append(det_ids)
    modes = []
    for i in all_det_ids:
        modes.append(mode(i))
    return modes

modes = modes()

def det_0_bins():
    det_0_bins = []
    for i in modes:
        if i == 0:
            det_0_bins.append(modes[i])
    return det_0_bins

def det_2_bins():
    det_2_bins = []
    for i in modes:
        if i == 2:
            det_2_bins.append(modes[i])
    return det_2_bins





#next, try making smaller rectangles to include more SL and less background





###SCRATCHWORK###
#rectangle region: obtain 4 corners
def rectangle_xcorners():
    xarray = []
    for i in region_bins:
        xarray.append((indices[i]// 36)*10)
    return xarray
xarray = rectangle_xcorners()


def rectangle_ycorners():
    yarray = []
    for i in region_bins:
        yarray.append((indices[i] % 36)*10)
    return yarray
yarray = rectangle_ycorners()

#FOM rectangle scratch
    lowerleft = (min(xarray)//10)*36 + (min(yarray)//10)
    upperleft = (min(xarray)//10)*36 + (max(yarray)//10)
    lowerright = (max(xarray)//10)*36 + (min(yarray)//10)
    upperright = (max(xarray)//10)*36 + (max(yarray)//10)

mincol = min(region_bins)//36
maxcol = max(region_bins)//36
minrow = min(region_bins)%36
maxrow = max(region_bins)%36

all_bins = range(1,1297,1)
def rectangle():
    rectangle = []
    for i in all_bins:
        if i//36 >= (mincol) and i//36 <= (maxcol) and i%36 >= (minrow-1) and i%36 <= (maxrow):
            rectangle.append(i)
    return rectangle
rectangle = rectangle()

def FOM():
    FOM = 0
    for i in rectangle:
        if i in region_bins:
            FOM += 1
        if i not in region_bins:
            FOM += -1
    return FOM

def bincount():
    bincount = 0
    for i in rectangle:
        if i in region_bins:
            bincount += 1
    return bincount

#find number of region bins per row to find rectangle top and bottom sides
def region_row(x):
    '''returns number of SL bins in given row, argument and output between 0 and 35'''
    bins_in_row = []
    for i in region_bins:
        if i % 36 == x:
            bins_in_row.append(i)
    return len(bins_in_row)

def all_region_rows():
    rows = []
    for x in range(0,36,1):
        rows.append(region_row(x))
    return rows

rows = all_region_rows()

#need to write code to find where values stop beign zero - this is bottom side
#looking at manually for now: bottom row is row 1
#need to write code to automatically find 2nd smallest value (don't care about 0s) - the index of this value is the top row
#looking at manually for now: 2nd smallest value is 1
rows.index(1)
# top row is row 18

#find number of region bins per column to find left and right sides
def region_column(x):
    '''returns number of SL bins in given column, argument and output between 0 and 35'''
    bins_in_column = []
    for i in region_bins:
        if i // 36 == x:
            bins_in_column.append(i)
    return len(bins_in_column)

def all_region_columns():
    columns = []
    for x in range(0,36,1):
        columns.append(region_column(x))
    return columns

columns = all_region_columns()

#first nonzero value is 3
columns.index(9)
#left side is column 9
#last nonzero is 5
columns.index(5)
#right side is column 34

def FOM():
    FOM = 0
    minbin = (36*9) + 1
    maxbin = (36*34) + 18
    rectangle = range(minbin, maxbin, 1)
    for i in rectangle:
        if i in region_bins:
            FOM += 1
        if i not in region_bins:
            FOM += -1
    return FOM





#plot binned image w/rectangle
indices = range(1, 1295, 1)
NBINS=(360,360)
fig,ax = plt.subplots(1)
ax.hist2d(data['DET1X'], data['DET1Y'], NBINS, range=[[1, 360], [1,360]], cmap='gray')
ax.set_xlabel('DET1X')
ax.set_ylabel('DET1Y')
x = min(xarray)
y = min(yarray)
width = (max(xarray)+10) - x
height = (max(yarray)+10) - y
rectangle_region = patches.Rectangle((x, y), width, height, edgecolor = 'green', facecolor ='none')
ax.add_patch(rectangle_region)
plt.show()











### MORE SCRATHWORK###
#create astropy region
from astropy.coordinates import SkyCoord
from astropy import units as u
from regions import PixCoord, RectangleSkyRegion, RectanglePixelRegion
width = (max(xarray)+10) - x
height = (max(yarray)+10) - y
centerx = (width / 2) + min(xarray)
centery = (height / 2) + min(yarray)
rectangle_pix = RectanglePixelRegion(center = PixCoord(centerx, centery), width = width, height = height)

ds9_objects_to_string(regions, coordsys = 'image')
filename = 'region.reg'
write_ds9(region, filename)


#practice saving as DS9 region file
from regions import DS9Parser
width = (max(xarray)+10) - x
height = (max(yarray)+10) - y
centerx = width / 2
centery = height / 2
print(centerx, centery, width, height)
#plug values into string
reg_string = 'image\nbox(130, 90, 260, 180, 0)  color=green'
parser = DS9Parser(reg_string)
print(parser.shapes[0])
regions = parser.shapes.to_regions()
print(regions[0])

from regions import read_ds9, write_ds9
filename = 'SL'
write_ds9(regions, filename)



#save as fits image?
hdu = fits.PrimaryHDU(fig,ax)
hdul = fits.HDUList([hdu])
hdul.writeto('region.fits')

all_bins = range(1,1297,1)

def col1():
    col1 = []
    for i in all_bins:
        if 36//i73 == 1:
            col1.append(i)
    return col1

#automate threshold guessing step
#eventually: eliminate visualizations, should be able to plug in event file and just return region file



#MULTIPLE RECTANGLES REGION
#finding coordinates of all region bins
def region_bin_cols():
    region_bin_cols = []
    for i in region_bins:
        region_bin_cols.append(i//36)
    return region_bin_cols

region_bin_cols = region_bin_cols()

def region_bin_rows():
    region_bin_rows = []
    for i in region_bins:
        region_bin_rows.append(i%36)
    return region_bin_rows

region_bin_rows = region_bin_rows()

'''run max(region_bin_rows) to find upper bound of first rectangle, run region_bin_rows.count(max-x) until value not negligible
for this file, best first upper bound looks like 17 (count = 9)

run region_bin_rows.index(max-x), find entry in regions list with this index: this is the leftmost bin for this rectangle

find column this rectangle is in for rectangle left bound

for this entry, it's bin 845 which is in column 23'''

'''so, the rightmost rectangle goes from columns 23 to 34, rows 1 to 17'''
indices = range(1, 1295, 1)
NBINS=(360,360)
fig,ax = plt.subplots(1)
ax.hist2d(data['DET1X'], data['DET1Y'], NBINS, range=[[1, 360], [1,360]], cmap='gray')
ax.set_xlabel('DET1X')
ax.set_ylabel('DET1Y')
rectangle_region = patches.Rectangle((230, 20), 120, 160, edgecolor = 'green', facecolor ='none')
ax.add_patch(rectangle_region)
plt.show()

'''leftmost bound of this rectangle is rightmost bound of next rectangle. try keeping width. just need new upper bound'''
indices = range(1, 1295, 1)
NBINS=(360,360)
fig,ax = plt.subplots(1)
ax.hist2d(data['DET1X'], data['DET1Y'], NBINS, range=[[1, 360], [1,360]], cmap='gray')
ax.set_xlabel('DET1X')
ax.set_ylabel('DET1Y')
rectangle_region_1 = patches.Rectangle((260, 20), 90, 160, edgecolor = 'green', facecolor ='none')
ax.add_patch(rectangle_region_1)
rectangle_region_2 = patches.Rectangle((170, 20), 90, 150, edgecolor = 'green', facecolor ='none')
ax.add_patch(rectangle_region_2)
rectangle_region_3 = patches.Rectangle((80, 20), 90, 110, edgecolor = 'green', facecolor ='none')
ax.add_patch(rectangle_region_3)
plt.show()

#write FOM for these rectangles
all_bins = range(1,1297,1)
def rectangle1():
    rectangle1 = []
    for i in all_bins:
        if i//36 >= (27) and i//36 <= (35) and i%36 >= (2) and i%36 <= (18):
            rectangle1.append(i)
    return rectangle1
rectangle1 = rectangle1()
def rectangle2():
    rectangle2 = []
    for i in all_bins:
        if i//36 >= (18) and i//36 <= (26) and i%36 >= (2) and i%36 <= (17):
            rectangle2.append(i)
    return rectangle2
rectangle2 = rectangle2()
def rectangle3():
    rectangle3 = []
    for i in all_bins:
        if i//36 >= (8) and i//36 <= (17) and i%36 >= (2) and i%36 <= (13):
            rectangle3.append(i)
    return rectangle3
rectangle3 = rectangle3()


def bincount():
    bincount = 0
    for i in rectangle1:
        if i in region_bins:
            bincount += 1
    for i in rectangle2:
        if i in region_bins:
            bincount += 1
    for i in rectangle3:
        if i in region_bins:
            bincount += 1
    return bincount

def FOM():
    FOM = 0
    for i in rectangle1:
        if i in region_bins:
            FOM += 1
        if i not in region_bins:
            FOM += -1
    for i in rectangle2:
        if i in region_bins:
            FOM += 1
        if i not in region_bins:
            FOM += -1
    for i in rectangle3:
        if i in region_bins:
            FOM += 1
        if i not in region_bins:
            FOM += -1
    return FOM
