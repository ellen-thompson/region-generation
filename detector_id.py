#define function that returns detector ids for events in a given bin:
def events_DET1X(b):
    bin = region_bins[b]
    xlower = (bin//36)*10
    xupper = xlower +10
    x = tab['DET1X']
    events_DET1X = []
    for i in x:
        if i >= xlower and i <= xupper:
            events_DET1X.append(i)
    return events_DET1X

def events_DET1Y(b):
    bin = region_bins[b]
    ylower = (bin%36)*10
    yupper = ylower + 10
    y = tab['DET1Y']
    events_DET1Y = []
    for i in y:
        if i >= ylower and i <= yupper:
            events_DET1Y.append(i)
    return len(events_DET1Y)

#ex: region bin 200, 250<x<260, 70<y<80
def ids(xlower, xupper, ylower, yupper):
    x = tab['DET1X']
    y = tab['DET1Y']
    events = np.where((x > xlower) & (x < xupper) & (y > ylower) & (y < yupper))
    return (np.size(events))

#define function that returns the detector id for all events
#events_indices is a list of arrays. convert to a list of lists and then to one list
def det_ids_for_bin(b):
    bin = region_bins[b]
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

##SCRATWHCWORK###
events_in_bin_200 = events_in_bin()
l = [arr.tolist() for arr in events_in_bin_200]
flat_list = []
for sublist in l:
    for item in sublist:
        flat_list.append(item)

def det_ids():
    det_ids = []
    for i in flat_list:
        det_ids.append(DET_ID[i])
    return det_ids

det_ids = det_ids()



#BAD bins
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
    return check_binsregionbin84
