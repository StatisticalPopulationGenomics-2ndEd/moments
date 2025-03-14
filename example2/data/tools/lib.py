## utilities used by other scripts

import gzip
import numpy as np


def regions_to_boolmask(regions, l=None):
    """
    Transform an array of regions into a boolean mask where `1` denotes sites
    outside the mask and `0` sites within it.
    """
    if l is None:
        l = regions[-1, 1]
    boolmask = np.ones(l, dtype=bool)
    for start, end in regions:
        if start >= l:
            break
        boolmask[start:end] = False

    return boolmask


def boolmask_to_regions(boolmask):
    """
    Transform a boolean mask where 1 is 'masked out' to an array of mask regions
    """
    jumps = np.diff(np.concatenate(([1], boolmask, [1])))
    starts = np.where(jumps == -1)[0]
    ends = np.where(jumps == 1)[0]
    regions = np.stack([starts, ends], axis=1)

    return regions


def read_bedfile(file):

    openfunc = gzip.open if file.endswith(".gz") else open 
    regions = []
    with openfunc(file, "rb") as fin:
        for lineb in fin:
            line = lineb.decode()
            if "#" in line:
                continue
            chromnum, start, end = line.split()[:3]
            regions.append([int(start), int(end)])
    region_arr = np.array(regions, dtype=np.int64)

    return region_arr, chromnum


def write_bedfile(file, regions, chromnum):
    """
    write a .bed file (supports .gz) from an array of mask regions.
    """
    openfunc = gzip.open if file.endswith(".gz") else open 
    with openfunc(file, "wb") as fout:
        for start, end in regions:
            lineb = f"{chromnum}\t{start}\t{end}\n".encode()
            fout.write(lineb)
    return


def intersect_regions(region_arrs):
    """
    Get a region array representing the intersection of sites in several input
    region arrays.
    """
    l = np.max([regs[-1, 1] for regs in region_arrs])
    counter = np.zeros(l, dtype=np.uint8)
    for regions in region_arrs:
        boolmask = regions_to_boolmask(regions, l=l)
        counter[~boolmask] += 1
    isec_mask = counter < len(region_arrs)  # recall that 1 is `masked`
    isec_regions = boolmask_to_regions(isec_mask)

    return isec_regions


def union_regions(region_arrs):
    """
    Obtain a region array representing the union of sites in input region arrays
    """
    l = np.max([regs[-1, 1] for regs in region_arrs])
    counter = np.zeros(l, dtype=np.uint8)
    for regions in region_arrs:
        boolmask = regions_to_boolmask(regions, l=l)
        counter[~boolmask] += 1
    union_mask = counter == 0  # recall that 1 is `masked`
    union_regions = boolmask_to_regions(union_mask)

    return union_regions


def flank_regions(region_arr, flank):
    """
    Extend region starts and ends by `flank` base pairs.
    """
    starts = region_arr[:, 0]  - flank 
    starts[starts < 0] = 0
    ends = region_arr[:, 1] + flank
    flanked = np.stack((starts, ends), axis=1)
    resolved = boolmask_to_regions(regions_to_boolmask(flanked))

    return resolved


def subtract_regions(regions, subtrahend):
    """
    Get an array of regions representing sites which belong to `regions` and
    not to `subtrahend`.
    """
    l = max(regions[-1, 1], subtrahend[-1, 1])
    boolmask = regions_to_boolmask(regions, l=l)
    submask = regions_to_boolmask(subtrahend, l=l)
    indicator = ~submask 
    counter = boolmask + indicator
    outmask = counter == 1
    ret = boolmask_to_regions(outmask)

    return ret
