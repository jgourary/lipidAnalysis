# -*- coding: utf-8 -*-
import numpy as np
import sys


# Parameters
verbose = False
num_bins = 100 # bin width = box height / num_bins


def getAtomicNumbers(file):
    prm = open(file, 'r')
    masses = dict()

    for line in prm:
        cleaned_line = line.rstrip()
        tokens = cleaned_line.split()
        if tokens[0] == "atom":
            # remove atom name in quotes
            tokens = cleaned_line.split("\"")
            del(tokens[1])
            cleaned_line = ''.join(tokens)
            cleaned_line = cleaned_line.replace("\"", "")

            # resplit into tokens
            tokens = cleaned_line.split()
            if tokens[0] == "atom" and len(tokens) > 5:
                masses[tokens[1]] = float(tokens[5])

    return masses


def getFrameCount(file):
    arc = open(file, 'r')

    num_frames = 0

    # Get first line to identify splits between frames
    first_line = arc.readline().rstrip()

    arc = open(file, 'r')
    for line in arc:
        cleaned_line = line.rstrip()

        # indicates we have found a new frame
        if cleaned_line == first_line:
            num_frames += 1

    return num_frames


class Atom:
    def __init__(self, num, name, coords, atom_type, bound_atoms):
        self.name = name
        self.num = num
        self.coords = coords
        self.atom_type = atom_type
        self.bound_atoms = bound_atoms


def get_box_dimensions(file):
    xyz = open(file, 'r')
    for line in xyz:
        if "90.00000000   90.00000000   90.00000000" in line:
            tokens = line.split()
            return [float(tokens[0]), float(tokens[1]), float(tokens[2])]

    return -1, -1, -1


def get_bin_ranges(file):
    dimensions = get_box_dimensions(file)
    bin_width = dimensions[2] / num_bins
    lower_bound = dimensions[2] / -2
    bin_cutoffs = np.empty((num_bins, 1))
    for i in range(0, num_bins):
        bin_cutoffs[i] = lower_bound + (i+1) * bin_width  # this is the upper cutoff
    bin_volume = bin_width * dimensions[0] * dimensions[1]
    return bin_cutoffs, bin_volume


def get_atom_bin(bin_cutoffs, z):
    i = 0
    for bin_cutoff in bin_cutoffs:
        if z < bin_cutoff:
            return i
        i += 1
    return -1


def return_atom(tokens):
    this_num = tokens[0]
    this_name = tokens[1]
    this_coords = [float(tokens[2]), float(tokens[3]), float(tokens[4])]
    this_atom_type = tokens[5]
    if len(tokens) >= 7:
        this_bound_atoms = tokens[6:]
    else:
        this_bound_atoms = []
    newAtom = Atom(this_num, this_name, this_coords, this_atom_type, this_bound_atoms)
    return newAtom
            
            
def single_frame_electron_density(atoms):
    single_frame_bins = np.empty((num_bins, 1))
    
    for atom_id in atoms:
        i = 0
        for cutoff in bin_cutoffs:
            if atoms[atom_id].coords[2] < cutoff:
                single_frame_bins[i] += atomic_numbers[atoms[atom_id].atom_type]
            i += 1
    i = 0
    for bin in single_frame_bins:
        single_frame_bins[i] = bin / bin_volume
        i += 1
    return single_frame_bins


def electron_density(first_frame, last_frame):
    arc = open(file, 'r')
    out = open(file, 'w')
    out.write(bin_cutoffs)
    # Get first line to identify splits between frames
    first_line = arc.readline().rstrip()
    
    # the number of frames we have read through
    frame_num = 0

    average_bins = np.empty((num_bins, 1), dtype=float)
    #  frames_to_process = last_frame - first_frame + 1
    atoms = {}
    
    arc = open(file, 'r')
    for line in arc:
        cleaned_line = line.rstrip()
        tokens = cleaned_line.split()
        
        # indicates we have found a new frame
        if cleaned_line == first_line:

            # Process last frame
            if 0 < frame_num <= last_frame and frame_num >= first_frame:
                single_frame_bins = single_frame_electron_density(atoms)
                out.write(single_frame_bins)
                for i in (0, len(single_frame_bins)):
                    average_bins[i] += single_frame_bins[i]

                if verbose:
                    print("Frame Number: " + str(frame_num))

            # iterate frame number counter
            frame_num += 1
             
            # start new frame
            atoms.clear()   
        
        # indicates line contains a valid atom that we should save to our frame data       
        elif len(tokens) > 5:
            
            newAtom = return_atom(tokens)
            atoms[newAtom.num] = newAtom

        for i in (0, len(average_bins)):
            average_bins[i] = average_bins[i] / bin_volume
            print("Bin Cutoff = " + str(bin_cutoffs[i]) + ", Density = " + str(average_bins[i]))
    return average_bins
            


file = sys.argv[1]
ref = sys.argv[2]
prm = sys.argv[3]
out = sys.argv[4]

bin_cutoffs, bin_volume = get_bin_ranges(file)
atomic_numbers = getAtomicNumbers(prm)
num_frames = getFrameCount(file)
print("Num Frames = " + str(num_frames))
bins = electron_density(file, 1, num_frames)

