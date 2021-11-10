# -*- coding: utf-8 -*-
import numpy as np
import sys


# Parameters
verbose = False
num_bins = 100 # bin width = box height / num_bins


def getAtomicNumbers():
    file = open(prm, 'r')
    masses = dict()

    for line in file:
        cleaned_line = line.rstrip()
        tokens = cleaned_line.split()
        if len(tokens) > 5 and tokens[0] == "atom":
            # remove atom name in quotes
            tokens = cleaned_line.split("\"")
            del(tokens[1])
            cleaned_line = ''.join(tokens)
            cleaned_line = cleaned_line.replace("\"", "")

            # resplit into tokens
            tokens = cleaned_line.split()
            if tokens[0] == "atom" and len(tokens) > 5:
                masses[tokens[1]] = float(tokens[4])
                if verbose:
                    print(tokens[1] + ", " + str(masses[tokens[1]]))
    print("Found atomic numbers for " + str(len(masses)) + " atom types")
    return masses


def getFrameCount():
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

    print("Found " + str(num_frames) + " frames in the archive")

    return num_frames


class Atom:
    def __init__(self, num, name, coords, atom_type, bound_atoms):
        self.name = name
        self.num = num
        self.coords = coords
        self.atom_type = atom_type
        self.bound_atoms = bound_atoms


class frameData:
    def __init__(self, frame_count, electrons, densities):
        self.frame_count = frame_count
        self.electrons = electrons
        self.densities = densities


def get_box_dimensions():
    arc = open(file, 'r')
    i = 0
    for line in arc:
        if i == 1:
            tokens = line.split()
            return [float(tokens[0]), float(tokens[1]), float(tokens[2])]
        i += 1

    return -1, -1, -1


def get_bin_ranges():
    print("Finding bin cutoffs")
    dimensions = get_box_dimensions()
    print("Box Dimensions: [" + str(dimensions[0]) + ", " + str(dimensions[1]) + ", " + str(dimensions[2]) + "] Angstroms")
    bin_width = dimensions[2] / num_bins
    print("Setting bin width to " + str(bin_width) + " Angstroms")
    lower_bound = dimensions[2] / -2
    print("Lower bound is " + str(lower_bound) + " Angstroms")
    bin_cutoffs = np.empty((num_bins, 1))
    for i in range(0, num_bins):
        bin_cutoffs[i] = lower_bound + (i+1) * bin_width  # this is the upper cutoff

    bin_volume = bin_width * dimensions[0] * dimensions[1]
    print("Setting bin volume to " + str(bin_volume) + " A^3")
    return bin_cutoffs, bin_volume


def get_atom_bin(z):
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
            
            
def single_frame_electron_density(frame_num, atoms):
    electrons = np.zeros((num_bins, 1))
    densities = np.zeros((num_bins, 1))

    # assign all atoms to relevant bins
    if verbose:
        print("Beginning atom bin assignment on frame " + str(frame_num))
    for atom_id in atoms:
        #  if verbose:
        #    print("Placing atom " + str(atoms[atom_id].num) + "...")
        i = 0
        for cutoff in bin_cutoffs:
            if atoms[atom_id].coords[2] < cutoff:
                electrons[i] += atomic_numbers[atoms[atom_id].atom_type]
                #if verbose:
                #    print("Placing atom " + str(atoms[atom_id].num) + " in bin " + str(i))
                break

            i += 1

    #  convert bin electrons to bin electron density
    i = 0
    for elec_count in electrons:
        densities[i] = elec_count / bin_volume
        if verbose:
            print(densities[i])
        i += 1
    newFrame = frameData(frame_num, electrons, densities)
    return newFrame


def electron_density(first_frame, last_frame, frame_interval):
    arc = open(file, 'r')
    out_file = open(out, 'w')
    # Get first line to identify splits between frames
    first_line = arc.readline().rstrip()
    
    # the number of frames we have read through
    frame_num = 0

    all_frame_data = []
    #  frames_to_process = last_frame - first_frame + 1
    atoms = {}
    
    arc = open(file, 'r')
    for line in arc:
        cleaned_line = line.rstrip()
        tokens = cleaned_line.split()

        # indicates we have found a new frame
        if cleaned_line == first_line:
            # Process last frame
            if verbose:
                print("Found frame no. " + str(frame_num))
            if 1 <= frame_num <= last_frame and frame_num >= first_frame and frame_num % frame_interval == 0:
                print("Finding densities for frame " + str(frame_num))
                this_frame_data = single_frame_electron_density(frame_num, atoms)
                all_frame_data.append(this_frame_data)

            # iterate frame number counter
            frame_num += 1
             
            # start new frame
            atoms.clear()   
        
        # indicates line contains a valid atom that we should save to our frame data       
        elif len(tokens) > 5:
            
            newAtom = return_atom(tokens)
            atoms[newAtom.num] = newAtom
    # Process last frame
    this_frame_data = single_frame_electron_density(frame_num, atoms)
    all_frame_data.append(this_frame_data)

    average_densities = np.zeros((num_bins, 1), dtype=frameData)
    for frame in all_frame_data:
        i = 0
        for density in frame.densities:
            if verbose:
                print(density)
            average_densities[i] += density
            i += 1
    i = 0
    for density in average_densities:
        average_densities[i] = density / len(all_frame_data)
        i += 1

    i = 0
    for cutoff in bin_cutoffs:
        out_file.write(str(cutoff[0]) + ", " + str(average_densities[i,0]) + "\n")
        i += 1

    return average_densities
            

file = sys.argv[1]
prm = sys.argv[2]
out = sys.argv[3]
frame_interval = 10

bin_cutoffs, bin_volume = get_bin_ranges()
atomic_numbers = getAtomicNumbers()
num_frames = getFrameCount()
bins = electron_density(1, num_frames, frame_interval)

