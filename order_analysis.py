# -*- coding: utf-8 -*-
import numpy as np
import sys


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


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
        
        
def getSingleFrameParam(atoms):
    angleSum = 0
    numAngles = 0
    normal = np.array([0,0,1])
    
    # iterate over all atoms
    for atom in atoms:
        
        # reset params
        hAtomNum = 0
        num_C = 0
        num_H = 0
        
        # check if atom is valid
        if "C" in atoms[atom].name and len(atoms[atom].bound_atoms) == 4:
            #print("Found carbon with 4 bonds")
            # print(atoms[atom].bound_atoms)
            # look through bonds in more detail
            for bonded_atom in atoms[atom].bound_atoms:
                #print(atoms[bonded_atom].name)
                if "C" in atoms[bonded_atom].name:
                    #print("found C")
                    num_C += 1
                elif "H" in atoms[bonded_atom].name:
                    #print("found H")
                    num_H += 1
                    hAtomNum = atoms[bonded_atom].num
                    
            # if still valid
            if num_C == 2 and num_H == 2:
                #print("also has 2C and 2H bound")
                v1 = np.array([atoms[atom].coords[0] - atoms[hAtomNum].coords[0], atoms[atom].coords[1] - atoms[hAtomNum].coords[1], atoms[atom].coords[2] - atoms[hAtomNum].coords[2]])
                angle = angle_between(normal, v1)
                #print("Angle = " + str(angle))
                angleSum += angle
                numAngles += 1
    
    return angleSum / numAngles
                


def orderParams(file, first_frame, last_frame):
    arc = open(file, 'r')
    
    # Get first line to identify splits between frames
    first_line = arc.readline().rstrip()
    
    # the number of frames we have read through
    frame_num = 0
        
    frames_to_process = last_frame - first_frame + 1
    angles = np.empty((frames_to_process, 1), dtype=float)
   
    
    atoms = {}
    
    arc = open(file, 'r')
    for line in arc:
        cleaned_line = line.rstrip()
        tokens = cleaned_line.split()
        
        # indicates we have found a new frame
        if cleaned_line == first_line:
            #print("found new frame") 
            
            # Process last frame
            if frame_num > 0 and frame_num >= first_frame and frame_num <= last_frame:
                avg_order = getSingleFrameParam(atoms)
                angles[frame_num] = avg_order
                print("Frame Number: " + str(frame_num) + ", Order: " + str(avg_order))
                
            
            # iterate frame number counter
            frame_num += 1
             
            # start new frame
            atoms.clear()   
        
        # indicates line contains a valid atom that we should save to our frame data       
        elif len(tokens) > 5:
            
            newAtom = return_atom(tokens)
            atoms[newAtom.num] = newAtom
            

file = sys.argv[1]

num_frames = getFrameCount(file)
print("Num Frames = " + str(num_frames))
avg_order = orderParams(file, 1, num_frames)