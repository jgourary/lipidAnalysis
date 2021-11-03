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
        
 
def getBackBoneAtoms(file):
    xyz = open(file, 'r')
    atoms = {}
    valid_id = []
    
    lineNum = 1
    for line in xyz:
        cleaned_line = line.rstrip()    
        tokens = cleaned_line.split()
        
        if len(tokens) > 5 and lineNum > 2:            
            newAtom = return_atom(tokens)
            atoms[newAtom.num] = newAtom
            
            
        lineNum += 1
        
    for atom in atoms:
        if "C" in atoms[atom].name and len(atoms[atom].bound_atoms) == 4:
            num_C = 0
            num_H = 0
            for bonded_atom in atoms[atom].bound_atoms:
                #print(atoms[bonded_atom].name)
                if "C" == atoms[bonded_atom].name[0]:
                    #print("found C")
                    num_C += 1
                elif "H" == atoms[bonded_atom].name[0]:
                    #print("found H")
                    num_H += 1
                    
            # if still valid
            if num_C == 2 and num_H == 2:
                valid_id.append(atom)
    print(valid_id)            
    return valid_id
            
            
def getSingleFrameParam(atoms):
    angleSum = 0
    orderSum = 0
    numAngles = 0
    normal = np.array([0,0,1])
    
    # iterate over all atoms
    for atom in atoms:
        
        # reset params
        hAtomNum = 0
        
        # check if atom is valid
        if atom in valid_id:
            
            # look through bonds in more detail
            for bonded_atom in atoms[atom].bound_atoms:
                if "H" in atoms[bonded_atom].name:
                    hAtomNum = atoms[bonded_atom].num
                    
            v1 = np.array([atoms[atom].coords[0] - atoms[hAtomNum].coords[0], atoms[atom].coords[1] - atoms[hAtomNum].coords[1], atoms[atom].coords[2] - atoms[hAtomNum].coords[2]])
            angle = angle_between(normal, v1)
            order = 3 * (np.cos(angle))**2 - 1
            #print(angle)
            #print(order)
            
            angleSum += angle
            orderSum += order
            numAngles += 1
                
    
    angleSum = angleSum / numAngles
    orderSum = 0.5 * orderSum / numAngles
    #print(angleSum)
    return angleSum, orderSum
                


def orderParams(file, first_frame, last_frame):
    arc = open(file, 'r')
    
    # Get first line to identify splits between frames
    first_line = arc.readline().rstrip()
    
    # the number of frames we have read through
    frame_num = 0
        
    frames_to_process = last_frame - first_frame + 1
    angles = np.empty((frames_to_process, 2), dtype=float)
   
    
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
                avg_angle, avg_order = getSingleFrameParam(atoms)
                angles[frame_num,0] = avg_angle
                angles[frame_num,1] = avg_order
                print("Frame Number: " + str(frame_num) + ", Angle: " + str(avg_angle) + ", Order: " + str(avg_order))
                
            
            # iterate frame number counter
            frame_num += 1
             
            # start new frame
            atoms.clear()   
        
        # indicates line contains a valid atom that we should save to our frame data       
        elif len(tokens) > 5:
            
            newAtom = return_atom(tokens)
            atoms[newAtom.num] = newAtom
    return angles
            

file = sys.argv[1]
ref = sys.argv[2]

valid_id = getBackBoneAtoms(ref)

num_frames = getFrameCount(file)
print("Num Frames = " + str(num_frames))

angles = orderParams(file, 1, num_frames)

final_angle = 0
final_order = 0
i = 0
for row in angles:
    final_angle += row[0]
    final_order += row[1]
    i += 1
final_angle = final_angle / i
final_order = final_order / i

print("Average Angle = " + str(final_angle))
print("Average Order = " + str(final_order))