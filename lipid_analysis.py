#import scipy
import numpy as np
import sys

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

def getFrameDistance(bottom_leaflet, top_leaflet, lipids_per_leaflet):
    
    top_pos = 0
    bottom_pos = 0
    #print(top_leaflet[0,0])
    for row in top_leaflet:
        top_pos += row[2]
    for row in bottom_leaflet:
        bottom_pos += row[2]   
    top_pos = top_pos / lipids_per_leaflet
    bottom_pos = bottom_pos / lipids_per_leaflet
            
    distance = top_pos - bottom_pos
       
    return top_pos, bottom_pos, distance




def bilayerDistance(file, head_atom_name, lipids_per_leaflet, first_frame, last_frame):
    arc = open(file, 'r')
    
    # Get first line to identify splits between frames
    first_line = arc.readline().rstrip()
    
    # the number of frames we have read through
    frame_num = 0
        
    frames_to_process = last_frame - first_frame + 1
    distances = np.empty((frames_to_process, 3), dtype=float)
    
    top_leaflet = np.empty((lipids_per_leaflet, 3), dtype=float)
    bottom_leaflet = np.empty((lipids_per_leaflet, 3), dtype=float)
    
    arc = open(file, 'r')
    for line in arc:
        cleaned_line = line.rstrip()
        tokens = cleaned_line.split()
        
        # indicates we have found a new frame
        if cleaned_line == first_line:
            #print("found new frame") 
            
            # Process last frame
            if frame_num > 0 and frame_num >= first_frame and frame_num <= last_frame:
                top_pos, bottom_pos, distance = getFrameDistance(bottom_leaflet, top_leaflet, lipids_per_leaflet)
                print("Frame Number: " + str(frame_num) + ", Distance: " + str(distance))
                distances[frame_num-1, 0] = top_pos
                distances[frame_num-1, 1] = bottom_pos
                distances[frame_num-1, 2] = distance
            
            # iterate frame number counter
            frame_num += 1
             
            # start new frame
            top_leaflet = np.empty((lipids_per_leaflet, 3), dtype=float)
            bottom_leaflet = np.empty((lipids_per_leaflet, 3), dtype=float)
            
            num_top_leaflet = 0
            num_bottom_leaflet = 0    
        
        # indicates line contains a valid atom that we should save to our frame data       
        elif len(tokens) > 5 and tokens[1] == head_atom_name:
            
            x = float(tokens[2])
            y = float(tokens[3])
            z = float(tokens[4])
            
            if z > 0:
                #print(frame_num, z)
                top_leaflet[num_top_leaflet, 0] = x
                top_leaflet[num_top_leaflet, 1] = y
                top_leaflet[num_top_leaflet, 2] = z
                num_top_leaflet += 1
            else:
                bottom_leaflet[num_bottom_leaflet, 0] = x
                bottom_leaflet[num_bottom_leaflet, 1] = y
                bottom_leaflet[num_bottom_leaflet, 2] = z
                num_bottom_leaflet += 1
                
    return distances


file = sys.argv[1]
lipids_per_leaflet = int(sys.argv[2])
num_frames = getFrameCount(file)
print("Num Frames = " + str(num_frames))
distances = bilayerDistance(file, "P1", lipids_per_leaflet, 1, num_frames)

avg_dist = 0
i = 0
for row in distances:
    avg_dist += row[2]
    i += 1
avg_dist = avg_dist / i

print("Average distance = " + str(avg_dist))
print("Std. Dev. = " + str(np.std(distances[:,2])))
            
                
                
            
        
        

