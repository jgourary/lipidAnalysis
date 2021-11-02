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
    this_coords = [tokens[2], tokens[3], tokens[4]]
    this_atom_type = tokens[5]
    if len(tokens) >= 7:
        this_bound_atoms = tokens[6:]
    else:
        this_bound_atoms = []
    newAtom = Atom(this_num, this_name, this_coords, this_atom_type, this_bound_atoms)
    return newAtom


def print_atom(atom):
    result = str(atom.num) + " " + atom.name + " " + str(atom.coords[0]) + " " + str(atom.coords[1])
    result += " " + str(atom.coords[2]) + " " + str(atom.atom_type)
    for bound_atom in atom.bound_atoms:
        result += " " + str(bound_atom)
    result += "\n"
    return result


def return_type(atom, lipid_size, lipid_type_dict):
    if atom.name == "POT":
        return 8
    elif atom.name == "CLA":
        return 15
    elif atom.name == "SOD":
        return 7
    elif atom.name == "OT":
        return 36
    elif atom.name == "HT":
        return 37
    else:
        numInLipid = int(atom.num) % lipid_size
        if numInLipid == 0:
            numInLipid = str(lipid_size)
        else:
            numInLipid = str(numInLipid)
        #print(numInLipid)
        return lipid_type_dict[numInLipid]


file1 = open(r'C:\Users\jtgou\OneDrive\Documents\UT_Austin\ren_lab\lipidAnalysis\lipid_references\dmpg.xyz', 'r')

lipid_size = 111

file2 = open(r'C:\Users\jtgou\OneDrive\Documents\UT_Austin\ren_lab\lipidAnalysis\bilayers\dmpg_production\step5_input.xyz', 'r')

file3 = open(r'C:\Users\jtgou\OneDrive\Documents\UT_Austin\ren_lab\lipidAnalysis\bilayers\dmpg_production\step5_output.xyz', 'w+')


# Parse model file
print("Reading model compound")
lines = file1.readlines()
lineNum = 1
lipid_type_dict = dict()
for line in lines:
    #print(line)
    tokens = line.split()
    if lineNum > 2 and len(tokens) >= 6:
        #print(line)
        newAtom = return_atom(tokens)
        if int(newAtom.num) <= lipid_size:
            lipid_type_dict[newAtom.num] = newAtom.atom_type
    lineNum += 1

# Parse input XYZ
print("Translating input XYZ")
lines = file2.readlines()
lineNum = 1
for line in lines:
    tokens = line.split()
    if lineNum <= 2:
        file3.write(line)
    elif len(tokens) >= 6:
        newAtom = return_atom(tokens)
        newAtom.atom_type = return_type(newAtom, lipid_size, lipid_type_dict)
        file3.write(print_atom(newAtom))
    lineNum += 1


