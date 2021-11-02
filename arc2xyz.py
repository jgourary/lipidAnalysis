import sys


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


class Atom:
    def __init__(self, num, name, coords, atom_type, bound_atoms):
        self.name = name
        self.num = num
        self.coords = coords
        self.atom_type = atom_type
        self.bound_atoms = bound_atoms


# Translate tokens to an atom data type
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


# convert an atom data type to a string
def print_atom(atom):
    result = str(atom.num) + " " + atom.name + " " + str(atom.coords[0]) + " " + str(atom.coords[1])
    result += " " + str(atom.coords[2]) + " " + str(atom.atom_type)
    for bound_atom in atom.bound_atoms:
        result += " " + str(bound_atom)
    result += "\n"
    return result

# open arc file
file1 = open(sys.argv[1], 'r')

# open output XYZ
file2 = open(sys.argv[2], 'w+')

lineNum = 1
header1 = file1.readline().rstrip()
header2 = ""
atoms = []
break_next_line = False
for line in reversed(file1.readlines()):
    tokens = line.split()

    if len(tokens) == 6 and isfloat(tokens[3]) and isfloat(tokens[4]) and isfloat(tokens[5]) and int(float(tokens[3])) == 90 and int(float(tokens[4])) == 90 and int(float(tokens[5])) == 90:
        print("breaking at " + str(lineNum) + ": " + line)
        header2 = line
        break_next_line = True
    elif break_next_line:
        header1 = line
        break
    elif len(tokens) >= 6:
        newAtom = return_atom(tokens)
        atoms.append(newAtom)

    lineNum += 1

file2.write(header1)
print(header1)
file2.write(header2)
print(header2)
atoms.reverse()
for atom in atoms:
    file2.write(print_atom(atom))

