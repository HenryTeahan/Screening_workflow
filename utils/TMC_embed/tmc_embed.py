from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import numpy as np
from itertools import combinations
import subprocess
import sys
import itertools
import os
import shutil
from rdkit.Geometry import Point3D
import py3Dmol

#sys.path.append("/home/maria/github/RMSD-PP25/HelperFunctions")
#sys.path.append("/home/maria/github/autoTS")
from .utils import draw3D, reorder_atoms_to_map, run_cmd, embed_smiles_far

#from pre_align_new import embed_smiles_far
pt = Chem.GetPeriodicTable()
TRANSITION_METALS_NUM = [21,22,23,24,25,26,27,57,28,29,30,39,40,41,42,43,44,45,46,47,48,71,72,73,74,75,76,77,78,79,80]
geom_dict = dict() #define default geometries
geom_dict[5] = "trigonal_bipyramidal"
geom_dict[4] = "square_planar"
geom_dict[3] = "trigonal_planar"
geom_dict[2] = "linear"
geom_dict[6] = "octahedral"


angle_dict = {
    "linear": [(0,1,180)],
    "trigonal_planar": [(0,1,120), (1,2,120)],
    "square_planar": [(0,2,180), (1,3,180), (0,1,90)],
    "trigonal_bipyramidal": [(0,1,180), (1,2,90), (0,3,90), (2,3,120), (3,4,120), (2,4,120)],
    "octahedral": [(0,1,180), (2,3,180), (4,5,180), (0,2,90), (1,4,90), (3,5,90)]
}

coord_permutations = {
    "square_planar": [
        [0,1,2,3],
        [0,2,1,3],
        [0,2,3,1]
    ],
    "trigonal_bipyramidal": [
        [0,1,2,3,4],
        [0,2,1,3,4],
        [0,3,1,2,4],
        [0,4,1,2,3],
        [1,2,0,3,4],
        [1,3,0,2,4],
        [1,4,0,2,3],
        [2,3,0,1,4],
        [2,4,0,1,3],
        [3,4,0,1,2]
    ],
    "octahedral": [
        [0,1,2,3,4,5],
        [0,1,2,4,3,5],
        [0,1,2,5,3,4],
        [0,2,1,3,4,5],
        [0,2,1,4,3,5],
        [0,2,1,5,3,4],
        [0,3,2,1,4,5],
        [0,3,2,4,1,5],
        [0,3,2,5,1,4],
        [0,4,2,3,1,5],
        [0,4,2,1,3,5],
        [0,4,2,5,3,1],
        [0,5,2,3,4,1],
        [0,5,2,4,3,1],
        [0,5,2,1,3,4]
    ]
}

def visualize_xyz_file(xyz_file):
    with open(xyz_file, "r") as f:
        xyz_data = f.read()
    #print(xyz_data)

    # Create a 3Dmol view
    view = py3Dmol.view(width=400, height=400)
    view.addModel(xyz_data, 'xyz')  # specify format as 'xyz'
    view.setStyle({'stick': {}})    # you can also use 'sphere', 'line', etc.
    view.zoomTo()
    view.show()

def setup_angle_constraints(geometry, coord_order, coord_dir, tmc_idx):
    angles = angle_dict[geometry]
    constraints = []
    for angle in angles:
        a1,a2,d = angle
        a1 = int(coord_order[a1])
        a2 = int(coord_order[a2])
        for element in itertools.product(*[coord_dir[a1], coord_dir[a2]]):
            aa1, aa2 = element

            constraints.append((aa1-1, int(tmc_idx), aa2-1, d))

    #print(constraints)
    return constraints

def setup_constrain_dict_pos(mhap, coord_order, coord_dir, geometry=None):
    """
    coord_order should contain the coordinating atoms in the order of position in the complex:
    metal always (0,0,0)
    square_planar: [(0,0,1),(0,1,0),(0,0,-1), (0,-1,0)]
    etc.
    """
    constraints = []
    tmc_idx = None
    mol = rdmolops.HapticBondsToDative(mhap)
    for a in mol.GetAtoms():
        if a.GetAtomicNum() in TRANSITION_METALS_NUM:
            tmc_idx = a.GetIdx()
    #print(tmc_idx)
    #coordinating_atoms = np.nonzero(Chem.rdmolops.GetAdjacencyMatrix(mol)[tmc_idx, :])[0]
    #print(coordinating_atoms)

    #mhap = rdmolops.DativeBondsToHaptic(mol)
    N_binding_sites = len(mhap.GetAtoms()[tmc_idx].GetBonds())
    if not geometry:
        geometry = geom_dict[N_binding_sites]
    print(f"using {geometry} as geometry template")


    angle_constraints = setup_angle_constraints(geometry, coord_order, coord_dir, tmc_idx)
    emol = Chem.RWMol(mol)


    for i in coord_dir:
        for j in coord_dir[i]:
            j -= 1
            atom = mol.GetAtomWithIdx(int(j))
            atomic_num = atom.GetAtomicNum()
            bond = mol.GetBondBetweenAtoms(int(j), tmc_idx)
            if bond:
                emol.RemoveBond(int(j), tmc_idx)
                ri = pt.GetRcovalent(atomic_num)
                rj = pt.GetRcovalent(mol.GetAtomWithIdx(tmc_idx).GetAtomicNum())
                constraints.append((int(j), tmc_idx, 0.95*(ri+rj)))

    return set(constraints), set(angle_constraints), emol.GetMol()

def write_constrain_file(constraints, constraints_angle, fconstant):
    constrain_str = ""
    for i, j, value in set(constraints):
        constrain_str += f"\n    {"distance"}: {', '.join([str(x) for x in [i+1,j+1]])}, {value}"
    for i, j, k, value in set(constraints_angle):
        #print(i,j,k,value)
        constrain_str += f"\n    {"angle"}: {', '.join([str(x) for x in [i+1,j+1,k+1]])}, {value}"

    detailed_input = f"""
    $constrain
        force constant={fconstant}
        {constrain_str}
    $end
    """.format()

    with open("constrain.inp", 'w') as _file:
        _file.write(detailed_input)

def assign_coordinates2mol(mol, xyz_file, constraints):
    atomic_symbols = []
    xyz_coordinates = []

    with open(xyz_file, "r") as file:
        for line_number,line in enumerate(file):
            if line_number == 0:
                num_atoms = int(line)
            elif line_number == 1:
                comment = line # might have useful information
            else:
                atomic_symbol, x, y, z = line.split()
                atomic_symbols.append(atomic_symbol)
                xyz_coordinates.append([float(x),float(y),float(z)])

    conf = mol.GetConformer()

    for i in range(mol.GetNumAtoms()):
        x,y,z = xyz_coordinates[i]
        conf.SetAtomPosition(i,Point3D(x,y,z))

    emol = Chem.RWMol(mol)
    for i,j,_ in set(constraints):
        emol.AddBond(i, j, order=Chem.rdchem.BondType.DATIVE)


    return emol.GetMol()

def get_hapt_coord_dir(mhap, coord_order):

    tmc_idx = None
    for a in mhap.GetAtoms():
        if a.GetAtomicNum() in TRANSITION_METALS_NUM:
            tmc_idx = a.GetIdx()

    coord_dir = {}
    for i,bond in enumerate(mhap.GetAtoms()[tmc_idx].GetBonds()):
        i = int(i)
        batom = bond.GetBeginAtom()
        eatom = bond.GetEndAtom()
        #print(batom.GetIdx()+1, eatom.GetIdx()+1)
        if batom.GetSymbol() == "*":
            hapt_atoms = bond.GetPropsAsDict()['_MolFileBondEndPts']
            hapt_atoms = [int(x) for x in hapt_atoms.strip('()').split()[1:]]
            if not coord_order:
                coord_dir[hapt_atoms[0]] = hapt_atoms
            else:
                for a in hapt_atoms:
                    if a in coord_order:
                        coord_dir[a] = hapt_atoms

        elif batom.GetIdx() == tmc_idx:
            coord_dir[eatom.GetIdx()+1] = [eatom.GetIdx()+1]
        else:
            coord_dir[batom.GetIdx()+1] = [batom.GetIdx()+1]

    return coord_dir


def get_possible_coord_permutations(smiles, geometry=None):
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    m = reorder_atoms_to_map(m)

    mhap = rdmolops.DativeBondsToHaptic(m)

    tmc_idx = None
    for a in mhap.GetAtoms():
        if a.GetAtomicNum() in TRANSITION_METALS_NUM:
            tmc_idx = a.GetIdx()


    N_binding_sites = len(mhap.GetAtoms()[tmc_idx].GetBonds())
    #print(N_binding_sites)
    coord_dir = get_hapt_coord_dir(mhap, None)
    coord_orders = []
    initial_coord_order = list(coord_dir.keys())
    if not geometry:
        geometry = geom_dict[N_binding_sites]
    for perm in coord_permutations[geometry]:
        coord_perm = [initial_coord_order[i] for i in perm]
        coord_orders.append(coord_perm)
    return coord_orders



def embed_tmc_complex(smiles, xtb_path, coord_order=None, geometry=None):
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    m = reorder_atoms_to_map(m)
    mhap = rdmolops.DativeBondsToHaptic(m)
    coord_dir = get_hapt_coord_dir(mhap, coord_order)
    if not coord_order: #choose random coord_order
        coord_order = list(coord_dir.keys())
        #print(coord_order)
    charge = Chem.GetFormalCharge(m)
    constraints, angle_constraints, new_mol = setup_constrain_dict_pos(mhap, coord_order, coord_dir)
    for i,atom in enumerate(new_mol.GetAtoms()):
        atom.SetAtomMapNum(i+1)
    sep_embed = embed_smiles_far(Chem.MolToSmiles(new_mol))
    Chem.MolToXYZFile(sep_embed, filename="sep_embed.xyz")
    os.mkdir("ff_embed")
    shutil.move("sep_embed.xyz", "ff_embed") # Move instead of copy
    os.chdir("ff_embed")
    write_constrain_file(constraints, angle_constraints, 0.0025)
    shutil.copy("constrain.inp", "../")


    try:
        output = run_cmd(f"{xtb_path} sep_embed.xyz --opt crude --input constrain.inp --gfnff --parallel 8 --chrg {charge}".format())
        #print(output)
        shutil.copy("xtbopt.xyz", "../sep_embed_ff.xyz")
        os.chdir("../")
        shutil.rmtree("ff_embed")
    except:
        print("stopped at ff embed part")
        os.chdir("../")
        shutil.rmtree("ff_embed")
        return -1, np.nan
    os.mkdir("gfn2_constrain")
    shutil.move("sep_embed_ff.xyz", "gfn2_constrain")
    os.chdir("gfn2_constrain")
    write_constrain_file(constraints, angle_constraints, 0.2)
    shutil.copy("constrain.inp", "../constrain2.inp")
    try:
        output = run_cmd(f"{xtb_path} sep_embed_ff.xyz --opt crude --input constrain.inp --parallel 8 --gfnff --chrg {charge}".format())
        #print(output)
        shutil.copy("xtbopt.xyz", "../sep_embed_final.xyz")
        os.chdir("../")
        shutil.rmtree("gfn2_constrain")
        mol = assign_coordinates2mol(sep_embed, "sep_embed_final.xyz", constraints)
        output = run_cmd(f"{xtb_path} sep_embed_final.xyz --opt verytight --iterations 1000 --input constrain2.inp --gfn2 --parallel 8 --chrg {charge}".format())

        with open("out.out", 'w') as _file:
            _file.writelines(output)
            #lines = _file.readlines()
        with open("out.out", 'r') as _file:
            lines = _file.readlines()
        for line in lines:
            #print(line)
            if "TOTAL ENERGY" in line:
                e = float(line.split()[3])
                #print(line)
        mol = assign_coordinates2mol(sep_embed, "xtbopt.xyz", constraints)
    except:
        print("failed at gfn2 constrain part")
        os.chdir("../")
        shutil.rmtree("gfn2_constrain")
        return -1, np.nan

    return mol, e


def get_tmc_mol(atom_mapped_smiles, xtb_path, coord_order=None, geometry=None, N_tries=30):
    """
    Do N_tries embeddings with final gfn2 optimizations and keep lowest energy
    """
    pos_mols = []
    pos_es = []
    for _ in range(N_tries):
        mol, e = embed_tmc_complex(atom_mapped_smiles, xtb_path, coord_order=coord_order, geometry=geometry)
        pos_mols.append(mol)
        pos_es.append(e)
       # print(f"energy = {e}")

    if np.all(pos_mols) == -1:
        print("embedding failed")
    else:
        #print(coord_order)
        print("")

    return pos_mols, pos_es  # NOTE: returning all the conformations, so I can run r2scan SP on each

if __name__ == "__main__":
    xtb_path = "/home/maria/bin/xtb-6.7.1/xtb-dist/bin/xtb"
    s = "[Rh:1](<-[C-:2]#[O+:6])(<-[C-:3]#[O+:7])(<-[C-:4]#[O+:8])[H:5]"
    possible_coord_orders = get_possible_coord_permutations(s)
    stereo_confs = []
    for coord_order in possible_coord_orders:
        mol = get_tmc_mol(s, xtb_path, coord_order=coord_order)
        stereo_confs.append(mol)

