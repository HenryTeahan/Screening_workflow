import subprocess
import itertools
import py3Dmol

import numpy as np

from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Geometry import Point3D
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers


from .xyz2mol import xyz2mol_local, xyz2mol_TMC
#import xyz2mol.xyz2mol_local
#import xyz2mol.xyz2mol_TMC

def draw3D(ms, p=None, confId=-1, removeHs=True,colors=('cyanCarbon','redCarbon','blueCarbon'), with_index=False):
        if p is None:
            p = py3Dmol.view(width=400, height=400)
            
        p.removeAllModels()
        if removeHs:
            m = Chem.RemoveHs(ms)
        else:
            m = deepcopy(ms)
        IPythonConsole.addMolToView(m,p,confId=confId)

        p.setStyle({'stick':{}})#{'colorscheme':colors[i%len(colors)]}})
        #p.addPropertyLabels("index","","")
        if with_index:
            coords = m.GetConformer().GetPositions()

            for idx, atom in enumerate(m.GetAtoms()):
                label = f"{idx+1}"  # Atomic index
                position = coords[idx]  # 3D coordinates of the atom
                p.addLabel(label,{"position":{"x":position[0],"y":position[1],"z":position[2]},"fontSize": 14})
        #    #p.addLabel(label, {"position": position.tolist(), "fontSize": 14, "color": "black"})
        p.zoomTo()
        return p.show()


def run_cmd(cmd):
    """
    Run command line
    """
    cmd = cmd.split()
    #print(cmd)
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    output, err = p.communicate()
    return output.decode('utf-8')


def reorder_atoms_to_map(mol):

    """
    Reorders the atoms in a mol objective to match that of the mapping
    """

    atom_map_order = np.zeros(mol.GetNumAtoms()).astype(int)
    for atom in mol.GetAtoms():
        map_number = atom.GetAtomMapNum()-1
        atom_map_order[map_number] = atom.GetIdx()
    mol = Chem.RenumberAtoms(mol, atom_map_order.tolist())
    return mol

def bonds_getting_formed_or_broken(rsmi, psmi):
    """
    Based on the reaction and product structure, the bonds that are
    fomed/broken are singled out for contraintment
    the difference in the afjacency matric tells whether bond has been formed
    (+1) or bond is broken (-1)
    """

    bond_pairs_changed = []
    rmol = Chem.MolFromSmiles(rsmi, sanitize=False)
    pmol = Chem.MolFromSmiles(psmi, sanitize=False)
    n_atoms = rmol.GetNumAtoms()
    rmol = reorder_atoms_to_map(rmol)
    pmol = reorder_atoms_to_map(pmol)
    Chem.SanitizeMol(rmol)
    Chem.SanitizeMol(pmol)

    p_ac = rdmolops.GetAdjacencyMatrix(pmol)
    r_ac = rdmolops.GetAdjacencyMatrix(rmol)

    difference_mat = p_ac - r_ac
    for combination in itertools.combinations(range(n_atoms), 2):
        combination = list(combination)
        bond_change = difference_mat[combination[0], combination[1]]
        if bond_change != 0:
            bond_pairs_changed.append(combination)

    return bond_pairs_changed

def get_core_atoms(active_atoms, reactant, product):
    core_atoms = []
    for combination in itertools.combinations(active_atoms, 2):
        ai, aj = list(combination)
        path_reactant = list(rdmolops.GetShortestPath(reactant, ai, aj))
        core_atoms += path_reactant
        path_product = list(rdmolops.GetShortestPath(product, ai, aj))
        core_atoms += path_reactant
    core_atoms = list(set(core_atoms))

    for ai in active_atoms:
        neighbors = [x.GetIdx() for x in reactant.GetAtomWithIdx(ai).GetNeighbors()]
        core_atoms += neighbors
        neighbors = [x.GetIdx() for x in product.GetAtomWithIdx(ai).GetNeighbors()]
        core_atoms += neighbors

    core_atoms = list(set(core_atoms))

    return core_atoms

def find_removable_atoms(core_atoms, reactant, product):
    start_atoms = core_atoms.copy()
    total_remove_atoms = []
    atoms_to_pad = []
    while start_atoms:
        for aidx in start_atoms:
            ai = reactant.GetAtomWithIdx(aidx)
            neighbors = [x.GetIdx() for x in ai.GetNeighbors()]
            #print(neighbors)
            #starting_atoms = []
            for x in neighbors:
                if x not in core_atoms:
                    ax = reactant.GetAtomWithIdx(x)
                    if ax.GetAtomicNum()==1:
                        core_atoms.append(x)
                    else:
                        rbond = reactant.GetBondBetweenAtoms(aidx, x)
                        pbond = product.GetBondBetweenAtoms(aidx, x)
                        if rbond.IsInRing() or pbond.IsInRing():
                            core_atoms.append(x)
                            start_atoms.append(x)
                            #print("atom in ring")
                        elif rbond.GetBondType() != Chem.BondType.SINGLE and pbond.GetBondType() != Chem.BondType.SINGLE:
                            core_atoms.append(x)
                            start_atoms.append(x)
                        else:
                            atoms_to_pad.append(aidx)
                            #remove_atoms = []
                            remove_neighbors = [aa.GetIdx() for aa in ax.GetNeighbors()]
                            #print(remove_neighbors)
                            remove_neighbors.remove(aidx)
                            #print(remove_neighbors)
                            remove_atoms = [x] + remove_neighbors.copy()
                            while remove_neighbors:
                                for j in remove_neighbors:
                                    #print(j)
                                    aj = reactant.GetAtomWithIdx(j)
                                    temp_neighbors = [y.GetIdx() for y in aj.GetNeighbors()]
                                    for k in temp_neighbors:
                                        if k not in remove_atoms and k not in core_atoms:
                                            remove_atoms.append(k)
                                            remove_neighbors.append(k)
                                    remove_neighbors.remove(j)
                            total_remove_atoms += remove_atoms
                            print("remove atoms:", remove_atoms)


                            print("should cut bond between", x, aidx)

                            print(rbond.GetBondType())
            start_atoms.remove(aidx)

            total_remove_atoms = list(set(total_remove_atoms))
    return total_remove_atoms, atoms_to_pad

def chiral_tags(in_mol):
    """
    Tag methylene and methyl groups with a chiral tag priority defined
    from the atom index of the hydrogens
    """
    mol = deepcopy(in_mol)
    li_list = []
    smarts_ch2 = '[!#1][*]([#1])([#1])([!#1])'
    atom_sets = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts_ch2))
    for atoms in atom_sets:
        atoms = sorted(atoms[2:4])
        prioritized_H = atoms[-1]
        li_list.append(prioritized_H)
        mol.GetAtoms()[prioritized_H].SetAtomicNum(9)
    smarts_ch3 = '[!#1][*]([#1])([#1])([#1])'
    atom_sets = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts_ch3))
    for atoms in atom_sets:
        atoms = sorted(atoms[2:])
        H1 = atoms[-1]
        H2 = atoms[-2]
        li_list.append(H1)
        li_list.append(H2)
        mol.GetAtoms()[H1].SetAtomicNum(9)
        mol.GetAtoms()[H2].SetAtomicNum(17)

    #Chem.AssignAtomChiralTagsFromStructure(mol, -1)
    #rdmolops.AssignStereochemistry(mol)
    rdmolops.AssignStereochemistry(mol, cleanIt=True, flagPossibleStereoCenters=True, force=True)

    mol = next(EnumerateStereoisomers(mol))
    rdmolops.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)

    for atom_idx in li_list:
        mol.GetAtoms()[atom_idx].SetAtomicNum(1)

    return mol

def choose_resonance_structure(mol):
    """
    This function creates all resonance structures of the mol object, counts
    the number of rotatable bonds for each structure and chooses the one with
    fewest rotatable bonds (most 'locked' structure)
    """
    #resonance_mols = rdchem.ResonanceMolSupplier(mol,
    #                                             rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION)
    resonance_mols = rdchem.ResonanceMolSupplier(mol)
    res_status = True
    new_mol = None
    if not resonance_mols:
        print("using input mol")
        new_mol = mol
        res_status = False
    for res_mol in resonance_mols:
        Chem.SanitizeMol(res_mol)
        n_rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(res_mol)
        if new_mol is None:
            smallest_rot_bonds = n_rot_bonds
            new_mol = res_mol
        if n_rot_bonds < smallest_rot_bonds:
            smallest_rot_bonds = n_rot_bonds
            new_mol = res_mol

    Chem.DetectBondStereochemistry(new_mol, -1)
    rdmolops.AssignStereochemistry(new_mol, flagPossibleStereoCenters=True,
                                   force=True)
    Chem.AssignAtomChiralTagsFromStructure(new_mol, -1)
    return new_mol, res_status

def extract_smiles(xyz_file, charge, allow_charge=True, asMol=False):
    """
    uses xyz2mol to extract smiles with as much 3d structural information as
    possible
    """
    atoms, _, xyz_coordinates = xyz2mol_local.read_xyz_file(xyz_file)
    try:
        input_mol = xyz2mol_local.xyz2mol(atoms, xyz_coordinates, charge=charge,
                                          use_graph=True,
                                          allow_charged_fragments=allow_charge,
                                          use_obabel=True, use_atom_maps=True,
                                          embed_chiral=True)
    except:
        input_mol = xyz2mol_local.xyz2mol(atoms, xyz_coordinates, charge=charge,
                                          use_graph=True,
                                          allow_charged_fragments=allow_charge,
                                          use_huckel=True, use_atom_maps=True,
                                          embed_chiral=True)

    input_mol = reorder_atoms_to_map(input_mol)
    structure_mol, res_status = choose_resonance_structure(input_mol)
    structure_mol = chiral_tags(structure_mol)
    rdmolops.AssignStereochemistry(structure_mol)
    if asMol:
        return structure_mol
    structure_smiles = Chem.MolToSmiles(structure_mol)

    return structure_smiles

def canonicalize_smiles(structure_smiles):
    """
    remove all structural info an atom mapping information
    """
    mol = Chem.MolFromSmiles(structure_smiles, sanitize=False)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    canonical_smiles = Chem.MolToSmiles(mol)

    return Chem.MolToSmiles(Chem.MolFromSmiles(canonical_smiles))

def get_indices(output, pattern):
    indices = []
    for i, l in enumerate(output):
        if pattern in l:
            indices.append(i)
    return indices

def extract_scf_energy_xtb(output, asFile=False):
    if asFile:
        with open(output, 'r') as _file:
            lines = _file.readlines()
    else:
        lines=output
    indices = get_indices(lines, "TOTAL ENERGY")
    line = lines[indices[-1]]
    E_hartree = float(line.split()[3])

    return E_hartree

def embed_smiles_far(smiles):#, RANDOM_SEED=1):
    """
    create 3D conformer with atom order matching atom mapping. If more than one
    fragment present, they are moved 0.5*n_atoms from each other
    """
    embedded_fragments = []
    mol_sizes = []

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = reorder_atoms_to_map(mol)

    fragments = Chem.GetMolFrags(mol, asMols=True)
    n_atoms = mol.GetNumAtoms()
    n_atom_fragments = [fragment.GetNumAtoms() for fragment in fragments]
    coordinates = np.zeros((n_atoms, 3))
    fragments = [x for _, x in sorted(zip(n_atom_fragments, list(fragments)),
                 key=lambda pair: pair[0], reverse=True)]

    for fragment in fragments:
        Chem.SanitizeMol(fragment)
        rdmolops.AssignStereochemistry(fragment)
        status = AllChem.EmbedMolecule(fragment, maxAttempts=1000)#, randomSeed=RANDOM_SEED)
        if status == -1:
            print('fragment could not be embedded')
            return -1

        AllChem.MMFFOptimizeMolecule(fragment)
        embedded_fragments.append(fragment)
        dm = AllChem.Get3DDistanceMatrix(fragment)
        biggest_distance = np.amax(dm)
        mol_sizes.append(biggest_distance)

    for n_frag, fragment in enumerate(embedded_fragments):
        if n_frag == 0:
            random_vector = np.zeros(3)
        else:
            random_vector = np.random.rand(3)*2-1
            random_vector = random_vector / np.linalg.norm(random_vector)

        conformer = fragment.GetConformer()
        translation_distance = 0.5*(mol_sizes[n_frag]+mol_sizes[0])+5

        for i, atom in enumerate(fragment.GetAtoms()):
            atom_id = atom.GetAtomMapNum()-1
            coord =  conformer.GetAtomPosition(i)
            coord += Point3D(*(translation_distance*random_vector))
            coordinates[atom_id, :] = coord

    Chem.SanitizeMol(mol)
    status = AllChem.EmbedMolecule(mol, maxAttempts=1000)
    if status == -1:
        return -1
    conf = mol.GetConformer()
    for i in range(n_atoms):
        x, y, z = coordinates[i, :]
        conf.SetAtomPosition(i, Point3D(x, y, z))

    return mol

def extract_TMC_smiles(
    xyz_file,
    charge,
    allow_charge=True,
    use_atom_maps=True
):
    """
    uses xyz2mol to extract smiles with as much 3d structural information as
    possible
    """
    # atoms, _, xyz_coordinates = xyz2mol_local.read_xyz_file(xyz_file)

    input_mol = xyz2mol_TMC.get_tmc_mol(xyz_file[:-4], charge, use_atom_maps=use_atom_maps)

    # input_mol = xyz2mol_local.xyz2mol(atoms, xyz_coordinates, charge=charge,
    #                                      use_graph=True,
    #                                      allow_charged_fragments=allow_charge,
    #                                      use_huckel=True, use_atom_maps=True,
    #                                      embed_chiral=True)
    # except:
    #    input_mol = xyz2mol_local.xyz2mol(atoms, xyz_coordinates, charge=charge,
    #                                      use_graph=True,
    #                                      allow_charged_fragments=allow_charge,
    #                                      use_huckel=False, use_atom_maps=True,
    #                                      embed_chiral=True)

    # print(f"Reordering and choosing resonance structure of {Chem.MolToSmiles(input_mol)}")
    # print(Chem.MolToSmiles(input_mol))
    input_mol = reorder_atoms_to_map(input_mol)
    structure_mol, res_status = choose_resonance_structure(input_mol)
    structure_mol = chiral_tags(structure_mol)
    rdmolops.AssignStereochemistry(structure_mol)
    structure_smiles = Chem.MolToSmiles(structure_mol)
    #structure_smiles = simplify_metal_center(structure_smiles, metal="Pd")

    return structure_smiles

def write_xyz_file(mol, file_name):
    """
    write xyz file with cooridnates optimized w ff based on input smiles
    """
    n_atoms = mol.GetNumAtoms()
    charge = Chem.GetFormalCharge(mol)

    symbols = [a.GetSymbol() for a in mol.GetAtoms()]

    with open(file_name, 'w') as _file:
        _file.write(str(n_atoms)+'\n\n')
        for atom, symbol in enumerate(symbols):
            coord = mol.GetConformers()[0].GetAtomPosition(atom)
            line = " ".join((symbol, str(coord.x), str(coord.y), str(coord.z),
                             "\n"))
            _file.write(line)

    return mol


