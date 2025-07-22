"""Protein structures

"""

import numpy as np
import math
from collections import defaultdict
from Bio.PDB import PPBuilder
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqIO.PdbIO import PdbSeqresIterator


def rgyrate(structure, chain_id=None, atom_selector=lambda a: True):
    """Compute the radius of gyration from a Biopython Structure object.

    Parameters:
        structure (Bio.PDB.Structure.Structure): The parsed structure.
        chain_id (str or None): If specified, restrict to a single chain.
        atom_selector (function): Function to filter atoms (default: all atoms).

    Returns:
        float: Radius of gyration in angstroms.
    """
    x = []
    mass = []
    for model in structure:
        for chain in model:
            if chain_id is not None and chain.id != chain_id:
                continue
            for residue in chain:
                for atom in residue:
                    if atom_selector(atom):
                        x.append(atom.coord)
                        mass.append(atom.mass)
    x = np.array(x)
    mass = np.array(mass)

    if len(x) == 0:
        raise ValueError("No atoms selected.")

    xm = [(m*i,m*j,m*k) for (i,j,k),m in zip(x, mass)]
    tmass = sum(mass)
    rr = sum(mi*i+mj*j+mk*k for (i,j,k),(mi,mj,mk) in zip(x, xm))
    mm = sum((sum(i)/tmass)**2 for i in zip(*xm))
    rg = np.sqrt(rr/tmass - mm)
    return rg


def classify_secondary_structure(phi, psi):
    """Classify residue as helix, sheet, or coil based on dihedral angles"""
    if phi is None or psi is None:
        return 'coil'
    
    phi_deg = math.degrees(phi)
    psi_deg = math.degrees(psi)
    
    # Alpha-helix range (-57,-47)
    if (-90 < phi_deg < -30) and (-77 < psi_deg < -17):
        return 'helix'
    # Beta-sheet range (-119,113)
    elif (-180 < phi_deg < -40) and (90 < psi_deg < 180) or (-180 < phi_deg < -40) and (-180 < psi_deg < -100):
        return 'sheet'
    else:
        return 'coil'
    

def get_secondary_structure_proportions(structure):
    """Calculate secondary structure proportions using dihedral angles.
    
    Args:
        structure: Protein structure.
    """
    ppb = PPBuilder()
    
    secondary_structure = defaultdict(int)
    total_residues = 0
    
    for pp in ppb.build_peptides(structure):
        # Get phi/psi angles for each residue
        angles = pp.get_phi_psi_list()
        for i, residue in enumerate(pp):
            total_residues += 1
            phi, psi = angles[i]
            ss_type = classify_secondary_structure(phi, psi)
            secondary_structure[ss_type] += 1
    
    if total_residues > 0:
        return {
            'helix': secondary_structure['helix'] / total_residues,
            'sheet': secondary_structure['sheet'] / total_residues,
            'coil': secondary_structure['coil'] / total_residues
        }
    return {'helix': 0., 'sheet': 0., 'coil': 0.}


def compute_sasa(structure, level="S"):
    """Compute surface accessibility surface area (SASA) using Shrake-Rupley.
    
    Args:
        structure: Protein structure.
        level (str): level at which to compute the SASA

    Returns:
        Average SASA value at the indicated level.
    """
    sr = ShrakeRupley()
    sr.compute(structure, level=level)
    return structure.sasa


def analyze_sequence(structure):
    """Wrapper for ProteinAnalysis to perform analysis given structure.

    Args:
        structure: Protein structure.

    Returns:
        Analyzed sequence object.
    """
    analyzed_seq = ProteinAnalysis(_struct2seq(structure))
    return analyzed_seq


def _struct2seq(structure):
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        protein_seq = pp.get_sequence()  # returns a Bio.Seq.Seq object
        break
    return str(protein_seq)
