"""A module automating the creation of PLUMED input files for local CVs.

This module provides a class-based interface to generate collective variable (CV)
definitions for large systems of identical molecules, such as molecular crystals.

Typical usage example:
```python
from plumed_tools.driver.creat_plumed_local_CVs import InputMaker as IM

paba = IM.pABA(no_of_mols=100)
water = IM(no_of_mols=100, no_of_atoms=3)

def main():
    # For pABA
    paba.write_distances("d3", 3, 2)
    paba.write_distances("d4", 1, 2)
    paba.write_distances("d8", 1, 16)

    # For Water
    water.write_distances("d1", 1, 1)
    water.write_distances("d2", 1, 2)

if __name__ == "__main__":
    main()
```
"""

from typing import List

class InputMaker:
    """Class defining the blueprint for generic local CV input scripts.

    This class encapsulates system-wide parameters to automate the generation
    of PLUMED input files. Filenames are automatically generated based on
    the CV type (e.g., 'plumed_distance.dat').

    Attributes:
        no_of_mols (int): The total number of molecules in the system.
        no_of_atoms (int): The number of atoms per molecule.
    """

    def __init__(self, no_of_mols: int, no_of_atoms: int) -> None:
        """Initializes InputMaker with system dimensions.

        Args:
            no_of_mols: Number of molecules in the simulation box.
            no_of_atoms: Number of atoms in a single molecule.
        """
        self.no_of_mols = no_of_mols
        self.no_of_atoms = no_of_atoms

    @classmethod
    def pABA(cls, no_of_mols: int) -> "InputMaker":
        """Factory method for para-aminobenzoic acid (pABA) systems.

        Sets the atoms per molecule to the standard 17 for pABA.

        Args:
            no_of_mols: Number of pABA molecules in the system.

        Returns:
            An instance of InputMaker configured for pABA (17 atoms/mol).
        """
        return cls(no_of_mols, no_of_atoms=17)

    def write_angles(self, species: str, atoms: List[int]) -> None:    
        labels = []    
        filename = "plumed_angles.dat"
        with open(filename, 'a') as f:
            for i in range(1, self.no_of_mols + 1):    
                current_atoms = [str(a + (i - 1) * self.no_of_atoms) for a in atoms]
                label = f"{species}_{i}"    
                labels.append(label)    
                atom_list_str = ",".join(current_atoms)
                f.write(f"{label}: ANGLE ATOMS={atom_list_str}\n")    
                
            f.write("\n")    
            args = ",".join(labels)    
            f.write(f"PRINT ARG={args} FILE=angle_{species}.dat\n")
            f.write("\n")

    def write_distances(self, species: str, start_index: int, offset: int) -> None:
        """Writes distance CVs to 'plumed_distance.dat'.

        Args:
            species: Identifier for this CV set (e.g., 'd1').
            start_index: The index of the first atom in the first molecule.
            offset: The relative index difference to the second atom.
        """
        if start_index + offset > self.no_of_atoms:
            print(
                f"ERROR: Incorrect definition of start_index ({start_index}) or offset ({offset})."
                f"\nCheck the no_of_atoms ({self.no_of_atoms}) in your molecules and try again."
                "\nThe start_index + offset must be less than or equal to the no_of_atoms."
            )
            return None
        labels = []
        filename = "plumed_distance.dat"
        with open(filename, 'a') as f:
            for i in range(1, self.no_of_mols + 1):
                mol_offset = (i - 1) * self.no_of_atoms
                m1, m2 = start_index + mol_offset, start_index + offset + mol_offset
                label = f"{species}_{i}"
                labels.append(label)
                f.write(f"{label}: DISTANCE ATOMS={m1},{m2}\n")
            args = ",".join(labels)
            f.write(f"\nPRINT ARG={args} FILE=distances_{species}.dat\n\n")

    def write_gyrations_list(self, species: str, atoms: List[int]) -> None:
        """Writes Radius of Gyration CVs to 'plumed_gyration.dat'.

        Args:
            species: Identifier for this CV set (e.g., 'rg1').
            atoms: List of atoms in the molecule over which to calculate the rg.
        """
        labels = []
        filename = "plumed_radii_of_gyration.dat"
        with open(filename, 'a') as f:
            for i in range(1, self.no_of_mols + 1):    
                current_atoms = [str(a + (i - 1) * self.no_of_atoms) for a in atoms]
                label = f"{species}_{i}"    
                labels.append(label)    
                atom_list_str = ",".join(current_atoms)
                f.write(f"{label}: GYRATION TYPE=RADIUS ATOMS={atom_list_str}\n")    
                
            f.write("\n")    
            args = ",".join(labels)    
            f.write(f"PRINT ARG={args} FILE=radii_of_gyration_{species}.dat\n\n")

    def write_gyrations_range(self, species: str, start_index: int, offset: int) -> None:
        """Writes Radius of Gyration CVs to 'plumed_gyration.dat'.

        Args:
            species: Identifier for this CV set (e.g., 'rg1').
            start_index: The first atom index of the molecule range.
            offset: The distance to the last atom index of the range.
        """
        if start_index + offset > self.no_of_atoms:
            print(
                f"ERROR: Incorrect definition of start_index ({start_index}) or offset ({offset})."
                f"\nCheck the no_of_atoms ({self.no_of_atoms}) in your molecules and try again."
                "\nThe start_index + offset must be less than or equal to the no_of_atoms."
            )
            return None
        labels = []
        filename = "plumed_radii_of_gyration.dat"
        with open(filename, 'a') as f:
            for i in range(1, self.no_of_mols + 1):
                mol_offset = (i - 1) * self.no_of_atoms
                m1, m2 = start_index + mol_offset, start_index + offset + mol_offset
                label = f"{species}_{i}"
                labels.append(label)
                f.write(f"{label}: GYRATION TYPE=RADIUS ATOMS={m1}-{m2}\n")
            args = ",".join(labels)
            f.write(f"\nPRINT ARG={args} FILE=radii_of_gyration_{species}.dat\n\n")

    def write_torsions(self, species: str, a1: int, a2: int, a3: int, a4: int) -> None:
        """Writes torsion CVs to 'plumed_torsion.dat'.

        Args:
            species: Identifier for this CV set (e.g., 't1').
            a1: Atom 1 of the torsion angle.
            a2: Atom 2 of the torsion angle.
            a3: Atom 3 of the torsion angle.
            a4: Atom 4 of the torsion angle.
        """
        labels = []
        filename = "plumed_torsion.dat"
        with open(filename, 'a') as f:
            for i in range(1, self.no_of_mols + 1):
                off = (i - 1) * self.no_of_atoms
                label = f"{species}_{i}"
                labels.append(label)
                f.write(f"{label}: TORSION ATOMS={a1+off},{a2+off},{a3+off},{a4+off}\n")
            args = ",".join(labels)
            f.write(f"\nPRINT ARG={args} FILE=torsions_{species}.dat\n\n")

