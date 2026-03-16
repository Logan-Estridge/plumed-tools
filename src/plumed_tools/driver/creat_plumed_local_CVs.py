"""Automates creation of PLUMED input files for local CVs.

Provides a class-based interface to generate collective variable (CV)
definitions for large systems of identical molecules.

Typical usage example:
    from plumed_tools.driver.creat_plumed_local_CVs import InputMaker

    paba = InputMaker.pABA(no_of_mols=100)
    water = InputMaker(no_of_mols=100, no_of_atoms=3)

    # For pABA
    paba.write_distances("d3", 3, 2)
    paba.write_distances("d4", 1, 2)
    paba.write_distances("d8", 1, 16)

    # For Water
    water.write_distances("d1", 1, 1)
    water.write_distances("d2", 1, 2)
"""

from typing import List

class InputMaker:
    """Blueprint for generic local CV input scripts.

    Attributes:
        no_of_mols (int): Total number of molecules.
        no_of_atoms (int): Number of atoms per molecule.
    """

    def __init__(self, no_of_mols: int, no_of_atoms: int) -> None:
        self.no_of_mols = no_of_mols
        self.no_of_atoms = no_of_atoms

    @classmethod
    def pABA(cls, no_of_mols: int) -> "InputMaker":
        """Factory for para-aminobenzoic acid (pABA) systems (17 atoms/mol)."""
        return cls(no_of_mols, no_of_atoms=17)

    def _get_atom_indices(self, atoms: List[int], mol_idx: int) -> List[str]:
        offset = (mol_idx - 1) * self.no_of_atoms
        return [str(a + offset) for a in atoms]

    def _validate_indices(self, start_index: int, offset: int) -> None:
        if start_index + offset > self.no_of_atoms:
            raise ValueError(
                f"Incorrect start_index ({start_index}) or offset ({offset}). "
                f"start_index + offset must be <= no_of_atoms ({self.no_of_atoms})."
            )

    def write_angles(self, species: str, atoms: List[int]) -> None:
        """Write angle CVs to 'plumed_angles.dat'."""
        filename = "plumed_angles.dat"
        labels = [f"{species}_{i}" for i in range(1, self.no_of_mols + 1)]
        with open(filename, "a") as f:
            for i, label in enumerate(labels, 1):
                atom_list_str = ",".join(self._get_atom_indices(atoms, i))
                f.write(f"{label}: ANGLE ATOMS={atom_list_str}\n")
            f.write(f"\nPRINT ARG={','.join(labels)} FILE=angles_{species}.dat\n\n")

    def write_distances(self, species: str, start_index: int, offset: int) -> None:
        """Write distance CVs to 'plumed_distances.dat'."""
        self._validate_indices(start_index, offset)
        filename = "plumed_distances.dat"
        labels = [f"{species}_{i}" for i in range(1, self.no_of_mols + 1)]
        with open(filename, "a") as f:
            for i, label in enumerate(labels, 1):
                mol_offset = (i - 1) * self.no_of_atoms
                m1 = start_index + mol_offset
                m2 = start_index + offset + mol_offset
                f.write(f"{label}: DISTANCE ATOMS={m1},{m2}\n")
            f.write(f"\nPRINT ARG={','.join(labels)} FILE=distances_{species}.dat\n\n")

    def write_gyrations_list(self, species: str, atoms: List[int]) -> None:
        """Write Radius of Gyration CVs to 'plumed_radii_of_gyration.dat'."""
        filename = "plumed_radii_of_gyration.dat"
        labels = [f"{species}_{i}" for i in range(1, self.no_of_mols + 1)]
        with open(filename, "a") as f:
            for i, label in enumerate(labels, 1):
                atom_list_str = ",".join(self._get_atom_indices(atoms, i))
                f.write(f"{label}: GYRATION TYPE=RADIUS ATOMS={atom_list_str}\n")
            f.write(f"\nPRINT ARG={','.join(labels)} FILE=radii_of_gyration_{species}.dat\n\n")

    def write_gyrations_range(
        self, species: str, start_index: int, offset: int
    ) -> None:
        """Write Radius of Gyration CVs to 'plumed_radii_of_gyration.dat'."""
        self._validate_indices(start_index, offset)
        filename = "plumed_radii_of_gyration.dat"
        labels = [f"{species}_{i}" for i in range(1, self.no_of_mols + 1)]
        with open(filename, "a") as f:
            for i, label in enumerate(labels, 1):
                mol_offset = (i - 1) * self.no_of_atoms
                m1 = start_index + mol_offset
                m2 = start_index + offset + mol_offset
                f.write(f"{label}: GYRATION TYPE=RADIUS ATOMS={m1}-{m2}\n")
            f.write(f"\nPRINT ARG={','.join(labels)} FILE=radii_of_gyration_{species}.dat\n\n")

    def write_torsions(self, species: str, a1: int, a2: int, a3: int, a4: int) -> None:
        """Write torsion CVs to 'plumed_torsions.dat'."""
        filename = "plumed_torsions.dat"
        labels = [f"{species}_{i}" for i in range(1, self.no_of_mols + 1)]
        with open(filename, "a") as f:
            for i, label in enumerate(labels, 1):
                off = (i - 1) * self.no_of_atoms
                f.write(f"{label}: TORSION ATOMS={a1 + off},{a2 + off},{a3 + off},{a4 + off}\n")
            f.write(f"\nPRINT ARG={','.join(labels)} FILE=torsions_{species}.dat\n\n")
