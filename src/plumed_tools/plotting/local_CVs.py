"""A plotting module, useful for analyzing PLUMED colvar files from local collective variables (CVs).

Contains the parent class 'PLUMEDAnalyzer' which initializes the appropriate
attributes common to all forms of local CV plotting. The @classmethods of
PLUMEDAnalyzer set up consistent settings for generating various kinds of plots
corresponding to CV type. E.g., for the distance CV type, a conversion_factor
of 10 is needed to convert the native PLUMED distance unit of nanometers to the
more desirable unit of Ångstroms, which is useful whether the generated plot is
a KDE, a scatterplot, etc. The various child classes of PLUMEDAnalyzer define
the particular plotting schemes (such as KDE, heatmap). 

Important Note: This module depends on the PLUMED colvar files being generated
in a particular format. Firstly, this module only applies to the plotting of
LOCAL CVs. The resulting colvar files must be dumped in the following format:
{prefix}_{species}.dat. With the file pertaining to only 1 species of local CV,
column 0 representing the timestep, or frame information, and columns 1 to n
representing the local CV data for molecules 1 to n in a system of n molecules.
A dictionary of {prefix} names defined in this module for all relevant CV types
is given below.

PLUMED_prefixes_definitions: dict[CV type (str), prefix (str)] = {
    "DISTANCE": "distances",
    "TORSION": "torsions"
}

Typical usage example: Creating a KDE and Heatmap plot for DISTANCE CV data
from species d1..d8.

See "Using plumed_tools" on the home page of the documentation website
(https://logan-estridge.github.io/plumed-tools/) to learn how to install
plumed_tools as an editable pip package. 

    Directory structure:
    ./location_of_your_MD_input_files/traj/distances_d[1..8].dat # colvar files
    ./my_plotting_script.py

    In 'my_plotting_script.py':
        import plumed_tools
        # or if you like to be more specific: 
        # from plumed_tools.plotting.local_CVs import KDEPlotter, HeatmapPlotter  

        def main():
            root = "location_of_your_MD_input_files"

            dist_species = [f"d{x}" for x in range(1, 9)]
            dist_jobs = [(KDEPlotter.distance(s, directory=root),
                        HeatmapPlotter.distance(s, directory=root)) for s in dist_species]

            for kde, heat in dist_jobs:
                kde.run()
                heat.run()

        if __name__ == "__main__":
            main()
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

class PLUMEDAnalyzer:
    """Base class for handling PLUMED CV data loading and metadata.

    Attributes:
        species (str): The specific identifier for the CV (e.g., 'd8').
        prefix (str): Name of the PLUMED colvar file.
        conversion_factor (float): Multiplier to adjust units (e.g., 10.0 for nm to Å).
        title (str): The main title string for the generated plot.
        x_label (str): The label for the horizontal axis.
        filename (str): The derived path to the input data file.
        directory (str): The directory containing the colvar file to be analyzed. 
    """

    def __init__(self, species: str, prefix: str, conversion_factor: float, title: str, x_label: str, directory: str = ".") -> None:
        """Initializes the instance's attributes for a given subclass with instances corresponding to a particular CV's metadata.

        Args: 
            species: The species of CV, e.g. d8 for distance 8.
            prefix: Base name of the colvar.dat file for the particular CV. Used to locate the colvar file.
            conversion_factor: Scaling factor applied to the raw data values.
            title: Title of the plot corresponding to the particular CV.
            x_label: X axis label of the plot corresponding to the particular CV (including units).
        """
        self.species = species
        self.prefix = prefix
        self.conversion_factor = conversion_factor
        self.title = title
        self.x_label = x_label
        self.directory = directory
        self.filename = os.path.join(self.directory, "traj", f"{prefix}_{species}.dat") # derivative of {prefix} and {species} thus not passed into __init__
        self.output_name = "" # output name of the saved PNG. Needs to be different 
        self.output_path = ""

    def load_data(self) -> np.ndarray:
        """Standardized data loading for all PLUMED plots."""
        if not os.path.exists(self.filename):
            print(f"Warning: {self.filename} not found.")
            return None

        df = pd.read_csv(self.filename, sep='\s+', comment='#', header=None)
        return df.iloc[:, 1:].values * self.conversion_factor

    def _run_safely(self, task_function) -> None:
        """Generic error handler wrapper for child class operations."""
        try:
            task_function()
        except pd.errors.EmptyDataError:
            print(f"Error: {self.filename} contains no data.")
        except pd.errors.ParserError:
            print(f"Error: {self.filename} is poorly formatted.")
        except ValueError as ve:
            print(f"Data Error on {self.species}: {ve}")
        except PermissionError:
            print(f"Permission Error: Cannot access {self.filename}.")
        except Exception as e:
            print(f"Unexpected Error on {self.species}: {e}")
        finally:
            # Ensure plots are closed even if an error occurs to save memory
            if 'plt' in globals() or 'plt' in locals():
                import matplotlib.pyplot as plt
                plt.close()

    @classmethod
    def distance(cls, species: str, directory: str = ".") -> "PLUMEDAnalyzer":
        """Factory for PLUMED distance CV analyses.

        Args:
            species: CV species instantiated by PLUMEDAnalyzer
            directory: The root directory where the colvar file is located. For example,
                        alpha/traj/distances_d8.dat <- the root directory is 'alpha'.
        """
        return cls(
            species=species,
            prefix="distances",
            conversion_factor=10,
            title="Distance",
            x_label="Distance (Å)",
            directory=directory
        )

    @classmethod
    def torsion(cls, species: str, directory: str = ".") -> "PLUMEDAnalyzer":
        """Factory for PLUMED torsion CV analyses.

        Args:
            species: CV species instantiated by PLUMEDAnalyzer
            directory: The root directory where the colvar file is located. For example,
                        alpha/traj/torsions_t1.dat <- the root directory is 'alpha'.
        """
        return cls(
            species=species,
            prefix="torsions",
            conversion_factor=1,
            title="Torsion Angle",
            x_label="Angle (Radians)",
            directory=directory
        )

class KDEPlotter(PLUMEDAnalyzer):
    """Generates KDE plots to see the distribution of the local CV over the system's molecules.

    The KDE for a local CV is plotted for each molecule in the system (the molecules that were
    analyzed with a particular CV and are present in the colvar file). If there are more than 20
    molecules, the legend is omitted. KDE plots are saved in the format {prefix}_kde_{species}.png.
    E.g., for a distance CV, "d8", the plot will be saved as distances_kde_d8.png. The y-axis represents
    the KDE, the x-axis represents the range of the CV data, for example bond distance, etc.
    """
    def run(self) -> None:
        def kde_plot_logic():
            data = self.load_data()
            if data is None: return

            num_frames, num_cols = data.shape

            df_plot = pd.DataFrame(data, columns=[f'Mol {i+1}' for i in range(num_cols)])

            plt.figure(figsize=(10, 6))
            colors = sns.color_palette("colorblind", num_cols)
            for i in range(num_cols):
                sns.kdeplot(
                    data[:, i], 
                    color=colors[i], 
                    label=f'Mol {i+1}',
                    alpha=0.9
                )
            plt.title(f"{self.title} ({self.species}): {num_cols} mols")
            plt.xlabel(self.x_label)
            plt.ylabel("Density")
            
            if num_cols <= 20:
                plt.legend()

            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()

            self.output_name = f"{self.prefix}_kde_{self.species}.png"
            self.output_path = os.path.join(self.directory, "traj", self.output_name)
            plt.savefig(self.output_path, dpi=300)
            plt.close()

        self._run_safely(kde_plot_logic)

class HeatmapPlotter(PLUMEDAnalyzer):
    """Generates Heatmap plots for analyzing CV change over simulation time.

    The threshold values for the colorbar are determined by the np.percentile() function
    to cut off the top and bottom 2 % of values from the range. Heatmap plots are saved in the
    format {prefix}_heatmap_{species}.png. E.g., for a distance CV, "d8", the plot will be saved
    as distances_heatmap_d8.png. The y-axis represents each molecule's CV data present in the 
    colvar file, the x-axis represents the frame index from the simulation trajectory. 
    """

    def run(self, low_thresh: float = None, high_thresh: float = None) -> None:
        def heat_plot_logic():
            data = self.load_data()
            if data is None: return

            _low = low_thresh if low_thresh is not None else np.percentile(data, 2)
            _high = high_thresh if high_thresh is not None else np.percentile(data, 98)

            num_frames, num_cols = data.shape
            plt.figure(figsize=(14,10))

            im = plt.imshow(data.T, aspect='auto', interpolation='nearest',
                            origin='lower', cmap='RdYlBu_r', vmin=_low,
                            vmax=_high)

            cbar = plt.colorbar(im, pad=0.02)
            cbar.set_label(self.x_label, rotation=270, labelpad=15)

            plt.title(f"{self.title} Heatmap {self.species}: {num_cols} Mols")
            plt.xlabel("Frame Index")
            plt.ylabel("Molecule Index")
            plt.yticks(np.arange(8, num_cols + 1, 8))
            
            plt.tight_layout()
            self.output_name = f"{self.prefix}_heatmap_{self.species}.png"
            self.output_path = os.path.join(self.directory, "traj", self.output_name)
            plt.savefig(self.output_path, dpi=300)
            plt.close()
        
        self._run_safely(heat_plot_logic)
