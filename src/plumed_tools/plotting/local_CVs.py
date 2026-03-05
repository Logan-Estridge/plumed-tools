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

PLUMED_prefixes_definitions: dict[str, str] = {
    "DISTANCE": "distances",
    "TORSION": "torsions",
    "ANGLE": "angles",
    "GYRATION": "radii_of_gyration"
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
        import plumed_tools as pt
        # or if you like to be more specific: 
        # from pt.plotting.local_CVs import KDEPlotter, HeatmapPlotter  

        def main():
            root = "location_of_your_MD_input_files"

            dist_species = [f"d{x}" for x in range(1, 9)]
            dist_jobs = [(pt.KDEPlotter.distance(s, directory=root),
                        pt.HeatmapPlotter.distance(s, directory=root)) for s in dist_species]

            for kde, heat in dist_jobs:
                kde.run()
                heat.run()

        if __name__ == "__main__":
            main()
"""
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
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

    def load_data(self) -> pd.DataFrame:
        """Standardized data loading for all PLUMED plots."""
        if not os.path.exists(self.filename):
            raise FileNotFoundError(f"File not found: {self.filename}")
            return None

        df = pd.read_csv(self.filename, sep='\s+', comment='#', header=None, engine='c', memory_map=True)
        
        if self.conversion_factor != 1.0:
            df.iloc[:, 1:] *= self.conversion_factor
        return df

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
            if plt.fignum_exists(1): plt.figure(1).clf()
            if plt.fignum_exists(2): plt.figure(2).clf()

    @classmethod
    def angle(cls, species: str, directory: str = ".") -> "PLUMEDAnalyzer":
        """Factory for PLUMED ANGLE CV analyses.

        Args:
            species: CV species instantiated by PLUMEDAnalyzer
            directory: The root directory where the colvar file is located. For example,
                        alpha/traj/torsions_t1.dat <- the root directory is 'alpha'.
        """
        return cls(
            species=species,
            prefix="angle",
            conversion_factor=1,
            title="Angle",
            x_label="Angle (Radians)",
            directory=directory
        )

    @classmethod
    def distance(cls, species: str, directory: str = ".") -> "PLUMEDAnalyzer":
        """Factory for PLUMED DISTANCE CV analyses.

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
    def rad_gyration(cls, species: str, directory: str = ".") -> "PLUMEDAnalyzer":
        """Factory for PLUMED GYRATION CV analyses.

        Args:
            species: CV species instantiated by PLUMEDAnalyzer
            directory: The root directory where the colvar file is located. For example,
                        alpha/traj/torsions_t1.dat <- the root directory is 'alpha'.
        """
        return cls(
            species=species,
            prefix="radii_of_gyration",
            conversion_factor=10,
            title="Distance",
            x_label="Distance (Å)",
            directory=directory
        )

    @classmethod
    def torsion(cls, species: str, directory: str = ".") -> "PLUMEDAnalyzer":
        """Factory for PLUMED TORSION CV analyses.

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
    def run(self, data: pd.DataFrame = None) -> None:
        """Generates KDE plot

        Args:
            data: Defining the data source here allows prevention of duplicate data loading if using 
                  for loop and plotting multiple kinds of plots from the same CSV file.
        """
        def kde_plot_logic():
            df = data if data is not None else self.load_data()
            if df is None: return

            df_plot = df.drop(columns=[0])
            num_mols = df_plot.shape[1]

            fig = plt.figure(1, figsize=(10, 6))
            fig.clf()
            ax = fig.add_subplot(111)

            sns.kdeplot(data=df_plot, ax=ax, palette="colorblind", alpha=0.9, legend=True)

            ax.set_title(f"{self.title} ({self.species}): {num_mols} mols")
            ax.set_xlabel(self.x_label)
            ax.set_ylabel("Density")
            ax.grid(axis='y', alpha=0.3)

            if num_mols > 20:
                legend = ax.gca().get_legend()
                if legend: legend.remove()

            self.output_name = f"{self.prefix}_kde_{self.species}.png"
            self.output_path = os.path.join(self.directory, "traj", self.output_name)
            fig.savefig(self.output_path, dpi=100, bbox_inches='tight')

        self._run_safely(kde_plot_logic)

class HeatmapPlotter(PLUMEDAnalyzer):
    """Generates Heatmap plots for analyzing CV change over simulation time.

    The threshold values for the colorbar are determined by the np.percentile() function
    to cut off the top and bottom 2 % of values from the range. Heatmap plots are saved in the
    format {prefix}_heatmap_{species}.png. E.g., for a distance CV, "d8", the plot will be saved
    as distances_heatmap_d8.png. The y-axis represents each molecule's CV data present in the 
    colvar file, the x-axis represents the frame index from the simulation trajectory. 
    """
    def run(self, data: pd.DataFrame = None, low_thresh: float = None, high_thresh: float = None) -> None:
        """Generates heatmap plot

        Args:
            data: Defining the data source here allows prevention of duplicate data loading if using 
                  for loop and plotting multiple kinds of plots from the same CSV file.
            low_thresh: Given x-range of data, values below this percentile will be excluded from the cbar
            high_thresh: Given x-range of data, values above this percentile will be excluded from the cbar
        """
        def heat_plot_logic():
            df = data if data is not None else self.load_data()
            if df is None: return

            vals = df.iloc[:, 1:].values
            num_mols = vals.shape[1]

            _low = low_thresh if low_thresh is not None else np.percentile(vals, 2)
            _high = high_thresh if high_thresh is not None else np.percentile(vals, 98)

            fig = plt.figure(2, figsize=(14,10))
            fig.clf()
            ax = fig.add_subplot(111)

            im = ax.imshow(vals.T, aspect='auto', interpolation='nearest',
                            origin='lower', cmap='RdYlBu_r', vmin=_low,
                            vmax=_high)

            cbar = fig.colorbar(im, ax=ax, pad=0.02)
            cbar.set_label(self.x_label, rotation=270, labelpad=15)

            ax.set_title(f"{self.title} Heatmap {self.species}: {num_mols} Mols")
            ax.set_xlabel("Frame Index")
            ax.set_ylabel("Molecule Index")

            if num_mols >= 8:
                ax.set_yticks(np.arange(7, num_mols, 8))
                ax.set_yticklabels(np.arange(8, num_mols + 1, 8))
            
            self.output_name = f"{self.prefix}_heatmap_{self.species}.png"
            self.output_path = os.path.join(self.directory, "traj", self.output_name)
            fig.savefig(self.output_path, dpi=100, bbox_inches='tight')
        
        self._run_safely(heat_plot_logic)
