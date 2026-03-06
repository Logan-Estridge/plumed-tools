"""
plumed-tools: Analysis and Visualization for PLUMED
====================================================

A collection of various modules/scripts relevant to PLUMED collective variable
analysis of molecular dynamics trajectories.

Main Components:
- driver/ : Contains modules/scripts relevant to PLUMED driver analysis.
- plotting/ : Contains modules/scripts relevant to creating plots from PLUMED colvar files.
    |- local_CVs.py : A plotting module, useful for analyzing PLUMED colvar files from local collective variables (CVs).
"""

from .driver import InputMaker
from .plotting import PLUMEDAnalyzer, KDEPlotter, HeatmapPlotter

__version__ = "0.1.0"
