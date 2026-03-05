"""A plotting script for local CVs"""
import os
import shutil
import time
from typing import Tuple, List, Dict, Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from plumed_tools.plotting.local_CVs import KDEPlotter, HeatmapPlotter

def process_species(task: Tuple[str, str, str]) -> str:
    """
    Worker function to process a single CV species in a separate process.

    Args:
        task: A tuple containing (species_name, root_directory, summary_directory).

    Returns:
        A status string indicating SUCCESS, SKIP, or ERROR for the task.
    """
    species, root, summary_dir = task

    # Map prefixes to factory methods: (KDE_Factory, Heatmap_Factory)
    factories: Dict[str, Tuple[Callable, Callable]] = {
        'd':  (KDEPlotter.distance,      HeatmapPlotter.distance),
        't':  (KDEPlotter.torsion,       HeatmapPlotter.torsion),
        'a':  (KDEPlotter.angle,         HeatmapPlotter.angle),
        'rg': (KDEPlotter.rad_gyration,  HeatmapPlotter.rad_gyration)
    }

    prefix = 'rg' if species.startswith('rg') else species[0]
    if prefix not in factories:
        return f"ERROR: Unknown species prefix for {species}"

    kde_factory, heat_factory = factories[prefix]

    # Initialize plotters
    kde = kde_factory(species, directory=root)
    heat = heat_factory(species, directory=root)

    # Validate file existence before loading
    if not os.path.exists(kde.filename):
        return f"SKIP: {kde.filename} not found."

    # Load data once to be shared by both plotters
    data_df = kde.load_data()
    if data_df is None:
        return f"ERROR: {species} file empty or corrupt."

    # Generate plots
    kde.run(data=data_df)
    heat.run(data=data_df)

    # Archive results to summary directory
    dest_summary = os.path.join(summary_dir, root)
    for plotter in [kde, heat]:
        if plotter.output_path and os.path.exists(plotter.output_path):
            shutil.copy2(plotter.output_path, os.path.join(dest_summary, plotter.output_name))

    return f"SUCCESS: {species} in {root}"

def main() -> None:
    """
    Main orchestration function to build tasks and execute them in parallel.
    """
    traj_roots: List[str] = ["alpha", "beta", "delta", "gamma"]
    summary_dir: str = "CV_plots"
    
    # Generate species identifiers
    dist_species: List[str] = [f"d{x}" for x in range(1, 9)]
    tors_species: List[str] = [f"t{x}" for x in range(1, 7)]
    ang_species: List[str] = [f"a{x}" for x in range(1, 7)]
    rad_species: List[str] = [f"rg{x}" for x in range(1, 7)]

    # Assemble task list
    tasks: List[Tuple[str, str, str]] = []
    for root in traj_roots:
        os.makedirs(os.path.join(summary_dir, root), exist_ok=True)
        all_species = dist_species + tors_species + ang_species + rad_species
        for species in all_species:
            tasks.append((species, root, summary_dir))

    print(f"Starting parallel processing of {len(tasks)} jobs...")

    missing_files: List[str] = []
    corrupted_files: List[str] = []
    success_files: List[str] = []
    
    # Execute parallel pool
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = [executor.submit(process_species, task) for task in tasks]

        for future in tqdm(as_completed(futures), total=len(tasks), desc="Plotting Progress"):
            result = future.result()
            if result.startswith("SKIP"):
                missing_files.append(result)
            if result.startswith("ERROR"):
                corrupted_files.append(result)
            if result.startswith("SUCCESS"):
                success_files.append(result)

    if missing_files:
        print(f"\nSummary: {len(missing_files)} files were missing and skipped.")
        print(missing_files)
    if corrupted_files:
        print(f"\nSummary: {len(corrupted_files)} files were empty or corrupted.")
        print(corrupted_files)

    print(f"\nSuccessfully generated and copied {len(success_files)}.")
    print(success_files)

if __name__ == "__main__":
    start_time = time.time()
    main()
    duration = time.time() - start_time
    print(f"\nTotal execution time: {duration:.2f} seconds")
