# Volumetric Averaging
This repository contains software that can be utilized to analyse how chemical patterns in molecule files are distributed in three-dimensional space.
Using this Software you can analyse how different chemical substructures are distributed in
three-dimensional space.
An SD File has to be supplied containing the molecules to be analysed and the results can
be shown as dummy-atoms in a PDB file or in the [GridDataFormat](https://griddataformats.readthedocs.io/en/latest/index.html).
Both File Formats are readable with current molecular viewers like [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)
or [Pymol](https://pymol.org/2/).

## Set up your own environment with conda
The repository contains a file named `scooder_dev.yml`. It specifies all requirements for the
*conda* environment to run the provided code.

For information on installation on different operating systems or troubleshooting, please refer to
the [Anaconda site](https://www.anaconda.com/distribution/#linux).

### Create a conda environment
After installing conda the `scooder_dev.yml` file comes into play. It lists all necessary packages needed for conda to install. Assuming you
set up your conda installation properly, you can create a new environment with all required packages
like this:

```bash
conda env create -f scooder_dev.yml
```

### Activate the environment
As soon as the virtual environment is created properly (see last step) you can activate it with:

```bash
conda activate scooder_dev
```
### Testing if everything works
After successfully creating and activating the environment, you can test that everything works with:
```bash
python -m unittest tests
```
in the volumetric averaging directory.
## Usage Examples
There are three different modes for the volumetric averaging software, query, nosc and analysis.
In addition to the description of all options in the Command Line Interface, usage examples for the
different modes are provided.
### Query
The query mode allows to directly calculate the grids/densities corresponding to a docking (SD File).
To use it an SD File has to be specified as well as the grid size, features and the wanted output.
Using the command line call below, the PDBs and densities of all features in the extended feature
list are calculated for a grid size of 0.5 Angstroem.
```bash
python /pathtovolav/VolumetricAveraging.py query -f 05_moe_prepared.sdf --features extended -g 0.5 --pdb --density
```
The nosc mode allows to only calculate the grids that are necessary for calculating the PDBs and densities.
These grids are saved to disk so they can be analysed at a later time.
```bash
python /pathtovolav/VolumetricAveraging.py nosc -f 05_moe_prepared.sdf -g 0.5 --features extended
```
The analysis mode allows to calculate the PDBs and densities for prior nosc runs. 
This is especially useful to combine the grids of multiple different docking runs, 
saving the time to tediously combine SDFs to analyse subsets of all docking runs conducted.
To combine nosc runs for multiple Dockings, both the features and the grid size used must be identical.
```bash
python /pathtovolav/VolumetricAveraging.py analysis --folders 05_moe_prepared --pdb --density
```