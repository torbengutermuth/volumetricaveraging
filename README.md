# Volumetric Averaging
This repository contains software that can be utilized to analyse how chemical patterns in molecules are distributed in three-dimensional space.
A typical application will be the analysis of the placement of certain functional groups after the docking of a library of small molecules to a protein binding site.

An SD File has to be supplied containing the molecules to be analysed and the results can
be shown as dummy-atoms in a PDB file or in the [GridDataFormat](https://griddataformats.readthedocs.io/en/latest/index.html).
Both File Formats are readable with current molecular viewers like [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)
or [Pymol](https://pymol.org/2/).

## Installation 
The Software was developed and tested on Ubuntu 20.04.2 LTS und should be executable on any standard computer.
The hardware requirements are dependent on the grid size chosen and the author recommends to not use 
grid sizes below 0.4 Angstrom. 
The install time is fully dependent on anaconda, to install the software you just have to download
the repository.

### Set up your own environment with conda
The repository contains a file named `scooder_dev.yml`. It specifies all requirements for the
*conda* environment to run the provided code.

For information on installation on different operating systems or troubleshooting, please refer to
the [Anaconda site](https://www.anaconda.com/distribution/#linux).

### Create a conda environment
After installing conda, the `scooder_dev.yml` file comes into play. It lists all necessary packages needed for conda to install. Assuming you
set up your conda installation properly, you can create a new environment with all required packages
like this:

```bash
conda env create -f scooder_dev.yml
```

### Activate the environment
As soon as the virtual environment has been created properly (see last step) you can activate it with:

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
There are three different modes for the volumetric averaging software: query, nosc and analysis.
In addition to the description of all options in the Command Line Interface, usage examples for the
different modes are provided. 
An SD File for trying out these usage examples can be found in the testdata directory.
### Query
The query mode allows one to directly calculate the grids/densities corresponding to a docking (SD File).
To use it, an SD File has to be specified as well as the grid size, features and the desired output.
Using the command line call below, the PDBs and densities of all features in the extended feature
list are calculated for a grid size of 0.5 Angstrom.
This calculation should take a couple seconds and a folder will be created containing the calculated
PDBs/densities.
```bash
python /pathtovolav/VolumetricAveraging.py query -f 05_moe_prepared.sdf --features extended -g 0.5 --pdb --density
```
### Nosc
The nosc mode enables a user to only calculate the grids that are necessary for calculating the PDBs and densities.
These grids are saved to disk so they can be analysed at a later time.
This calculation should take a couple seconds and a folder will be created containing the calculated
grids in the npy format including a json with settings.
```bash
python /pathtovolav/VolumetricAveraging.py nosc -f 05_moe_prepared.sdf -g 0.5 --features extended
```
### Analysis
The analysis mode allows a user to calculate the PDBs and densities for prior nosc runs. 
This is especially useful to combine the grids of multiple different docking runs, 
saving the time to tediously combine SDFs to analyse subsets of all docking runs conducted.
To combine nosc runs for multiple docking calculations, both the features and the grid size used must be identical.
This calculation should take a couple seconds and a folder will be created containing the calculated
PDBs/densities.
```bash
python /pathtovolav/VolumetricAveraging.py analysis --folders 05_moe_prepared --pdb --density
```