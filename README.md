# Molecular Dynamics Simulation of Nano-Sphere Endocytosis

## Overview
This repository contains all the necessary code and supplementary information to perform molecular dynamics simulations of a nano-sphere being endocytosed. The project is based on the findings published in "Valency of Ligand−Receptor Binding from Pair Potentials" by William Morton, Robert Vácha, and Stefano Angioletti-Uberti, showcasing a novel approach to understanding the dynamics of nanoparticle uptake by cell membranes through ligand−receptor interactions.

## Citation
If you use any part of this code or the supplementary materials in your research, please cite our paper:
> Morton, W., Vácha, R., & Angioletti-Uberti, S. (2024). Valency of Ligand−Receptor Binding from Pair Potentials. Journal of Chemical Theory and Computation. https://doi.org/10.1021/acs.jctc.4c00112

![Abstract Image](/Publication_Images/graphical_abs.png)

## Content Description
- **Code**: Contains the LAMMPS scripts for running molecular dynamics simulations of nano-sphere endocytosis. These scripts are configured based on the coarse-grained model detailed in our publication.
- **Supplementary Information**:
  - **Alternative Potentials**: Lists the different potentials used in our simulations, providing a broader understanding of the ligand−receptor interactions.
  - **Graphing Notebook**: A Jupyter notebook is included for graphing the potentials, aiding in the visualization and analysis of the simulations' results.

## Getting Started
1. **Prerequisites**: Ensure you have [LAMMPS](https://lammps.sandia.gov/) installed on your system for running the simulations.
2. **Running Simulations**: Edit _start.py to choose your parameters and then run the LAMMPS '.in' file. 

## Support
For any questions or issues, please open an issue on this GitHub repository, and we will try to assist you as soon as possible.

## License
This project is licensed under the CC-BY 4.0. See the LICENSE file for details.

## Acknowledgements
We thank the contributors and reviewers for their valuable input and discussions. The simulations were performed at the Imperial College Research Computing Service, [DOI:10.14469/hpc/2232](https://doi.org/10.14469/hpc/2232).
