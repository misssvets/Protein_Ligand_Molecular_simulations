MD Simulation Protocols for GROMACS
1. MD Simulation and Analysis for Protein Only

This protocol outlines the steps for a complete molecular dynamics (MD) simulation of a single protein using the GROMACS software package.

Prerequisites 🔑

    Access to an HPC (High-Performance Computing) cluster.

    GROMACS installed and configured on the HPC.

Simulation Steps 🏃

    Process the PDB File: Use gmx pdb2gmx to prepare your protein's PDB file for simulation. This command generates a processed structure (

    .gro), a topology file (.top), and a position restraints file (.itp). 

```Bash
gmx pdb2gmx -f protein.pdb -o protein.gro -water spce -ignh
```
Select option 9 for protein simulation 

Define the Simulation Box: Use gmx editconf to center the protein and define a cubic box around it. 
```
Bash

gmx editconf -f protein.gro -o protein.gro -c -d 1.0 -bt cubic
```

Add Solvent: Fill the box with water molecules using gmx solvate. 
```
Bash

gmx solvate -cp protein_newbox.gro -cs spc216.gro -o solv.gro -p topol.top
```

Add Ions: Add ions to neutralize the system. First, download 

ions.mdp and use 

gmx grompp to prepare the run file, then use gmx genion to add the ions. 
```
Bash

gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```
# Select solvent to replace with ions [cite: 37]

Energy Minimization: Minimize the system's potential energy to resolve any clashes between atoms. Download 

minim.mdp , rename it to 

em.mdp, and run the following commands. 
```
Bash

gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

NVT Equilibration: Equilibrate the system's temperature and pressure. Download 

nvt.mdp and run the 

gmx grompp and gmx mdrun commands. 
```
Bash

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
```

NPT Equilibration: Proceed to NPT equilibration, which stabilizes the pressure. Download 

npt.mdp and run the commands. 
```
Bash

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt
```
Production MD: Run the final production simulation. Download 

md.mdp , change the 

nsteps for a 100 ns simulation , and execute the run. 
```
Bash

    gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1000
    gmx mdrun -v -deffnm md
```
Analysis Steps 📊

    Center and Remove PBC: Use gmx trjconv to center the protein and remove periodic boundary conditions (PBC) for a clean trajectory. 
```
Bash

gmx trjconv -f md_0_1.xtc -s md_0_1.tpr -pbc mol -ur compact -center -o nopbc.xtc
```
Select group 1 (Protein) and 0 (System) 

Calculate RMSD: Compute the Root Mean Square Deviation to assess the protein's stability over time. 
```
Bash

gmx rms -f nopbc.xtc -s md_0_1.tpr -o native_rmsd.xvg -tu ns
```
Select group 4 (backbone) 

Calculate RMSF: Analyze the Root Mean Square Fluctuation to identify flexible regions of the protein. 
```
Bash

gmx rmsf -f nopbc.xtc -s md_0_1.tpr -o native_rmsf.xvg -res
```
Select group 1 (Protein) 

Calculate Radius of Gyration: Use gmx gyrate to measure the compactness of the protein. 
```
Bash

gmx gyrate -f nopbc.xtc -s md_0_1.tpr -o native_gyrate.xvg
```
Select group 1 (Protein) 

Calculate SASA: Determine the Solvent Accessible Surface Area (SASA) to see how the protein's surface area changes. 
```
Bash

gmx sasa -f nopbc.xtc -s md_0_1.tpr -o native_sasa.xvg -tu ns
```
Select Type: 1 

PCA Analysis: Perform Principal Component Analysis to identify the most significant motions of the protein. 
```
Bash

gmx covar -f nopbc.xtc -s md_0_1.tpr -o eigenval.xvg -v eigenvec.trr -xpma covar.xpm
gmx anaeig -f nopbc.xtc -s md_0_1.tpr -first 1 -last 2 -2d pca-native.xvg
```
Select group 1 (Protein) for both commands 

Potential Energy: Analyze the potential energy of the system to ensure stability. 
```
Bash

gmx energy -f md_0_1.edr -s md_0_1.tpr -o potential.xvg
```
Select 11 (Potential energy) 

Create Movie File: Generate a PDB movie file from the trajectory. 
```
Bash

gmx_mpi trjconv -f nopbc.xtc -s md_0_1.tpr -o movie.pdb -dt 500
```
select protein 

