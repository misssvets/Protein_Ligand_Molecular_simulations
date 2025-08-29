# MD Simulation and Analysis for Protein-Ligand Complexes

This protocol covers the essential steps for setting up and analyzing a molecular dynamics simulation of a protein bound to a ligand.

Prerequisites 🔑

    Access to an HPC cluster with GROMACS installed.

    Web-based topology builders like SwissParam (Old version for charmm).

Simulation Steps 🏃

    Extract Protein and Ligand: Separate the protein and ligand from the complex PDB file using grep. 


grep UNK ligand.pdb > ligand.pdb
grep ATOM protein.pdb > Protein.pdb

Generate Protein Topology: Use gmx pdb2gmx to create a processed protein structure and a topology file. 

```bash
gmx pdb2gmx -f protein.pdb -o protein.gro -water tip4p -ignh
```
Select option 8 (Charmm36) for protein-ligand simulation 


Generate Ligand Topology: Use a web tool like Swissparam/ CGenFF to generate the ligand's topology (.gro and .itp files). 

    Upload the ligand PDB file (

    ligand.pdb) to the swissparam website. 

Download and unzip the file, then rename the 

.gro and .itp files to drg.gro and drg.itp. 

Merge Protein and Ligand Files: Manually combine the protein and ligand files to create a single system. 

    Open 

    drg.gro and MurD_processed.gro to find the number of atoms in the ligand. 

Add this number to the total atom count in 

protein.gro. 

Copy and paste the atom entries from drg.gro to the bottom of protein.gro. 

Edit the Topology File (topol.top): Add the ligand topology by including the drg.itp file. Also, update the 

[ molecules ] section to include the ligand.     

In your topo.top file, add this line:

#include "drg.itp"

In the [molecules] section, add this:

UNK 1

Add Position Restraints for Ligand: Create a position restraints file for the ligand to keep it from moving during the early stages of the simulation. 

```bash
gmx genrestr -f drg.gro -o posre_drg.itp -fc 1000 1000 1000
```

Select 2 for UNK 

Edit topol.top (again): Add a line to include the new position restraints file for the ligand. 

In your topo.top file, add this line after including drg.itp:
#include "posre_drg.itp"

Merge Protein and Ligand Indices: Create an index file that merges the protein and ligand groups. 
```bash
gmx make_ndx -f em.gro -o index.ndx
```
Select '1 | 13' to merge Protein and UNK (ligand). 
Type 'q' and press Enter to exit. 

Continue with Standard MD Steps as present in the protein format. 

.mdp files and include the index.ndx file in your grompp commands. 

# Analysis Steps 📊

Center and Remove PBC: Use gmx trjconv to clean the trajectory.
```bash
gmx trjconv -f md.xtc -s md.tpr -center -pbc nojump -o nopbc.xtc
```
Select group 22 (Protein_UNK) and 0 (System) 

Calculate RMSD: Compute the RMSD for the complex. 
```bash
gmx rms -f nopbc.xtc -s md.tpr -o rmsd.xvg -n index.ndx -tu ns
```

Select group 22 (Protein_UNK) twice / or backbone 

Calculate RMSF: Measure the fluctuations of the complex. 
```bash
gmx rmsf -f nopbc.xtc -s md.tpr -o rmsf.xvg -n index.ndx -res
```
Select group 22 (Protein_UNK) 

Calculate Radius of Gyration: Check the compactness of the complex. 
```bash
gmx gyrate -f nopbc.xtc -s md.tpr -o gyrate.xvg -n index.ndx
```
Select group 22 (Protein_UNK) 

Analyze Hydrogen Bonds: Analyze the number of hydrogen bonds between the protein and the ligand over time. 
```bash
gmx hbond -f nopbc.xtc -s md.tpr -dist distance.xvg -ac bac.xvg -life hblife.xvg -hbn hbd.xvg  -num hbond.xvg -n index.ndx -tu ns
```
Select group 1 (Protein) and 13 (UNK)

Calculate SASA: Determine the Solvent Accessible Surface Area of the complex. 
```bash
gmx sasa -f nopbc.xtc -s md.tpr -o sasa.xvg -n index.ndx -tu ns
```
Select group 22 (Protein_UNK) 

PCA Analysis: Perform PCA on the protein-ligand complex to identify primary motions. 
```bash
gmx covar -f nopbc.xtc -s md.tpr -o eigenval.xvg -v eigenvec.trr -xpma covara.xpm -l covar.log -n index.ndx
```

Select group 22 (Protein_UNK)
```bash
gmx anaeig -f nopbc.xtc -s md.tpr -first 1 -last 2 -2d pca.xvg -n index.ndx
```
Select group 22 (Protein_UNK)

Potential Energy: Analyze the potential energy of the complex to ensure stability. 
```bash
gmx energy -f md.edr -s md.tpr -o potential.xvg
```
Select 11 (Potential energy) 

