#!/bin/bash

# This script runs a continuous GROMACS MD simulation from start to finish.
# Make sure to have your initial MurD.pdb file in the same directory.

# Exit on first error
set -e

# --- 1. Generate topology and position-restraint file ---
echo "Running pdb2gmx..."
printf "8\n" | gmx pdb2gmx -f VanA.pdb -o VanA.gro -water tip3p -ignh

# --- 2. Define the simulation box ---
echo "Running editconf..."
gmx editconf -f VanA.gro -o VanA_newbox.gro -c -d 1.0 -bt cubic

# --- 3. Solvate the system with water ---
echo "Running solvate..."
gmx solvate -cp VanA_newbox.gro -cs tip3p.gro -o solv.gro -p topol.top

# --- 4. Add ions to neutralize the system ---
echo "Running grompp for ions..."
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 10

echo "Running genion..."
printf "13\n" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# --- 5. Energy Minimization ---
echo "Running grompp for energy minimization..."
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1

echo "Running mdrun for energy minimization..."
gmx mdrun -v -deffnm em

# --- 6. NVT Equilibration ---
echo "Running grompp for NVT equilibration..."
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1

echo "Running mdrun for NVT equilibration..."
gmx mdrun -v -deffnm nvt

# --- 7. NPT Equilibration ---
echo "Running grompp for NPT equilibration..."
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1

echo "Running mdrun for NPT equilibration..."
gmx mdrun -v -deffnm npt

# --- 8. Final MD Production Run ---
echo "Running grompp for MD production..."
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1

echo "Running mdrun for MD production..."
gmx mdrun -v -deffnm md
