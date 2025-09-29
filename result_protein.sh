#!/bin/bash

# --- GROMACS MD ANALYSIS SCRIPT ---
# This script automates common post-MD analysis steps.
# Assumes md.xtc, md.tpr, and md.edr files are in the current directory.

# Exit on first error
set -e

# Define input files
TRAJ="md.xtc"
TOPOL="md.tpr"
ENERGY="md.edr"
OUTPUT_DIR="analysis_results"

# --- 0. Setup and Cleaning ---
echo "--- Setting up workspace ---"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
# Link original files for cleaner script structure (optional but good practice)
ln -sf ../"$TRAJ" .
ln -sf ../"$TOPOL" .
ln -sf ../"$ENERGY" .

# --- 1. Center and Remove PBC (Periodic Boundary Conditions) ---
echo "--- 1. Centering and Removing PBC ---"
# Input: 1 (Protein) for center, 0 (System) for output
printf "1\n0\n" | gmx trjconv -f "$TRAJ" -s "$TOPOL" -pbc mol -ur compact -center -o nopbc.xtc

# --- 2. Calculate RMSD (Root Mean Square Deviation) ---
echo "--- 2. Calculating RMSD ---"
# Input: 4 (backbone) for alignment, 4 (backbone) for output
printf "4\n4\n" | gmx rms -f nopbc.xtc -s "$TOPOL" -o native_rmsd.xvg -tu ns

echo "--- 2.1. Calculating RMSD protien ---"
# Input: 4 (backbone) for alignment, 4 (backbone) for output
printf "1\n1\n" | gmx rms -f nopbc.xtc -s "$TOPOL" -o protein_rmsd.xvg -tu ns

# --- 3. Calculate RMSF (Root Mean Square Fluctuation) ---
echo "--- 3. Calculating RMSF ---"
# Input: 1 (Protein)
printf "1\n" | gmx rmsf -f nopbc.xtc -s "$TOPOL" -o native_rmsf.xvg -res

# --- 4. Calculate Radius of Gyration (Rg) ---
echo "--- 4. Calculating Radius of Gyration ---"
# Input: 1 (Protein)
printf "1\n" | gmx gyrate -f nopbc.xtc -s "$TOPOL" -o native_gyrate.xvg

# --- 5. Calculate SASA (Solvent Accessible Surface Area) ---
echo "--- 5. Calculating SASA ---"
# Input: 1 (Protein) for area calculation
printf "1\n" | gmx sasa -f nopbc.xtc -s "$TOPOL" -o native_sasa.xvg -tu ns 

# --- 6. Principal Component Analysis (PCA) ---
echo "--- 6. Performing PCA (Covariance Matrix) ---"
# Input: 1 (Protein)
printf "1\n1\n" | gmx covar -f nopbc.xtc -s "$TOPOL" -o eigenval.xvg -v eigenvec.trr -xpma covar.xpm

echo "--- Performing PCA (Anaeig) ---"
# Input: 1 (Protein)
printf "1\n" | gmx anaeig -f nopbc.xtc -s "$TOPOL" -first 1 -last 2 -2d pca-native.xvg

# --- 7. Potential Energy Analysis ---
echo "--- 7. Analyzing Potential Energy ---"
# Input: 11 (Potential energy)
printf "11\n" | gmx energy -f "$ENERGY" -s "$TOPOL" -o potential.xvg

# --- 8. Create Movie File (PDB) ---
echo "--- 8. Generating PDB Movie File ---"
# Input: 1 (Protein)
printf "1\n" | gmx trjconv -f nopbc.xtc -s "$TOPOL" -o movie.pdb -dt 500

echo "--- 9. Gmx sham ---"
# Input: 1 (Protein)
printf "1\n" | gmx sham -f pca-native.xvg -ls output.xpm -g minima.log -notime

# --- 2. AUTOMATED PARSING OF THE GLOBAL MINIMUM FRAME ---
# This is the crucial step. We find the frame number associated with the
# lowest free energy structure from the minima.log file.
# Note: This logic assumes the minima.log clearly lists the frame index.
# The common pattern is to find the minimum point's index number.

echo "--- 2. Extracting Frame Index of Global Minimum ---"

# This complex line uses AWK to find the *lowest free energy value* in the log
# and extracts the corresponding frame index.
# (Assuming columns are Frame | PC1_Coord | PC2_Coord | Free_Energy)
MIN_FRAME=$(
    awk '/Global/ {next} {if($4 < min || min == "") {min=$4; frame=$1}} END{print frame}' minima.log
)

if [ -z "$MIN_FRAME" ]; then
    echo "ERROR: Could not automatically find the minimum energy frame index in minima.log."
    exit 1
fi

echo "SUCCESS: Global minimum found at Frame Index: $MIN_FRAME"

# --- 3. CONVERSION TO ACTUAL DATA (.pdb) ---
echo "--- 3. Extracting PDB structure for Frame $MIN_FRAME ---"
# gmx trjconv uses the -dump flag to extract a single frame index.
# We pipe the selection (SELECTION) to the trjconv command.

printf "%s\n" "$SELECTION" | gmx trjconv -s "$TOPOL" -f "$TRAJ" -o minimum_energy_structure.pdb -dump "$MIN_FRAME"

echo "--- SUCCESS: Minimum energy structure saved to minimum_energy_structure.pdb ---"

echo "--- Analysis Complete! Results are in the '$OUTPUT_DIR' folder. ---"
