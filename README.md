# md_protein
MOlecular dynamics data processing and analysis 
gmx convert-trj -f nopbc.xtc -s md.tpr -o md.dcd -n index.ndx
gmx distance -f nopbc.xtc -s md.tpr -n index.ndx -oall dis4.xvg
gmx angle -f nopbc.xtc -n index.ndx -type dihedral -od angle4.xvg
