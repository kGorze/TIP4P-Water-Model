#!/bin/bash

# Create a single TIP4P water molecule
cat > single_water.gro << EOF
Single TIP4P water molecule
    4
    1SOL     OW    1   0.000   0.000   0.000
    1SOL    HW1    2   0.095   0.000   0.000
    1SOL    HW2    3  -0.024   0.093   0.000
    1SOL     MW    4   0.015   0.015   0.000
   1.00000   1.00000   1.00000
EOF

# Convert to PDB format
gmx editconf -f single_water.gro -o single_water.pdb

echo "Created TIP4P water molecule in PDB format" 