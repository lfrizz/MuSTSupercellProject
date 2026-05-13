This code was written to be added to the MuST (Multiple Scattering Theory) code suite, to simulate thermal displacements of atoms in a pure solid or a simple alloy.
The greater MuST code suite can be found at the following link: https://github.com/mstsuite/MuST

Our files can be run using the following commands, assuming one has already fully installed MuST:

1. Place "testSupercell.f90" in the directory "~/MuST/MST/bin".
2. Place "SupercellModule.f90" in the directory "~/MuST/MST/src".
3. Place the file "Makefile" from this github repository in the directory "~/MuST/MST/bin", then run this file using the command "make" in this same directory.

You now have an executable called "disp" that you can run from any directory using the command "~/MuST/MST/bin/disp", which is fully MuST-compatible, and can generate
realisitc thermal displacements for any given temperature (in Kelvin) in seconds.
