# findrings

Compilation:

gfortran -O3 -ffree-form -ffast-math findrings.f -o findrings

Input files:
   findrings.control // control file
   input.xyz // xyz input file. Name can be set in control file.

findrings.control:
   line1: name of input file (string)
   line2: convergence radius for ring-finding. (12 A is usually sufficient.)

input.xyz:
   first line: number of atoms
   second line: 'ANGLE' (keyword), amag, bmag, cmag, alpha, beta, gamma
      where amag,bmag,cmag,alpha,beta,gamma are the unit cell vectors in Angstroms
   next natomns lines: atomic_name ("O" or "H"), x,y,z (Angstroms)
      assume an order O1,H1a,H1b,O2,H2a,H2b,etc.

output:
    Should print out a list of the number of rings per unit cell by size of ring.
    Will also print ring indices and coordinates to ringx.dat, where x = 3..10

*******
NOTES
*******

The basic problem of identifying rings in a H-bonded graph is complicated somewhat by periodic boundary conditions, which,
for small unit cells in particular, allows rings to 'interfere' with each-other, making it impossible from the graph alone to
correctly count the rings that exist in the periodic system.

To get round this, FINDRINGS.F will attempt to replicate the unit cell, to create a supercell, large enough that this interference
effect is no longer a problem.  The size of the supercell is determined by the tolerance radius in the control file, with the
supercell adjusted to be large enough to encompass a sphere defined by this radius. In practice, a tolerance radius of about
12 Angstroms is sufficient, but it should be checked that the ring-counting is properly converged such that whatever tolerance
radius is chosen, increasing the radius further doesn't change the ring-count.

Although the algorithm tries to find all the rings in the supercell, the ring-count is given *per unit cell*, that is, it divides
its total number of rings for each ring-size by the number of unit cells in the supercell.

It's unfortunately not easy with this approach to identify the indices and coordinates of the *unique* rings for each ring-size.
The files ringx.dat for x = 3..10 give the indices and coordinates for *every* ring in the supercell, which means that, for
small unit cells, it's likely to contain many duplicate rings.

It may be thought that it should be easy to identify duplicates, but, for small unit cells, two rings can have exactly the same
indices, but represent structurally distinct rings. Again, this is due to the periodic boundary conditions problem. However,
users are welcome to try modifying the code, or writing their own code to find the unique rings from the output.
