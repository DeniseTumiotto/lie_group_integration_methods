# Numerical time integrators on Lie groups
The original code and content of this repository was developed by Stefan Hante and can be found at [StHante Homepage](https://github.com/StHante).
The code in the present repository has been extended and/or modified.

The repository is organized as follow.

Integrators:

 * The Lie group generalized-α method `gena` (Stefan Hante)
 * The Lie group BDF method `BLieDF` (Stefan Hante)
 * The Lie group RATTLE method `RATTLie` (Stefan Hante)
 * The Lie group SHAKE method `SHAKELie` (Stefan Hante)
 * The nonholonomic RATTLie method `RATTLie_nonhol` (Stefan Hante)
 * The half-explicit Runge-Kutta Lie group methods `half_explicit`

Test problems:

 * The (multiple) heavy top example `heavy_top`
 * The constrained Cosserat beam model `crmS3R3`
 * The rolling disk example `rolling_disk`

Miscellaneous:

 * Implementation of Lie group functions `liegroup`
 * Expand a config file with different configurations to several files `expandconfig`
 * Read lua files in Matlab and Octave `readLua`

Third party projects:

 * [Reading lua files in Fortran `aotus`](https://geb.sts.nt.uni-siegen.de/doxy/aotus/)
 * [GFortran](https://gcc.gnu.org/fortran/)
 * [GNU Parallel](https://www.gnu.org/software/parallel/)
