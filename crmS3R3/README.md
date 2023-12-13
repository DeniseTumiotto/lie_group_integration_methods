# Cosserat beam model discretized in the semidirect Lie group 𝕊³⋉ℝ³

This project can be used to simulate a Cosserat beam model. It may be constrained to an (extensible) Kirchhoff beam.
The Cosserat beam is modeled to be, at each time instance 𝑡, a curve 𝑠↦𝑞(𝑠,𝑡) in the semidirect product Lie group 𝕊³⋉ℝ³. Here, 𝕊³ is the Lie group of unit quaternions representing orientation (more precisely: rotation from a reference orientation). For each time 𝑡, 𝑞(𝑠,𝑡)=(𝑝(𝑠,𝑡),𝑥(𝑠,𝑡)) represents the orientation 𝑝(𝑠,𝑡)∈𝕊³ and position 𝑥(𝑠,𝑡)∈ℝ³ of a cross section at the spatial parameter 𝑠∈[0,𝐿] of the beam, where 𝐿 is the length of the beam. Each cross section is assumed to stay rigid.

For the time integration, one of the following three integration methods can be used:

 * A generalized-α Lie group method [`gena`](https://github.com/StHante/gena)
 * A Lie group generalization of the RATTLE method [`RATTLie`](https://github.com/StHante/RATTLie) (In the case of an ODE, RATTLE is just the Störmer-Verlet method.)
 * A generalization of the BDF method to (constrained) differential equations of second order on Lie groups [`BLieDF`](https://github.com/StHante/BLieDF).

The implementation of this project is done in modern Fortran. It was only tested with `gfortran` on Linux, but ports to different Fortran compilers and platforms should be easily possible.

## Prerequisites
In order to build this project the following other projects are required. Make sure that they can be found by the makefile.
This need the following projects:

 * For the time integration, at least one of [`gena`](https://github.com/StHante/gena), [`RATTLie`](https://github.com/StHante/RATTLie) or [`BLieDF`]https://github.com/StHante/BLieDF) is needed. The path and the name of the integrator should be saved as the variable `INTEGRATORP` and `INTEGRATOR` in the makefile.
 * [`liegroup`](https://github.com/StHante/liegroup): This is a project that implements a lot of the functions related to quaternions and the Lie groups 𝕊³ and 𝕊³⋉ℝ³. Use `git clone https://github.com/StHante/liegroup.git` in order to clone this project. The path should be saved as the variable `LIEFUNP` in the makefile.
 * [`aotus`](https://geb.sts.nt.uni-siegen.de/doxy/aotus/): This project makes it possible to read lua files in Fortran. Use `hg clone https://hg.osdn.net/view/apes/aotus` in order to clone this project by using the mercurial command `hg`. The path should be saved as the variable `AOTP` in the makefile.
 * [`expandconfig`](https://github.com/StHante/expandconfig): This project is used to preprocess the lua-files. Use `git clone https://github.com/StHante/expandconfig.git` in order to clone this project. The path should be saved as the variable `EXPANDCONFIG` in the makefile.
 * [`readLua`](https://github.com/StHante/readLua-for-Matlab-and-Octave): This project makes it possible to read lua files in Matlab. It is not needed to compile the executable, but is used to analyze the test results in Matlab. Use `git clone https://github.com/StHante/readLua-for-Matlab-and-Octave.git` in order to clone this project. The path should be saved as the variable `READLUAP` in the makefile. Also make sure that the variable `MATLAB` is a valid command that starts Matlab. (Note that in order to compile `readLua` your Matlab must be set up to compile mex.) If you can not get this to work, you can manually download the `.m` and `.mex*` files from the repository and copy them to `test/als`.

## Usage
After getting all prerequisites, the command `make` can be used to build the executable. There are a few special goals in the makefile:

 * `make try`: This will build the executable and execute the first test that is described in `test/config.lua`.
 * `make test`: Builds the executable and executes all tests described in `test/config.lua`.
 * `make cleanup`: Removes all files within the directory `obj`.
 * `make clean`: Like `cleanup`, but also removes the executables.
 * `make cleantest`: Removes all test results.
 * `make mrproper`: Like `clean` and `cleantest` together.

Note that the integrator defined in the makefile can be overridden by adding `INTEGRATOR=<integratorname>` and `INTEGRATORP=<integratorpath>` to the `make` command. Also note that before changing the integrator, `make cleanup` is necessary and `make cleanintegrator` is recommended.

## The configuration file
The main configuration file is `test/config.lua`. It is a lua file that additionally may use `expandconfig` syntax. (Short example: `a= [[ 0 || 1 ]]` in the configuration file would result in two expanded configuration files with `a=  0` and `a=  1`. When commenting out such a line, remember to replace `[[` by `[--[` or something similar otherwise you will get a lot of expanded files that only differ in comments.)
It follows a list of all relevant options with short explanation. Refer also to the documentation of the integrators, where some of the options concerning the integrator are described in more detail. Note that not all integrator options apply to all integrators.
Besides the file `test/config.lua`, there are some other lua files in `test` which can serve as templates for lua configuration files.

 * `t0`: The initial value for the time.
 * `te`: End time of the integration.
 * `steps`: Amount of steps between `t0` and `te` to be taken. The step size can be calculated by `h=(te-t0)/steps`. below.
 * `problem_name`: A string that describes the problem that is given in the configuration file. This will not influence the integration and is just useful for the analysis of the results.
 * `kirchhoff`: Set this to `1`, if the full Cosserat beam model should be constrained to an extensible Kirchhoff beam model. In that case, the problem becomes constrained. Set this to `0` to use the full Cosserat beam model. In that case, the problem is unconstrained.
 * `nr_subdiag`, `nr_superdiag`: Refer to the documentation of the integrator. This is should be calculated automatically in the configuration file, depending on the constraints.
 * `n`: Number of spatial discretization intervals. The beam will be discretized in space using the equidistant discretization points 𝑞₀(𝑡), 𝑞₁(𝑡), …, 𝑞ₙ(𝑡).
 * `L`: Length of the beam.
 * `CGam`: Array with three elements containing the diagonal of the stiffness matrix corresponding to the material strain vector Γ. Often, `CGam = { 𝐺𝐴₁, 𝐺𝐴₂, 𝐸𝐴 }` and sometimes, 𝐺𝐴ₖ is the product of the G-modulus 𝐺, the area of cross section 𝐴 and the Timoshenko shear correction factor ϰₖ, 𝐸𝐴 the product of the E-modulus 𝐸 and the area 𝐴.
 * `CK`: Array with three elements containing the diagonal of the stiffness matrix corresponding to the material curvature vector 𝐾. Often, `CK = { 𝐸𝐼₁, 𝐸𝐼₂, 𝐺𝐼₃ }` and sometimes, 𝐸𝐼ₖ is the product of 𝐸 and 𝐼ₖ and 𝐺𝐼₃ is the product of 𝐺 and 𝐼₃, where 𝐼ₖ are the diagonal elements of the inertial tensor density.
 * `CGamd`: Array with three elements containing the diagonal of the dissipation (damping) matrix corresponding to the time derivative ̇Γ of material strain vector.
 * `CKd`: Array with three elements containing the diagonal of the dissipation (damping) matrix corresponding to the time derivative ̇𝐾 of material curvature vector.
 * `ds`: Length of a discretization interval. Usually `ds=L/n`. This is just for computing following parameters.
 * `m`: Mass of a beam segment. Usually `m = A*rho*ds`, where `A` is the area of the cross section of the beam and `rho` is the beam's density.
 * `mI`: Diagonal elements of the inertial tensor of a beam segment. Usually, the inertial tensor is the product of the density `rho`, the length `ds` and the inertial tensor density of the beam.
 * `p0`: This should be a function of one scalar parameter `sigma` and `p0(sigma)` should return the orientation of the cross section of the beam as a unit quaternion (array with four entries) at the spatial parameter `sigma*L`. Note that `sigma` must be between `0` and `1` (inclusive), whereas the spatial parameter 𝑠∈[0,𝐿].
 * `x0`: This should be a function of one scalar parameter `sigma` and `x0(sigma)` should return the position of the cross section of the beam (array with three entries) at the spatial parameter `sigma*L`. Note that `sigma` must be between `0` and `1` (inclusive), whereas the spatial parameter 𝑠∈[0,𝐿].
 * `Om0`: This should be a function of one scalar parameter `sigma` and `Om0(sigma)` should return the angular velocity of the cross section of the beam (array with three entries) at the spatial parameter `sigma*L`. Note that the angular velocity must be expressed with respect to the body-fixed frame. Note that `sigma` must be between `0` and `1` (inclusive), whereas the spatial parameter 𝑠∈[0,𝐿].
 * `V0`: This should be a function of one scalar parameter `sigma` and `V0(sigma)` should return the velocity of the cross section of the beam (array with three entries) at the spatial parameter `sigma*L`. Note that the velocity must be expressed with respect to the body-fixed frame. Note that `sigma` must be between `0` and `1` (inclusive), whereas the spatial parameter 𝑠∈[0,𝐿].
 * `external`, `external_parameters`: With these two variables, the external generalized forces (moments and forces) that are applied to the beam are described. Note that all generalized forces that are not due to the internal deformation of the beam are considered "external"; this applies e.g. to gravity. The variable `external` should be a specific string and `external_parameters` a table with specific keys depending on `external`. Here is a list of all possible values of `external` with description of all associated keys of `external_parameters`:
   * `external='gravity'`: This will apply gravity to the beam. In this case, `external_parameters` only needs to have one key: `g`, that should hold the amount of gravity that acts in the positive z-direction. (Often `g=-9.81`.)
   * `external='roll-up'`: This will apply a moment around y-axis to the last beam element and its negative to the first beam element. In this case, `external_parameters` only needs to have one key: `factor`, that should hold amout of moment. (Often `factor=math.pi*CK[1]/L` (in case the beam is free-free) or `factor=2*math.pi*CK[1]/L` (in case the beam is clamped at exactly one end), because in that case, the beam will roll-up to a full circle (by applying heavy damping). This test would fail if locking is present, which is not the case for this model.)
   * `external='flying_spaghetti'`: This will apply a moment and force to the first beam element. We require that `t0=0`. The moment and force will start out as zero vectors, then increase linearly for a time span of `increasing_time` until they reach `maximum_height*moment_factors` and `maximum_height*force_factors` respectively, then decrease linearly for a time span of `decreasing_time` until they vanish. After that, no more moments and forces are applied. The keys of `external_parameters` are the aforementioned `increasing_time`, `decreasing_time`, `maximum_height`, `moment_factors` and `force_factors`.
 * `fixed_x0`: Set this to `0` if the left end of the beam is free to move translationally. Set this to `1`, if the left end of the beam should be fixed.
 * `fixed_x0_position`: If `fixed_x0=1` this should be the position that the left end is supposed to be fixed at. Usually `fixed_x0_position=x0(0)` is a good choice.
 * `fixed_xn`: Set this to `0` if the right end of the beam is free to move translationally. Set this to `1`, if the right end of the beam should be fixed.
 * `fixed_xn_position`: If `fixed_xn=1` this should be the position that the right end is supposed to be fixed at. Usually `fixed_xn_position=x0(1)` is a good choice.
 * `fixed_p0`: Set this to `0` if the left end of the beam is free to rotate. Set this to `1`, if the left end of the beam should be fixed rotationally.
 * `fixed_p0_orientation`: If `fixed_p0=1` this should be the orientation that the left end is supposed to be fixed at. Usually `fixed_x0_position=p0(0)` is a good choice.
 * `fixed_pn`: Set this to `0` if the right end of the beam is free to rotate. Set this to `1`, if the right end of the beam should be fixed rotationally.
 * `fixed_pn_orientation`: If `fixed_pn=1` this should be the orientation that the right end is supposed to be fixed at. Usually `fixed_xn_position=p0(1)` is a good choice.
 * `output_t_at`: If this is `0`, output will be given after each successful integration step. If this is `1`, output will only be given if after a successful integration step the time `t` is a multiple of `t_output_at_multiples_of`.
 * `t_output_at_multiples_of`: see above
 * `output_s_at`: Set this to `0` to output all configurations of all beam elements. Set this to `1` to interpolate the beam in space and output at defined spatial parameters `sigma*L`, where `sigma` is between `0` and `1`. The `sigma`s used are defined in the array `output_s`.
 * `output_s`: see above.
 * `alpha_m`, `alpha_f`, `beta`, `gamma`: Only applies if the integrator is `gena`. Algorithmic parameters for the generalized-α method `gena`.
 * `k_bdf`: Only applies if the integrator is `BLieDF`. Number of previous steps that the method uses.
 * `const_mass_matrix`: Set this to `1` to assume the mass matrix is constant, which is the case in this example. Otherwise, use `0`.
 * `diag_mass_matrix`: Set this to `1`, since the mass matrix is diagonal. Otherwise, use `0`.
 * `banded_iteration_matrix`: Set this to `1`, in order to exploit the band structure of the problem. Otherwiese, use `0`.
 * `recalc_iteration_matrix`: Set this to `0` if the iteration matrix should not be recalculated after each Newton step. Set this to `1` to do so (increases runtime but might be more robust.)
 * `perturb`: Only applies if the integrator is `gena` or `BLieDF` and only in the index-3 case. Refer to documentation of `gena` and `BLieDF`. Note that this is set automatically if `problemset>0`.
 * `perturb_s`: Only applies if the integrator is `gena` or `BLieDF`. Refer to documentation of `gena` and `BLieDF`. Note that this is set automatically if `problemset>0`.
 * `use_num_K`, `use_num_D`: Set this to `1` since the analytic tangent stiffness and damping matrices are not implemented. They are approximated using finite differences.
 * `no_K`, `no_D`: Set this to `1` to omit tangent stiffness and tangent damping matrices in the iteration matrix. Otherwise, set this to `0`.
 * `rtol`, `atol`: Relative and absolute tolerances of the Newton iteration.
 * `imax`: The number of Newton iteration steps after which the Newton method is considered to not converge.
 * `stab2`: Only applies in the constrained case. Set this to `0` to use the index-3 formulation, set this to `1` to use the stabilized index-2 formulation.

You have to make sure that the initial positions, orientations, angular velocities and velocities must be compatible with the Kirchhoff constraints, if `kirchhoff=1`. Note that fixing the orientation but not the position has never been tested.

## Output of the tests
When `make test` is called, the executable `crm` will be built, the configuration file `config.lua` will be fed through `expandconfig` and the output files (that should be plain lua files) are written to `test/cfg_exp/`. Then, `crm` is executed multiple times with every lua file in `test/cfg_exp/`. Each instance of `crm` will write output to the directory `test/out/`. For each test (ie. lua file in `test/cfg_exp/`) three files will be created:

 * `.bin` file: Containing the output of `crm` as binary data.
 * `.lua` file: Contains the original lua file file that was fed to `crm` with additional details like the time and date of compilation, which integrator was used and some stats like runtime, details of the Newton method and number of calls of some functions.
 * `.misc` file: Usually the time when the output procedure wrote something to the `.bin` file is written to this file in human-readable form. For long running tests, one can use `tail -f test/out/*.misc` to see the ends of all `.misc` files.
 * `.err` file: Everything that is written by `crm` to standard output will be redirected to this file. (Does not apply to `make try`.)

## Analysis of the test results
The analysis of the test results is done in Matlab due to its easy way of manipulating data and visualizing it. All Matlab files are situated in `test/als/`. Here is a list with the most important functions that can and should be used:

 * `sol = load_latest_config_and_bin()`: Will look for the most recent duo of `.bin` and `.lua` files and will load them as a struct containing the most important variables from the lua file and the output from the binary file. The results are then available in the struct as `sol.rslt`, which is itself a struct with fields `t`, `q`, `v` and possibly `l` (for the Lagrange multipliers). Note that this will need `readLua`, see also the section Prerequisites in this readme. Note that you can also load the oldest file by passing the parameter `1`, the second to latest file by passing `-1` and so on.
 * `load_all_config_and_bin()`: Like above but will load all duos of `.bin` and `.lua` files in a cell of structs. You can pass a struct that can act as a pattern (eg. saying `pattern.steps = 10; load_all_config_and_bin(pattern)` will only load test results with exactly 10 steps.)
 * `solcell = calc_errors(solcell, refpattern)`: Calculates the errors with respect to a reference solution. The reference solution should be a part of the cell `solcell` and is determined by the pattern `refpattern`. There should be exactly one solution in `solcell` that matches `refpattern`. The relative and absolute errors are then available in the structs `solcell{k}.err.rel` and `solcell{k}.err.abs`.
 * `delete_all_config_and_bin`, `delete_unfinished_config_and_bin`, `deleta_duplicate_config_and_bin`: These functions will delete all, unfinished, or duplicate duos of `.bin` and `.lua`. Use with caution.
 * `makexyplot(solcell, pattern, xname, yname, byname, varargin)`: Use this to plot stuff. For example `makexyplot(solcell, struct(), 'steps', 'err.abs.q', 'constrained')` would open a figure where the absolute error in 𝑞 is plotted over the number of steps. There will be two lines in the diagram distinguished by the value of `constrained`. (The `varargin` is passed to `plot`.)
 * `makexyyyplot(solcell, pattern, xname, ynames, varargin)`: This function is similar to `makexyplot`, but while `makexyplot` plots the same y-value sorted by a different field, `makexyyyplot` can plot different y-values. Here, `ynames` is supposed to be a cell with strings. For example `makexyyyplot(solcell, struct(), 'steps', {'err.abs.q', 'err.abs.v'})` would open a figure where the absolute error in 𝑞 and in 𝑣 are plotted over the number of steps.
 * `matlab2csv(path, ax, nsteps)`: Will take the axes `ax` (if omitted the current axes `gca`) and convert the data to a csv-file that will be saved to the directory `path`. The parameter `nsteps` will take only export each `nsteps`th data point. (If omitted, `nsteps=1`. Choose higher `nsteps` for plots with a _lot_ of data points.) These csv files can be easily imported with `pgfplots` to produce beautiful plots in LaTeX.
 * `visualize_beam(sol)`: Will produce a moving 3D plot of the solution. (There are some options that are defined inside the function; eg for creating an actual video.)
 * `snapshotplot(sol, ts, stride)`: Creates a 3D plot of snapshots of the beam. The vector `ts` defines the time instances at which snapshots are to be displayed. It also plots the trajectories of the ends of the beam. Increase the natural number `stride` in order to decrease the resolution of these trajectories.
 * `sol = calc_energy(sol)`: Calculates several energies of the system available as `sol.rslt.potential_energy`, `sol.rslt.shearing_energy`, `sol.rslt.torsion_energy`, `sol.rslt.bending_energy`, `sol.rslt.extension_energy`, `sol.rslt.kinetic_energy` and the overall mechanical energy `sol.rslt.energy`.


There are a few more functions that work in the background, but might be very important, eg. `load_config_and_bin` (that might be changed if `crm` is changed) and `structmatch`.


## Source files and what they do
Here is a list with all source files and a short overview of what they contain and do:

 * `main.F90`: This Fortran source file contains the `program main`, which will be the program that is executed, when the compiled executable is invoked. Here is a list of what is done in which order:
   1. Get the first command line argument, which should be the path to a valid lua configuration file.
   3. Load all relevant parameters from the lua configuration file, check for some errors and print them to standard output.
   5. Get the second command line argument, which hsould be the path to a output file. (All folders must exist, they will not be created and the file must not exist.)
   6. Write the lua configuration file to the lua output file.
   7. Write some additional information (time/date of compilation, integrator) to the lua output file.
   8. Open binary and misc output files.
   9. Start the integration.
   9. Close binary and misc output files.
   9. Write statistics to lua output file and standard output.
   9. Clean up the problem object.
 * `crm.F90`: This Fortran source file contains the module `heavy_top` which defines a type `heavy_top_t` which extends the abstract problem type from the integrator. Here, all deferred procedures are implemented. See the documentation of the integrator for details.
 * `get_line_of_variable_length.F90`: Defines a module with the same name that contains a function that can read a line of an opened file without knowing the length of thte line in advance.
 * `makefile`: A makefile that can be executed by calling `make`.
   1. Variables get defined. Modify them to your needs, especially the paths for the different projects and executable names, eg. how to call Matlab.
   2. Goals with prerequisites and how to make them. In this part, only values of the earlier defined variables are used. Usually, no modifications are needed here.

## Preprocessor functions and macros
In the source code, there are some preprocessor functions and macros used. Here is a list of most of them with some explanation:

 * `GL`: At the top, there is the macro `GL` defined as `#define GL(x) x`. The name is short for "glue" and can be used to "glue" the value of preprocessor variables to other text. As an example, let's say there is a preprocessor variable `VAR` and we want to use its value as part of a name. Simply using `VAR_this_text` will not work, because the preprocessor recognizes `VAR_this_text` as one entity. This can be circumvented by saying `GL(VAR)_this_text`. Now, the preprocessor knows that `GL(VAR)` is independent of the rest. If the value of `VAR` would be `example`, then `GL(VAR)_this_text` will become `example_this_text`, whereas `VAR_this_text` would stay `VAR_this_text`.
 * `INTEGRATOR`: This variable contains the name of the integrator (`gena`, `BLieDF`, `RATTLie` etc.)
 * `INT_RATTLie`, `INT_SHAKELie`, `INT_gena`, `INT_BLieDF`, `INT_varint4lie`: These get defined when RATTLie, SHAKELie, BLieDF, gena or varint4lie are used as an integrator. (SHAKELie is seldom used, since RATTLie is superior and varint4lie is unfinished and broken.)
 * `for(i,f,l)` and `endfor`: These two macros are defined in `crm.F90` serve as a wrapper to loops. There are two different scenarios:
   * `crm` is compiled with the preprocessor variable definition `pure=''`. This will make all otherwise `pure` procedures non-`pure`, because `pure` gets replaced by the empty string by the preprocessor. Refer also to the integrator documentation. Since we now expect that all procedures may print data (or have side-effects), all these loops are realized as regular `do` loops. This means, that the preprocessor replaces `for(i,f,l)` by `do i=f,l` and `endfor` by `end do`.
   * Otherwise. Since a lot of procedures are `pure`, meaning they may have no side-effects, we can realize a lot of loops by using Fortran's `forall` construct. This means, that the preprocessor replaces `for(i,f,l)` by `forall(i=f:l)` and `endfor` by `end forall`.
 * `LINEAR_W`: If this preprocessor variable is defined, a linearization of the logarithm in the calculation of the spatial velocity vectors 𝑤 is used. This might be faster, but will probably yield less accurate results.
 * `VARIABLE_STEPS`: Only applies when `RATTLie` is used as an integrator, see also its documentation. If `crm` (and `RATTLie` as well) is compiled with `VARIABLE_STEPS` defined, an additional lua variable `tspan` must be defined in the configuration lua file that holds all time instances. They don't have to be equidistant.
 * Note that `print`ing (or `write`ing) the value of a preprocessor variable is not easily possible, but there is a really ugly way to do it. The following segment will print the value of the preprocessor variable `VARIABLE` to standard output. It makes use of the fact that continuing lines in possible even in strings.
```
print *, '&
VARIABLE'
```

## Related projects
Integrators:

 * [The Lie group generalized-α method `gena`](https://github.com/StHante/gena)
 * [The Lie group BDF method `BLieDF`](https://github.com/StHante/BLieDF)
 * [The Lie group RATTLE method `RATTLie`](https://github.com/StHante/RATTLie)
 * [The Lie group SHAKE method `SHAKELie`](https://github.com/StHante/SHAKELie)
 * [The nonholonomic RATTLie method `RATTLie_nonhol`](https://github.com/StHante/RATTLie_nonhol)

Test problems:

 * [The heavy top example `heavy_top`](https://github.com/StHante/heavy_top)
 * [The constrained Cosserat beam model `crmS3R3`](https://github.com/StHante/crmS3R3)
 * [The rolling disk example `rolling_disk`](https://github.com/StHante/rolling_disk)

Miscellaneous:

 * [Implementation of Lie group functions `liegroup`](https://github.com/StHante/liegroup)
 * [Expand a config file with different configurations to several files `expandconfig`](https://github.com/StHante/expandconfig)
 * [Read lua files in Matlab and Octave `readLua`](https://github.com/StHante/readLua-for-Matlab-and-Octave)

Third party projects:

 * [Reading lua files in Fortran `aotus`](https://geb.sts.nt.uni-siegen.de/doxy/aotus/)
 * [GFortran](https://gcc.gnu.org/fortran/)
 * [GNU Parallel](https://www.gnu.org/software/parallel/)
