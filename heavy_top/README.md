# Heavy top problem
This project can be used to simulate a heavy top. There are seven different Lie group formulations implemented:

 1. ODE on SO(3),
 2. ODE on 𝕊³,
 3. DAE on SO(3)×ℝ³,
 4. DAE on 𝕊³×ℝ³,
 5. DAE on SO(3)⋉ℝ³,
 6. DAE on 𝕊³⋉ℝ³,
 7. DAE on unit dual quaternions.

For the time integration, one of the following three integration methods can be used:

 * A generalized-α Lie group method [`gena`](https://github.com/StHante/gena)
 * A Lie group generalization of the RATTLE method [`RATTLie`](https://github.com/StHante/RATTLie) (In the case of an ODE, RATTLE is just the Störmer-Verlet method.)
 * A generalization of the BDF method to (constrained) differential equations of second order on Lie groups [`BLieDF`](https://github.com/StHante/BLieDF).

The implementation of this project is done in modern Fortran. It was only tested with `gfortran` on Linux, but ports to different Fortran compilers and platforms should be easily possible.

## Prerequisites
In order to build this project the follwing other projects are required. Make sure that they can be found by the makefile.
This need the following projects:

 * For the time integration, at least one of [`gena`](https://github.com/StHante/gena), [`RATTLie`](https://github.com/StHante/RATTLie) or [`BLieDF`](https://github.com/StHante/BLieDF) is needed. The path and the name of the integrator should be saved as the variable `INTEGRATORP` and `INTEGRATOR` in the makefile.
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
 * `gena`: This goal itself does nothing, but when it is added to the make goals, the integrator `gena` will be used.
 * `RATTLie`: This goal itself does nothing, but when it is added to the make goals, the integrator `RATTLie` will be used.
 * `BLieDF`: This goal itself does nothing, but when it is added to the make goals, the integrator `BLieDF` will be used.

Note that if none of `gena`, `RATTLie` or `BLieDF` are given, the integrator is determined by what the makefile defines. Also note that before changing the integrator, `make cleanup` is necessary.

## The config file
The main configuration file is `test/config.lua`. It is a lua file that additionally may use `expandconfig` syntax. (Short example: `a= [[ 0 || 1 ]]` in the config file would result in two expanded config files with `a=  0` and `a=  1`. When commenting out such a line, remember to replace `[[` by `[--[` or something similar otherwise you will get a lot of expanded files that only differ in comments.)
It follows a list of all relevant options with short explanation. Refer also to the documentation of the integrators, where some of the options concerning the integrator are described in more detail. Note that not all integrator options apply to all integrators.

 * `t0`: The initial value for the time.
 * `te`: End time of the integration.
 * `steps`: Amount of steps between `t0` and `te` to be taken. The step size can be calculated by `h=(te-t0)/steps`.
 * `problemset`: Here we can easily select the Lie group and formulation to be used to simulate the heavy top. If `problemset=-1` then the options can be chosen indepentently by variables later in the config file. See the table below.
 * `mass`: The mass of the top.
 * `gravity`: Strength and direction of the gravity field as an array with three elements.
 * `inerJ`: The diagonal entries of the inertial tensor as an array with three elements. Since we always describe the angular velocity of the top with respect to the body fixed frame (with principal axes), the inertial tensor must be a diagonal matrix.
 * `x0`: Initial position of the top as an array with three elements.
 * `p0`: Initial orientation of the top as a unit quaternion (array with four elements). Note that `x0` and `p0` have to be chosen in such a way that if we rotate `x0` with the inverse of `p0`, we have to end up with `{0,1,0}`, otherwise the inertial tensor is wrong (or the constraints are violated).
 * `Om0`: The initial angular velocity of the top (measured in the body-fixed frame) as an array with three elements. Note that there is no way to specify the initial velocity (`d/d𝑡 x0`), because it will be calculated accordingly from the constraints.
 * `liegroup`: Which Lie group to use, see enumeration at the very beginning of this document. Note that this is set automatically if `problemset>0`.
 * `output_type`: Can be either `0` for no output or `1` for regular output in terms of the chosen Lie group. Option `2` for output in terms of SO(3)×ℝ³ was not implemented (yet).
 * `output_t_at`: If this is `0`, output will be given after each successful integration step. If this is `1`, output will only be given if after a successful integration step the time `t` is a multiple of `t_output_at_multiples_of`.
 * `t_output_at_multiples_of`: see above
 * `alpha_m`, `alpha_f`, `beta`, `gamma`: Only applies if the integrator is `gena`. Algorithmic parameters for the generalized-α method `gena`.
 * `k_bdf`: Only applies if the integrator is `BLieDF`. Number of previous steps that the method uses.
 * `const_mass_matrix`: Set this to `1` to assume the mass matrix is constant, which is the case in this example.
 * `diag_mass_matrix`: Set this to `0`, since we are working with a full mass matrix. (Although the mass matrix is actually diagonal in this case. Maybe implement this in the future.)
 * `banded_iteration_matrix`: Set this to `0`, since there is no useful band structure of the iteration matrix in any case.
 * `recalc_iteration_matrix`: Set this to `0` if the iteration matrix should not be recalculated after each Newton step. Set this to `1` to do so (increases runtime but might be more robust.)
 * `perturb`: Only applies if the integrator is `gena` or `BLieDF` and only in the index-3 case. Refer to documentation of `gena` and `BLieDF`. Note that this is set automatically if `problemset>0`.
 * `perturb_s`: Only applies if the integrator is `gena` or `BLieDF`. Refer to documentation of `gena` and `BLieDF`. Note that this is set automatically if `problemset>0`.
 * `use_num_K`, `use_num_D`: Set this to `0` to use the implemented tangent stiffness and tangent damping matrices. Set this to `1` to use finite difference approximation (not recommended).
 * `no_K`, `no_D`: Set this to `1` to omit tangent stiffness and tangent damping matrices in the iteration matrix. Otherwise, set this to `0`. Note that this is set automatically if `problemset>0`.
 * `rtol`, `atol`: Relative and absolute tolerances of the Newton iteration.
 * `imax`: The number of Newton iteration steps after which the Newton method is considered to not converge.
 * `stab2`: Set this to `0` to use the index-3 formulation, set this to `1` to use the stabilized index-2 formulation. Note that this is set automatically if `problemset>0`.

| `problemset=`? | SO(3) | 𝕊³ | SO(3)×ℝ³ | 𝕊³×ℝ³ | SO(3)⋉ℝ³ | 𝕊³⋉ℝ³ | unit dual quaternions |
| ------------- | ----- | --- | -------- | ------ | --- | --- | --- | --- |
| unconstrained | 10    | 20 |          |        |     |    |    |     |
| index-3 formulation (no perturbation) | | | 31 | 41 | 51 | 61 | 71 |
| index-3 formulation (with perturbation) | | | 32 | 42 | 52 | 62 | 72 |
| stabilized index-2 formulation (omitting stiffness and damping matrix) | | | 33 | 43 | 53 | 63 | 73 |
| stabilized index-2 formulation | | | 34 | 44 | 54 | 64 | 74 |

## Output of the tests
When `make test` is called, the executable `heavy_top` will be built, the configuration file `config.lua` will be fed through `expandconfig` and the output files (that should be plain lua files) are written to `test/cfg_exp/`. Then, `heavy_top` is executed multiple times with every lua file in `test/cfg_exp/`. Each instance of `heavy_top` will write output to the directory `test/out/`. For each test (ie. lua file in `test/cfg_exp/`) three files will be created:
 
 * `.bin` file: Containing the output of `heavy_top` as binary data.
 * `.lua` file: Contains the original lua file file that was fed to `heavy_top` with additional details like the time and date of compilation, which integrator was used and some stats like runtime, details of the Newton method and number of calls of some functions.
 * `.misc` file: Usually the time when the output procedure wrote something to the `.bin` file is written to this file in human-readable form. For long running tests, one can use `tail -f test/out/*.misc` to see the ends of all `.misc` files.
 * `.err` file: Everything that is written by `heavy_top` to standard output will be redirected to this file. (Does not apply to `make try`.)

## Analysis of the test results
The analysis of the test results is done in Matlab due to its easy way of manipulating data and visualizing it. All Matlab files are situated in `test/als/`. Here is a list with the most important functions that can and should be used:

 * `sol = load_latest_config_and_bin()`: Will look for the most recent duo of `.bin` and `.lua` files and will load them as a struct containing the most important variables from the lua file and the output from the binary file. The results are then available in the struct as `sol.rslt`, which is itself a struct with fields `t`, `q`, `v` and possibly `l` (for the Lagrange multipliers). Note that this will need `readLua`, see also the section Prerequisites in this readme. Note that you can also load the oldest file by passing the parameter `1`, the second to latest file by passing `-1` and so on.
 * `load_all_config_and_bin()`: Like above but will load all duos of `.bin` and `.lua` files in a cell of structs. You can pass a struct that can act as a pattern (eg. saying `pattern.steps = 10; load_all_config_and_bin(pattern)` will only load test results with exactly 10 steps.)
 * `solcell = calc_errors(solcell, refpattern)`: Calculates the errors with respect to a reference solution. The reference solution should be a part of the cell `solcell` and is determined by the pattern `refpattern`. There should be exactly one solution in `solcell` that matches `refpattern`. The relative and absolute errors are then available in the structs `solcell{k}.err.rel` and `solcell{k}.err.abs`.
 * `calc_SO3xR3_errors(solcell, refcell)`: Like above but the reference solution should be passed as a cell with exactly one entry. All errors are calculated by first converting the configuration to an element of SO(3)×ℝ³. The relative and absolute errors are then available in the structs `solcell{k}.SO3xR3_err.rel` and `solcell{k}.SO3xR3_err.abs`.
 * `delete_all_config_and_bin`, `delete_unfinished_config_and_bin`, `deleta_duplicate_config_and_bin`: These functions will delete all, unfinished, or duplicate duos of `.bin` and `.lua`. Use with caution.
 * `makexyplot(solcell, pattern, xname, yname, byname, varargin)`: Use this to plot stuff. For example `makexyplot(solcell, struct(), 'steps', 'err.abs.q', 'liegroup')` will open a figure where the absolute error in 𝑞 is plotted ove the number of steps. There will be several lines in the diagram distinguished by the value of `liegroup`. (The `varargin` is passed to `plot`.)
 * `matlab2csv(path, ax, nsteps)`: Will take the axes `ax` (if omitted the current axes `gca`) and convert the data to a csv-file that will be saved to the directory `path`. The parameter `nsteps` will take only export each `nsteps`th data point. (If omitted, `nsteps=1`. Choose higher `nsteps` for plots with a _lot_ of data points.) These csv files can be easily imported with `pgfplots` to produce beautiful plots in LaTeX.
 * `visualize_heavy_top(sol)`: Will produce a moving 3D plot of the solution. (Pass `1` as a second argument in order to calculate the position from the rotation instead of directly using the rotation.)
 * `sol = calc_energy(sol)`: Calculates the kinetic, potential and mechanical energies of the system available as `sol.rslt.energy_kinetic`, `sol.rslt.energy_potential` and `sol.rslt.energy`.

There are a few more functions (also in the subfolder `test\als\private`) that work in the background, but might be very important, eg. `load_config_and_bin` (that might be changed if `heavy_top` is changed) and `structmatch`.


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
 * `heavy_top.F90`: This Fortran source file contains the module `heavy_top` which defines a type `heavy_top_t` which extends the abstract problem type from the integrator. Here, all deferred procedures are implemented. See the documentation of the integrator for details.
 * `lie_group_functions.F90`: This Fortran source file contains functions related to Lie groups. Here, the functions that are implemented in the project `liegroup` are not implemented again. Rather, functions related to the Lie group of SO(3) and dual quaternions are implemented, which are not covered in the project `liegroup`. Note that since the integrators assume the configuration variable `q` to be a vector rather than a matrix, there are some functions that perform regular matrix operations but their inputs and/or outputs are matrices stored column-wise as a vector.
 * `get_line_of_variable_length.F90`: Defines a module with the same name that contains a function that can read a line of an opened file without knowing the length of thte line in advance.
 * `makefile`: A makefile that can be executed by calling `make`.
   1. Variables get defined. Modify them to your needs, especially the paths for the different projects and executable names, eg. how to call Matlab.
   2. Goals with prerequisites and how to make them. In this part, only values of the earlier defined variables are used. Usually, no modifications are needed here.

## Preprocessor functions and macros
In the source code, there are some preprocessor functions and macros used. Here is a list of most of them with some explanation:

 * `GL`: At the top, there is the macro `GL` defined as `#define GL(x) x`. The name is short for "glue" and can be used to "glue" the value of preprocessor variables to other text. As an example, let's say there is a preprocessor variable `VAR` and we want to use its value as part of a name. Simply using `VAR_this_text` will not work, because the preprocessor recognizes `VAR_this_text` as one entity. This can be circumvented by saying `GL(VAR)_this_text`. Now, the preprocessor knows that `GL(VAR)` is independent of the rest. If the value of `VAR` would be `example`, then `GL(VAR)_this_text` will become `example_this_text`, whereas `VAR_this_text` would stay `VAR_this_text`.
 * `INTEGRATOR`: This variable contains the name of the integrator (`gena`, `BLieDF`, `RATTLie` etc.)
 * `INT_RATTLie`, `INT_SHAKELie`, `INT_gena`, `INT_BLieDF`, `INT_varint4lie`: These get defined when RATTLie, SHAKELie, BLieDF, gena or varint4lie are used as an integrator. (SHAKELie is seldom used, since RATTLie is superior and varint4lie is unfinished and broken.)
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
