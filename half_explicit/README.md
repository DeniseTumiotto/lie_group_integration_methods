# Half-explicit Runge-Kutta methods `half_explicit` for holonomic constrained systems on Lie group structured configuration spaces

These are a time integration methods that can be used to solve differential algebraic equations of the form

        ̇𝑞(𝑡) = d𝐿_𝑞(𝑒) ̃𝑣(𝑡)
    𝐌 ⋅ ̇𝑣(𝑡) = -𝑔(𝑡,𝑞,𝑣) - 𝐁ᵀ(𝑞(𝑡))⋅λ(𝑡)
           𝟎 = Φ(𝑞(𝑡))

on a Lie group 𝐺 of dimension 𝑛. Here 𝑞(𝑡)∈𝐺 is the configuration of the system and 𝑣(𝑡)∈ℝⁿ the associated velocity vectors. 𝐿 is the left translation and d𝐿 its derivative. We have the tilde operator which maps 𝑛-vectors to the Lie algebra of 𝐺. The 𝑛×𝑛-matrix 𝐌 is a constant mass matrix, 𝑔 is a function that gives the negative of all generalized forces including Coriolis (or inertial) forces. Moreover, Φ is the constraint function with Φ(𝑞)∈ℝᵏ and 𝐁 is its derivative in the sense
    𝐁(𝑞)⋅𝑤 = dΦ(𝑞) d𝐿_𝑞(𝑒) ̃𝑤  for all  𝑤∈ℝⁿ.
Using the derivative operator 𝐃, we could also write in short 𝐁=𝐃Φ. Furthermore, λ(𝑡)∈ℝᵏ are the Lagrange multipliers.

Alternatively, it can also solve unconstrained equations of motion of the form

        ̇𝑞(𝑡) = d𝐿_𝑞(𝑒) ̃𝑣(𝑡)
    𝐌 ⋅ ̇𝑣(𝑡) = -𝑔(𝑡,𝑞,𝑣)

the half-explicit Runge-Kutta methods `half_explicit` were implemented in modern Fortran for Linux. (Compiling it on Windows, however should not be too hard.)

## Content
`half_explicit` is a Fortran module that defines the following types:

 * The abstract type `half_explicit_problem`. When implementing a problem, this type needs to be extended to a non-abstract type by implementing all deferred procedures.
 * The type `half_explicit_options`, which contains parameters of the integration scheme.
 * The type `half_explicit_statistics`, which will contain information about the runtime after a performed integration step.

In the following, `prob` or `this` will be the variable of class `half_explicit_problem`.

### Deferred procedures to be implemented
Here is a list of all deferred procedures that must be implemented when extending `half_explicit_problem` including a short description. Note that it is helpful if a procedure that is never called or a dummy, it should contain `error stop "<this procedure> should never be called"`.

 * `half_explicit_M`: This function should return the mass matrix 𝐌 in dependence of the configuration 𝑞. Note however, that mass matrices that are not constant have _not_ been tested and the results might not be correct. This function will never be called, if `prob%opts%diag_mass_matrix==1`.
  * `half_explicit_diag_M`: This function should return the diagonal of the mass matrix 𝐌 as a vector. (See also `half_explicit_M`) This function will never be called, if `prob%opts%diag_mass_matrix==0`.
 * `half_explicit_g`: This function should return 𝑔(𝑡,𝑞,𝑣), the negative of of all applied generalized forces. Note that the result must include any Coriolis forces (inertial forces) but exclude, of course, any constraint forces.
 * `half_explicit_qlpexpDqtilde`: The cryptic name means "Q Lie Product EXP of Delta Q TILDE". It should return the result of 𝑞∘exp(ℎ⋅tilde(Δ𝑞)), where ∘ is the Lie group product, exp:𝑇ₑ𝐺→𝐺 the exponential map and tilde:ℝⁿ→𝑇ₑ𝐺 the tilde operator. We could also write this as 𝑞∘expt(ℎ⋅Δ𝑞).
 * `half_explicit_itlbtvtw`: The cryptic name means "Inverse Tilde of Lie Bracket of Tilde V and Tilde W". This is a very complicated way of saying that this function should return ̂𝑣⋅𝑤, where ̂𝑣 is the application of the hat operator to the vector 𝑣. This function will never be called if the system is unconstrained, ie. `prob%opts%constrained==0`.
 * `half_explicit_tilde`: This function exists only for backwards compatibility. It may be a dummy function and may/can/will(?) be removed in the future. (The reason being that nowhere in the algorithm the actual 𝑇ₑ𝐺 is used.)
 * `half_explicit_Ct`: This function should return the tangent damping matrix 𝐂, which is the derivative wrt 𝑣 of the residual 𝐌⋅̇𝑣 + 𝑔(𝑡,𝑞,𝑣) + 𝐁ᵀ(𝑞)⋅λ. This means that 𝐂=d/d𝑣 𝑔(𝑡,𝑞,𝑣). This function will never be called if `prob%opts%use_num_Ct==1` or `prob%opts%use_no_Ct==1`.
 * `half_explicit_Kt`: This function should return the tangent stiffness matrix 𝐊 in the unconstrained case (ie. `prob%opts%constrained==0`). That means 𝐊 is the derivative wrt 𝑞 of the residual 𝐌⋅̇𝑣 + 𝑔(𝑡,𝑞,𝑣), so we have 𝐊=𝐃_𝑞(𝐌⋅̇𝑣 + 𝑔(𝑡,𝑞,𝑣)) or in other words: For all vectors 𝑤 it holds 𝐊⋅𝒘=d_𝑞(𝐌⋅̇𝑣 + 𝑔(𝑡,𝑞,𝑣)) d𝐿_𝑞(𝑒) ̃𝑤. Note that if the mass matrix 𝐌 is independent of 𝑞, 𝐊 is just the derivative of 𝑔 wrt. 𝑞. This function will never be called if `prob%opts%use_num_Ct==1` or `prob%opts%use_no_Ct==1` and also if `prob%opts%constrained==1`.
 * `half_explicit_Kt_lambda`: This function should return the tangent stiffness matrix 𝐊 in the constrained case (ie. `prob%opts%constrained==1`). That means 𝐊 is the derivative wrt 𝑞 of the residual 𝐌⋅̇𝑣 + 𝑔(𝑡,𝑞,𝑣) + 𝐁ᵀ(𝑞)⋅λ, so we have 𝐊=𝐃_𝑞(𝐌⋅̇𝑣 + 𝑔(𝑡,𝑞,𝑣) + 𝐁ᵀ(𝑞)⋅λ) or in other words: For all vectors 𝑤 it holds 𝐊⋅𝒘=d_𝑞(𝐌⋅̇𝑣 + 𝑔(𝑡,𝑞,𝑣) + 𝐁ᵀ(𝑞)⋅λ) d𝐿_𝑞(𝑒) ̃𝑤. Note that if the mass matrix 𝐌 might be independent of 𝑞. This function will never be called if `prob%opts%use_num_Ct==1` or `prob%opts%use_no_Ct==1` and also if `prob%opts%constrained==0`.
 * `half_explicit_Tg`: This functions should return the tangent operator 𝐓(ℎ⋅Δ𝑞). The tangent operator should be implemented in such a way that it works reliably, even if the norm of ℎ⋅Δ𝑞 is very small.
 * `half_explicit_d_Tg_inv_T`: This function should return the Jacobi matrix of 𝐓⁻ᵀ(𝑣)⋅𝑤 with respect to 𝑣.
 * `half_explicit_norm`: This function may be a dummy function. It is only here for backwards compatibility. It may/should/will(?) be removed in the future.
 * `half_explicit_outputFunction`: Despite its name, this is not a function but rather a subroutine. It will be called after each successful integration step. Note, that all arguments are `intent(in)`, meaning that the problem object may not be altered. There are currently three possible values of `info`:
   * `info==0`: "Initialization": The `half_explicit_outputFunction` is called once after the problem was initialized (ie. after `half_explicit_init` was called and, in the constrained case, initial values for ̇𝑣, 𝑎 and λ were calculated.) Note that here, we should not open files, because there is no way of storing the identifier. Open files in `half_explicit_init` or before calling the integration routine.
   * `info==1`: "Normal output": The `half_explicit_outputFunction` was called after a successful integration step.
   * `info==99`: "Termination": The `half_explicit_outputFunction` is called once at the very end of the full integration run.
 * `half_explicit_init`: This subroutine should initialize the problem. There are a few things that must happen in this subroutine:
   * Set the sizes `this%sizeq`, `this%sizev` and in the constrained case also `this%sizel`.
   * Allocate `this%q`, `this%v`, `this%vd`, `this%a` and in the constrained case also `this%l`. If we use the stabilized index-2 formulation, `this%eta` has to be allocated as well. Not that it is a good idea to check, whether the variables are already allocated and if so, deallocate them before allocating.
   * Set initial values:
     * `this%t = this%opts%t0`
     * `this%q`
     * `this%v`
   * Allocate and fill `this%opts%jour` if `this%opts%banded_iteration_matrix==1`.
 * `half_explicit_Phi`: This function should return the constraint function Φ(𝑞). It will never be called in the unconstrained case (ie. `prob%opts%constrained==0`).
 * `half_explicit_B`: This function should return the derivative of the constraint function 𝐁(𝑞)=𝐃Φ(𝑞). In other words we have 𝐁(𝑞)⋅𝑤 = d/d𝑞 Φ(𝑞) d𝐿_𝑞(𝑒) ̃𝑤 for all vectors 𝑤. This function will never be called in the unconstrained case (ie. `prob%opts%constrained==0`).
 * `half_explicit_Z`: This function should return the curvature term 𝐃_𝑞(𝐁(𝑞)⋅𝑣)⋅𝑣=d_𝑞(𝐁(𝑞)⋅𝑣) d𝐿_𝑞(𝑒) ̃𝑣. This function will never be called in the unconstrained case (ie. `prob%opts%constrained==0`) and only is used in the calculation of the initial values for ̇𝑣, 𝑎 and λ.
 * `half_explicit_matZ`: This function should return a matrix 𝐙=𝐃_𝑞(𝐁(𝑞)⋅𝒗)⋅𝐓 or in other words, it holds 𝐙⋅𝑤=d_𝑞(𝐁(𝑞)⋅𝒗) d𝐿_𝑞(𝑒) tilde(𝐓⋅𝑤) for all vectors 𝑤. Note that here, 𝐓 is a matrix argument. This function will never be called in the unconstrained case (ie. `prob%opts%constrained==0`) and only is used in the stabilized index-2 case (ie. `prob%opts%stab2==1`). Note that if the preprocessor macro `NOZ` is set, the occurence of 𝐙 is omitted and thus this function will not be called.

### Public subroutines
The following subroutines are meant to be called by the problem file.

 * `half_explicit_integrate`: Call this subroutine in order to start the integration process after all important variables and integrator options have been set. This routine will call `half_explicit_init` once and `half_explicit_outputFunction` after each successful integration step. Note that `half_explicit` does not take care of saving intermediate results. This has to be done by the `half_explicit_outputFunction`.
 * `half_explicit_print_stats`: Prints the contents of `prob%half_explicit_stats` to the standard output.
 * `half_explicit_cleanup`: Resets most of the internal variables and integrator options and deallocates most internal allocatable variables. This will, of course, not reset any variables that are added to `half_explicit_problem` by its extension.

### Integrator options
Most of the integrator options are found in `prob%opts`, except of the Runge-Kutta coefficients, which are directly part of the problem object for compatibility reasons.
Here is a list of all other integrator options in `prob%opts`:

 * `constrained`: Set this to `1` if the system is constrained, for unconstrained systems set this to `0`.
 * `stab2`: Only applies in the constrained case. Set this to `1`, if the stabilized index-2 formulation should be used. Set this to `0` in order to use the index-3 formulation. Note that the stabilized index-2 formulation is usually the better choice.
 * `a_baumgarte`: Only applies in the constrained index-2 stabilized case. Set this `real(8)` variable to the value of the Baumgarte coefficient.
 * `stab_proj`: Only applies in the constrained index-2 stabilized case. Set this to `1` if the stabilization is performed via projection of the solution on the submanifold generated by the constraints.
 * `const_mass_matrix`: Set this to `1` if the mass matrix 𝐌 does not depend on the configuration 𝑞. Set this to `0` if 𝐌 depends on 𝑞. Note that `half_explicit` was not tested with nonconstant mass matrices, results may not be accurate.
 * `diag_mass_matrix`: Set this to `1` if the mass matrix 𝐌 is a diagonal matrix. In this case `half_explicit_diag_M` will be used instead of `half_explicit_M`. Set this to `0` for full mass matrices.
 * `banded_iteration_matrix`: Set this to `1` if the iteration matrix has band structure or can be rearranged to a matrix with band structure. This is usually the case when integrating a system with several bodies that are chained but otherwise don't interact directly. Set this to `0` if the iteration matrix can be full.
 * `nr_subdiag`: Only applies to the case with banded iteration matrix. Set this to the number of subdiagonals (diagonals below the main diagonal) of the band structure not counting the main diagonal. Note that the Newton method may fail or still converge if `nr_subdiag` is smaller than the actual numer of subdiagonals. If chosen too large, the integration will be slower than it has to be.
 * `nr_superdiag`: See `nr_subdiag`, but for the superdiagonals (diagonals above the main diagonal) instead of the subdiagonals.
 * `jour`: Only applies to the case with banded iteration matrix. This should be an integer vector such that `St(jour,jour)` is a banded matrix, where `St` is the iteration matrix. The size should be:
   * Unconstrained case: `prob%sizev`
   * Constrained index-3 case: `prob%sizev+prob%sizel`
   * Constrained stabilized index-2 case: `prob%sizev + 2*prob%sizel`
 * `pertube`: Only applies in the constrained index-3 case. Set this to `1` if the initial values should be perturbed in order to minimize spurious oscillations in the first integration steps. Set this to `0` if the initial values should not be perturbed. Note that these spurious oscillations do generally not appear in the stabilized index-2 case. Also, I'm sorry, but for compatibility reasons, this variable is misspelled.
 * `use_num_Ct`: Set this to `1` if the tangent damping matrix should be calculated via finite differences. This is often slower and less accurate, but then `prob%half_explicit_Ct` is never called. Set this to `0` if the damping matrix should be calculated by `prob%half_explicit_Ct`.
 * `use_num_Kt`: Set this to `1` if the tangent stiffness matrix should be calculated via finite differences. This is often slower and less accurate, but then `prob%half_explicit_Kt` and `prob%half_explicit_Kt_lambda` is never called. Set this to `0` if the damping matrix should be calculated by `prob%half_explicit_Kt` or `prob%half_explicit_Kt_lambda`.
 * `no_Ct`: Set this to `1` if the tangent damping matrix should be completely omitted in the iteration matrix. Then, `half_explicit_Ct` is never called. This speeds up the integration a lot, but may lead to problems with convergence of the Newton method. Set this to `0` if the damping matrix should be included in the iteration matrix.
 * `no_Kt`: Set this to `1` if the tangent stiffness matrix should be completely omitted in the iteration matrix. Then, `half_explicit_Kt` and `half_explicit_Kt_lambda` are never called. This speeds up the integration a lot, but may lead to problems with convergence of the Newton method. Set this to `0` if the stiffness matrix should be included in the iteration matrix.
 * `atol`: Set this `real(8)` variable to the absolute tolerance to be used in the Newton method.
 * `rtol`: Set this `real(8)` variable to the relative tolerance to be used in the Newton method.
 * `imax`: Set this variable to the maximum number of iterations after which the Newton iteration is considered to not converge. If this variable is to low, integration might not suceed.
 * `t0`: Set this `real(8)` variable to the beginning of the time integration interval 𝑡₀. Don't forget to set `this%t = this%t0` in `half_explicit_init`.
 * `te`: Set this `real(8)` variable to the end of the time integration interval 𝑡ₑ.
 * `nsteps`: Set this variable to the number of integration steps (of equal length) to be made between 𝑡₀ and 𝑡ₑ. This means the step size can be calculated by ℎ=(𝑡ₑ-𝑡₀)/`nsteps`.
 * `step_size_control`: Set this to `.true.` if the time step size should be controlled by the estimate of the local error.
 * `local_error_control`: Set this to `.true.` if the local error estimate should be evaluated and saved in the output file.
 * `update_a`: Set this to `.true.` if the Baumgarte coefficient needs to be updated together with the new time step size.
 * `facmax`: Set this `real(8)` variable to the maximum factor in the evaluation of the new time step size, usually in the interval `1.5-5.0`.
 * `facmin`: Set this `real(8)` variable to the minimum factor in the evaluation of the new time step size, default set to `0.1`.
 * `fac`: Set this `real(8)` variable to the factor in the evaluation of the new time step size, usually `0.8, 0.9` or in dependance of the order of the numerical method.
 
### Compiler flags
This project has been written for `gfortran` on Linux, although porting to code to different compilers or different platforms should not be too hard.
The makefile of `half_explicit` will look for the variable `EXTRAFFLAGS` that is supposed to contain compiler flags. In order to make the variable `EXTRAFFLAGS` visible to the makefile, you should put `export EXTRAFFLAGS` in the problem makefile after defining it.
Here is a list of helpful compiler flags:

 * `-O`: Pass this flag to "optimize" the code. It will produce code that may run (a lot) faster. The optimization can be controlled in levels: `-O0` to turn it off completely, `-O1` and `-O3` exist, see `gfortran`s manual, and `-O2` is equivalent to `-O`.
 * `-Wall`: Turn on "all" warnings while compiling. Helpful for debugging.
 * `-Wextra`: Turn on even more warnings than "all". Helpful for debugging.
 * `-Dpure=''`: This defines the preprocessor variable `pure` to be the empty string. This will cause the Fortran keyword `pure` to be removed from the whole code, making all otherwise `pure` procedures non-`pure`. This is very helpful for debugging, because in a `pure` procedure, no side effects such as printing are allowed. (Exept right before `error stop`). _Note that this is black magic and probably everybody will tell you not to do such things in a program that should produce sensible output._
 * `-g`: Turn on debug mode. Extremely helpful for debugging, obviously. Slows down the program on the other hand.
 * `-fbounds-check`: Check bounds of vectors and such. May help to discover errors.
 * `-fimplicit-none`: Automatically uses `implicit none` everywhere. Use only for debugging and be sure to write `implicit none` everywhere it belongs.
 * `-fbacktrace`: Backtrace. Probably useful for debugging.
 * `-fcheck=all`: Check stuff. Probably useful for debugging.
 * `-finit-real=NaN`: Initialize every `real` variable with `NaN` (not a number). Useful for finding variables that are used before they were defined. (Usually undefined variables have a random value that happens to be in the memory before.)
 * `-finit-integer=-77`: See above. If this flag is given you find an integer with value `-77` you probably forgot to define it.
 * `-ffpe-trap=zero,invalid`: Find the most common floating point errors (fpe) such as: dividing by zero and other invalid stuff. (Usually would result in "undefined behaviour", whatever that means.)

It usually makes sense to define `EXTRAFFLAGS = -O` for producing test results and defining `EXTRAFFLAGS` to be all other mentioned flags for debugging.

### Preprocessor variables
There are a few preprocessor variables that can be defined when compiling `half_explicit`. Some of which are purely for debugging purposes, but some might be really useful. A preprocessor variable `MYVAR` can be defined by giving the compiler (`gfortran`) the compiler flag `-DMYVAR`. Sometimes it is meaningful to give the preprocessor variable a value, say `0`. This can be done by passing the compiler glag `-DMYVAR=0`.
Here is a list of some preprocessor variables:

 * `NOZ`: Only applies to the constrained index-2 case. If this preprocessor variable is defined, then the function `half_explicit_matZ` will never be called and the curvature term is omitted in the iteration matrix. Note that this preprocessor variable is defined in the makefile of `half_explicit`. If you actually want `half_explicit_matZ` to be used you have to edit the makefile.
 * `ZEROINIT`: If this preprocessor variable is defined, the Newton method is initialized with zero instead of with the values from the last integration step. This usually makes the integrator take longer and often, it will not converge.
 * `STNUM`: If this preprocessor variable is defined, the iteration matrix will be determined by finite differences completely. This is usually a lot slower and less reliable.
 * `DEBUG_PRINT_ITERATION_MATRIX_AT`: This variable should have a valid numeric value if it is defined. In that case, integration will be stopped at the first integration step that is equal or larger than the value and the iteration matrix is printed to standard output.
 * `ITERATIVE_RELAXATION`: If this variable is defined a iterative relaxation approach to the Newton method will be used. This will make the intergrator a little slower but more reliable. _Note that this has never been successfully tested and this functionality is probably completely broken._
 * All other variables are purely for special debugging tasks or their meaning is not clear to me any more.

All of the aforementioned preprocessor variables, except for `NOZ`, unlock features meant for debugging or that are unfinished.


## Usage
In order to implement a problem that should be integrated with `half_explicit` the following files are probably needed to be created:
 
 * `problem.F90` which contains a module `problem` which used the module `half_explicit` and extends the abstract type `half_explicit_problem` to a non-abstract type.
 * `main.F90` which contains the `program main` and uses the module `problem` implemented in `problem.F90`.
 * `makefile` which contains recipes to compile and maybe run the program.

The contents of these files could look like this:

### `main.F90`
```fortran
! implement the program "main" which is the program that will be executed, when the compiled executable is called
program main
   ! use the module "problem" implemented in "problem.F90"
   use problem

   ! this goes without saying
   implicit none

   ! define variables
   type(problem_t)  :: prob

   ! we have a variable "prob" of type "problem_t"
   ! first, we will set all important variables, that were added to "problem_t" in the
   ! extension of "half_explicit_problem":
   prob%some_really_important_variable = 5.3e-2_8
   prob%might_be_multidimensional = [0.0_8, 1.0_8, 2.0_8, 3.0_8, 4.0_8]
   prob%or_just_an_int = 4

   ! Now we set all integrator options like:
   prob%opts%t0 = 0.0_8
   prob%opts%te = 1.0_8
   prob%opts%const_mass_matrix = 1
   !etc etc

   ! Here you would maybe open some files that the output function of the integrator
   ! can write to

   ! Start the actual integration
   call prob%half_explicit_integrate()

   ! Here you would maybe close all files that were opened

   ! Maybe you want to show some statistics
   call prob%half_explicit_print_stats()

   ! Clean up after yourself! This is very important when you want to call
   ! "half_explicit_integrate" twice on the same problem variable
   call prob%half_explicit_cleanup()

   ! You are done
   print *, "Done"
end program main
```

### `problem.F90`
```fortran
! define the module "problem"
module problem
   ! we need to use the integrator module "half_explicit"
   use half_explicit

   ! this goes without saying in Fortran
   implicit none

   ! extend the abstract type "half_explicit_problem"
   type, extends(half_explicit_problem) :: problem_t
      ! here come some internal problem-specific variables
      ! <type of variable> :: <name of variable>
      real(8)   :: some_really_important_variable
      real(8)   :: might_be_multidimensional(5)
      integer   :: or_just_an_int
      ! etc etc
   contains
      ! referencing the former deferred procedures
      ! procedure :: <name of procedure in the abstract interface of half_explicit>  =>  <name of the actual procedure in this module>
      procedure :: half_explicit_M       => problem_M
      ! etc etc
      ! of course, we need to reference ALL deferred procedures here
   end type problem_t

contains
   ! now the actual procedures have to be implemented as an example, here is the mass matrix:
   pure function problem_M(this, q) result(rslt)
      ! input
      class(problem_t), intent(in)  :: this
      real(8),          intent(in)  :: q(:)
      ! output
      real(8)                       :: rslt(this%sizev, this%sizev)
      !
      ! here should be some fancy code:
      rslt = 0.0_8
      ! or here could also be an error if this function should never be called 
      ! (eg in this case because only the diagonal mass matrix should be used)
      error stop "Full mass matrix not supported, use diagonal mass matrix instead"
   end function problem_M

   ! now all other procedures have to be implemented. As a template, use the definition
   ! of the procedures in the abstract interface of the "half_explicit" module. But keep in mind
   ! that the "import" statement needs to be removed and all occurences of "half_explicit_problem"
   ! have to be replaced by "problem_t".
end module problem
```
           
### `makefile`
```makefile
# This is a makefile, so make sure all indentations are tabs not spaces!

# Here we could define some variables, eg. the Fortran compiler to be used
FC = gfortran
# Variables can be referenced by saying eg $(FC)
# To keep this more readable I will not use variables, although it would be
# the right way to do this!

# Define the variable EXTRAFFLAGS and export it. 
# FFLAGS stands for Fortran FLAGS and is supposed to contain compiler flags.
# The makefile of half_explicit will look for the exported variable EXTRAFFLAGS and
# use the contents to build itself.
EXTRAFFLAGS = -O
export EXTRAFFLAGS

# Tell make which of the recipies are not actual files
.PHONY: default clean

# Here comes the recipies
# <name of the recipie>: <ingredient1> <ingredient2> <etc>
# The first recipe is the one that will be made, when make is called without arguments
default: problem makefile
    # Do nothing. This will still cause all ingredients to be made!
    # The ingredient "makefile" is there so everything will be recompiled
    # if the makefile is changed

# First recipe that is an actual file name
problem: half_explicit.o problem.o main.o
    # the $@ is a variable that will be the name of the recipe
    # the $+ gives all prerequisites (with duplicates and in order)
    gfortran -o $@ $+ -llapack
    # the option "-o" means "output"
    # in this case, what is called will be
    #gfortran -o problem half_explicit.o problem.o main.o -llapack
    # which will output the executable file "problem" by using the object files
    # half_explicit.o, prblem.o and main.o
    # The "-llapack" means that it should link the library LAPACK, which is
    # needed for solcing systems of linear equations in half_explicit

# Build the object files
problem.o: problem.F90 half_explicit.o
    # the $< is a variable that will be the name of the first ingredient
    # the other ingredients are only there, because they already need to be done
    gfortran -c -o $@ $<
    # the option "-c" means "only compile, no linking"

half_explicit.o:
   # Call this version of make by using $(MAKE)
   $(MAKE) -C ../half_explicit/
   # the -C flag means "go to the following directory"
   # this assumes, that the "half_explicit" project is located at "../half_explicit/"
   #
   # Then copy the important files half_explicit.o (object file) and half_explicit.mod (module file)
   # in the current directory (which is ".")
   cp ../half_explicit/obj/half_explicit.o .
   cp ../half_explicit/obj/half_explicit.mod .

main.o: main.F90
   gfortran -c -o $@ $<

clean:
   # Remove all generated files
   # The "-" before the "rm" means "ignore if there is an error executing"
   # (An error would occur, eg. when there are no files to remove)
   -rm *.o
   -rm *.mod
   -rm problem
```

If everything is done right a call of `make` will build the executable `problem` which then can be executed by using `./problem`.

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
