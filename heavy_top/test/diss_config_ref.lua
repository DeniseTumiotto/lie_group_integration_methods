-- In order to make this a little easier to control

-- Integration interval
t0 = 0
te = 1
-- number of integration steps
--steps = math.ceil((te-t0)*2^20)

problemset = 33
--problemset = [--[ 31 || 32 || 33 || 34 ]]
--problemset = [--[ 10
--             || 20
--             || 31 || 32 || 33 || 34
--             || 41 || 42 || 43 || 44
--             || 51 || 52 || 53 || 54
--             || 61 || 62 || 63 || 64
--             || 71 || 72 || 73 || 74
--             ]]
--problemset = [--[ 53 || 63 || 73 ]]
-- For no predefined problem set use -1
-- Lie groups:
--  10: SO(3)
--  20: S^3
--  30: SO(3) x R^3
--  40: S^3 x R^3
--  50: SE(3)
--  60: S^3 |x R^3 (semidirect product)
--  70: Unit dual quaternions
-- Formulations and options:
--  +0: unconstrained (only applies to SO(3) and S^3)
--  +1: index-3 formulation, no perturbation
--  +2: index-3 formulation, with perturbation
--  +3: stabilized index-2 formulation, using no_K and no_D
--  +4: stabilized index-2 formulation


-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Problem options   -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- Mass of the top
mass = 15

-- Gravity including direction
gravity = {0, 0, -9.81}

-- Diagonal elements of inertial tensor wrt. the center of mass
inerJ = {0.234375 , 0.46875 , 0.234375}

-- -- -- Initial values (will be adapted to the chosen Lie group formulation -- -- --
-- Initial positions (in R^3)
x0 = {0, 1, 0}
-- Initial rotation (in S^3)
p0 = {1, 0, 0, 0}
-- Initial angular velocity
Om0 = {0 , 150 , -4.61538}
-- Initial velocity is chosen to be consistent with Om0:
-- u0 = R(p0)*cross(Om0,X)
-- X is the reference point calculated by
-- X = kpvp(p0,x0) = R^T(p0)*X

-- Lie group formulation
--  1: SO(3)
--  2: S^3
--  3: SO(3) x R^3
--  4: S^3 x R^3
--  5: SE(3)
--  6: S^3 |x R^3 (semidirect product)
--  7: Unit dual quaternions
if problemset > 0 then
   liegroup = math.floor(problemset/10)
else
   -- this may be changed when problemset == -1
   liegroup = 1
end

-- -- -- Output options -- -- --
-- How to output the calculated data:
--  0: No output
--  1: Output in terms of the Lie group formulation
--  2: Output in terms of SO(3) x R^3
output_type = 1
-- Output only at certain times
output_t_at = 1
t_output_at_multiples_of = 1/2^8

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Integrator options   -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- Algorithmic paramters for gena (for rho_infinity = 0.9)
alpha_m =   8/ 19
alpha_f =   9/ 19
beta    = 100/361
gamma   =  21/ 38

-- Algorithmic paramters for BLieDF
k_bdf = 2

-- Use constant mass matrix
const_mass_matrix = 1
-- Use diagonal mass matrix
diag_mass_matrix = 0
-- Use banded solvers for the iteration matrix
banded_iteration_matrix = 0
-- Recalculate the iteration matrix in ever Newton step
recalc_iteration_matrix = 0
-- Perturb initial values (only applies to the constrained case with direct integration of the index-3 formulation)
if problemset > 0 then
   if problemset % 10 == 2 then
      perturb = 1
   else
      perturb = 0
   end
else
   -- this may be changed when problemset == -1
   perturb = 1
end
perturb_s = 1.0
-- Use numerical approximation for stiffness matrix K
use_num_K = 0
-- Use numerical approximation for damping matrix D
use_num_D = 0
-- Omit stiffness matrix K and damping matrix D in the iteration matrix
if problemset > 0 then
   if problemset % 10 == 3 then
      no_K = 1
      no_D = 1
   else
      no_K = 0
      no_D = 0
   end
else
   -- this may be changedn when problemset == -1
   no_K = 0
   no_D = 0
end

-- Relative error bound for the Newton-Raphson method
rtol = [[((tol)) 1.0e-8 ]]
-- Absolute error bound for the Newton-Raphson method
atol = [[((tol)) 1.0e-10 ]]
-- Maximum unsuccessful iteration steps after which the method is considered not to converge
imax = 100

-- Integration interval and step size (see top of the file)

-- Use stabilized index-2 formulation (only applies to the constrained case)
if problemset > 0 then
   if problemset % 10 == 3 or problemset % 10 == 4 then
      stab2 = 1
   else
      stab2 = 0
   end
else
   -- this may be changed when problemset == -1
   stab2 = 0
end
