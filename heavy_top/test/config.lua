-- Integration interval
t0 = 0
te = 10
-- number of integration steps
steps = math.ceil((te-t0)*2^[[ 10 || 11 || 12 || 13 || 14 ]])
-- steps = math.ceil((te-t0)*2^10)
h = steps^-1 * (te-t0)
-- h  = [--[ 0.01 || 0.001 || 0.0001 || 0.00001 ]]
-- steps = math.ceil((te-t0)*h^-1)

-- problemset = [--[ 31 || 33 || 51 || 53 ||
--                 41 || 43 || 61 || 63 ]]
problemset = [[ 33 || 43 || 53 || 63 ]]
-- problemset = [--[ 10
--              || 20
--              || 31 || 32 || 33 || 34
--              || 41 || 42 || 43 || 44
--              || 51 || 52 || 53 || 54
--              || 61 || 62 || 63 || 64
--              || 71 || 72 || 73 || 74
--              ]]
-- problemset = [--[ 53 || 63 || 73 ]]
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

-- Use stabilized index-2 formulation (only applies to the constrained case)
if problemset > 0 then
   if problemset % 10 == 3 or problemset % 10 == 4 then
      stab2 = 1
      -- Baumgarte parameter
      -- my_C = [--[ 1. || 0.01 || 0.001 ]]
      my_C = 1.
      -- a_baumgarte = my_C / h
      a_baumgarte = [[ 1. || 10. || 100. ]]
      -- Use iterative projection
      stab_proj = 0
   else
      stab2 = 0
      -- Baumgarte parameter
      my_C = 0
      a_baumgarte = 0.
      -- Use iterative projection
      stab_proj = 0
   end
else
   -- this may be changed when problemset == -1
   stab2 = 0
   -- Baumgarte parameter
   my_C = 0
   a_baumgarte = 0.
   -- Use iterative projection
   stab_proj = 0
end

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Problem options   -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- Number of heavy tops
number_heavy_tops = 1

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
output_t_at = 0
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

-- Algorithmic paramters for half-explicit Lie Runge-Kutta

-- A = { 0.,  0.,
--       1.,  0.,
--      1/2, 1/2}
-- c = {0., 1., 1.}
-- d = {1/2, 1/2}
-- order  = 2
-- stages = 2
-- stages_bar = 2
-- local_error_control = false
-- step_size_control = false

-- A = { 0.,  0.,  0.,
--      1/2,  0.,  0.,
--      2/9, 4/9,  0.,
--      1/4,  0., 3/4}
-- c = {0., 1/2, 2/3, 1.}
-- d = {0., -2., 3.}
-- order  = 3
-- stages = 3
-- stages_bar = 3
-- local_error_control = false
-- step_size_control = false

-- -- Only for constrained problems
-- A = { 0.,   0.,  0.,  0., 0.,
--       1.,   0.,  0.,  0., 0.,
--      3/8,  1/8,  0.,  0., 0.,
--     -1/2, -1/2,  2.,  0., 0.,
--      1/6,  0., 2/3, 1/6, 0.,
--      1/6 - math.sqrt(3)/108,  0., 1/3 - 4*math.sqrt(3)/27,  -7*math.sqrt(3)/108, math.sqrt(3)/18}
-- c = {0., 1., 1/2, 1., 1., 0.2113}
-- d = {0., 0., 0., 0., 1.}
-- order  = 4
-- stages = 4
-- stages_bar = 5
-- local_error_control = false
-- step_size_control = false

-- -- Only for unconstrained problems
-- A = { 0.,  0.,  0., 0.,
--      1/2,  0.,  0., 0.,
--       0., 1/2,  0., 0.,
--       0.,  0.,  1., 0.,
--      1/6, 1/3, 1/3, 1/6}
-- c = {0., 1/2, 1/2, 1., 1.}
-- d = {0., 0., 0., 1.}
-- order  = 4
-- stages = 4
-- stages_bar = 4
-- local_error_control = false
-- step_size_control = false

-- -- Automatic step size control DOPRI5 (kirchhoff = 0 --> unconstrained)
-- A = {   0.,          0.,         0.,         0.,          0.,    0., 0.,
--        1/5,          0.,         0.,         0.,          0.,    0., 0.,
--       3/40,        9/40,         0.,         0.,          0.,    0., 0.,
--      44/45,      -56/15,       32/9,         0.,          0.,    0., 0.,
-- 19372/6561, -25360/2187, 64448/6561,   -212/729,          0.,    0., 0.,
--  9017/3168,     -355/33, 46732/5247,     49/176, -5103/18656,    0., 0.,
--     35/384,          0.,   500/1113,    125/192,  -2187/6784, 11/84, 0.,
--     35/384,          0.,   500/1113,    125/192,  -2187/6784, 11/84, 0.}
-- c = {0., 1/5, 3/10, 4/5, 8/9, 1., 1., 1.}
-- b = {5179/57600, 0., 7571/16695, 393/640, -92097/339200, 187/2100, 1/40}
-- d = {0., 0., 0., 0., 0., 0., 1.}
-- order  = 5
-- order_step_control = 4
-- stages = 6
-- stages_bar = 7
-- local_error_control = false
-- step_size_control = false

-- Automatic step size control HEDOP5 (kirchhoff = 1 --> constrained)
A = {   0.,          0.,         0.,         0.,          0.,    0., 0.,
       1/5,          0.,         0.,         0.,          0.,    0., 0.,
      3/40,        9/40,         0.,         0.,          0.,    0., 0.,
     44/45,      -56/15,       32/9,         0.,          0.,    0., 0.,
19372/6561, -25360/2187, 64448/6561,   -212/729,          0.,    0., 0.,
 9017/3168,     -355/33, 46732/5247,     49/176, -5103/18656,    0., 0.,
    35/384,          0.,   500/1113,    125/192,  -2187/6784, 11/84, 0.,
    -18611506045861/19738176307200, 59332529/14479296,   -2509441598627/893904224850,    2763523204159/3289696051200,  -41262869588913/116235927142400, 46310205821/287848404480, -3280/75413}
c = {0., 1/5, 3/10, 4/5, 8/9, 1., 1., 19/20}
b = {5179/57600, 0., 7571/16695, 393/640, -92097/339200, 187/2100, 1/40}
d = {0., 0., 0., 0., 0., 0., 1.}
order  = 5
order_step_control = 4
stages = 6
stages_bar = 7
local_error_control = false
step_size_control = false

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



