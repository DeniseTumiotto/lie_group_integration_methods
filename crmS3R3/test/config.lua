math.randomseed(26031990)
function randgrid(start, stop, length)
   maxh = (stop - start)/length
   minh = maxh/4
   out = {}
   out[1] = start
   h = 3/4*maxh
   smaller = false
   bigger = false
   outind = 1
   while true do
      repeat
         if smaller then
            newh = h*(math.random()*(1.0 - 0.5) + 0.5)
         elseif bigger then
            newh = h*(math.random()*(1.5 - 1.0) + 1.0)
         else
            newh = h*(math.random()*(1.5 - 0.5) + 0.5)
         end
         if newh > maxh then
            smaller = true
            break -- = continue
         end
         if newh < minh then
            bigger = true
            break -- = continue
         end
         if out[outind] + 1.5*newh > stop then
            out[outind+1] = out[outind] + (stop - out[outind])/2
            out[outind+2] = stop
            out[outind+3] = stop + (stop - out[outind])/2
            return out
         end
         out[outind+1] = out[outind] + newh
         outind = outind + 1
         h = newh
         smaller = false
         bigger = false
         break
      until true
   end
end

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Integrator options   -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- Algorithmic paramters for gena(for rho_infinity = 0.9)
rho_infinity = 0.9
alpha_m = (2*rho_infinity - 1)/(rho_infinity + 1)
alpha_f = rho_infinity/(rho_infinity + 1)
beta    = 1/4 * (1 - (alpha_m - alpha_f))^2
gamma   = 1/2 - (alpha_m - alpha_f)
--alpha_m =   8/ 19
--alpha_f =   9/ 19
--beta    = 100/361
--gamma   =  21/ 38

--Algorithmic parameters for BLieDF
k_bdf = 2

-- Algorithmic paramters for half-explicit Lie Runge-Kutta

-- [--[
A = { 0.,  0.,
     1/2,  0.,
      0.,  1.}
c = {0., 1/2, 1.}
d = {0.,  1.}
order  = 2
stages = 2
stages_bar = 2
local_error_control = false
step_size_control = false
-- ||
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
-- ||
-- A = { 0.,  0.,
--       2.,  0.,
--      3/4, 1/4}
-- c = {0., 2., 1.}
-- d = {3/4, 1/4}
-- order  = 2
-- stages = 2
-- stages_bar = 2
-- local_error_control = false
-- step_size_control = false
-- ]]

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

-- -- Automatic step size control HEDOP5 (kirchhoff = 1 --> constrained)
-- A = {   0.,          0.,         0.,         0.,          0.,    0., 0.,
--        1/5,          0.,         0.,         0.,          0.,    0., 0.,
--       3/40,        9/40,         0.,         0.,          0.,    0., 0.,
--      44/45,      -56/15,       32/9,         0.,          0.,    0., 0.,
-- 19372/6561, -25360/2187, 64448/6561,   -212/729,          0.,    0., 0.,
--  9017/3168,     -355/33, 46732/5247,     49/176, -5103/18656,    0., 0.,
--     35/384,          0.,   500/1113,    125/192,  -2187/6784, 11/84, 0.,
--     -18611506045861/19738176307200, 59332529/14479296,   -2509441598627/893904224850,    2763523204159/3289696051200,  -41262869588913/116235927142400, 46310205821/287848404480, -3280/75413}
-- c = {0., 1/5, 3/10, 4/5, 8/9, 1., 1., 19/20}
-- b = {5179/57600, 0., 7571/16695, 393/640, -92097/339200, 187/2100, 1/40}
-- d = {0., 0., 0., 0., 0., 0., 1.}
-- order  = 5
-- order_step_control = 4
-- stages = 6
-- stages_bar = 7
-- local_error_control = false
-- step_size_control = false

-- Use constant mass matrix
const_mass_matrix = 1
-- Use diagonal mass matrix
diag_mass_matrix = 1
-- Use banded solvers for the iteration matrix
banded_iteration_matrix = 1
-- Recalculate the iteration matrix in ever Newton step
recalc_iteration_matrix = 0
-- Perturb initial values (only applies to the constrained case with direct integration of the index-3 formulation)
perturb   = 1
perturb_s = 1.0
-- Use numerical approximation for stiffness matrix K
use_num_K = 1
-- Use numerical approximation for damping matrix D
use_num_D = 1
-- Omit stiffness matrix K in the iteration matrix
no_K = 1
-- Omit damping matrix D in the iteration matrix
no_D = 1

-- Relative error bound for the Newton-Raphson method
rtol = [[((tol)) 1.0e-6 ]]
-- Absolute error bound for the Newton-Raphson method
atol = [[((tol)) 1.0e-8 ]]

-- Integration interval and step size
t0 = 0
te = 10
h  = 0.001
-- steps = math.ceil((te-t0)*2^[--[ 7 || 8 || 9 || 10 || 11 || 12 || 18 ]])
-- steps = 2^([--[ 5 || 6 || 7 || 8 || 9 ]])
-- steps = math.ceil((te-t0)*2^10)
steps = math.ceil((te-t0) * h^-1)

-- Use stabilized index-2 formulation (only applies to the constrained case)
stab2 = 0
-- Baumgarte parameter
my_C = 0.
update_a = false
a_baumgarte = my_C / h
-- Use iterative projection
stab_proj = 1
-- Maximum unsuccessful iteration steps after which the method is considered not to converge
imax = 1

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Problem options   -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
problem_name = 'roll-up'

-- Kirchhoff model
kirchhoff =  1
---- inextensible model
inextensible = 0

-- Calculate number of subdiagonals and superdiagonals of the iteration matrix
additional_subdiag = 0
if kirchhoff then
   additional_subdiag = 2
end
if inextensible then
  additional_subdiag = additional_subdiag + 1
end
if stab2 then
   additional_subdiag = 2*additional_subdiag
end
-- In the unconstrained case there are 11 subdiagonals
nr_subdiag = 11 + additional_subdiag
-- The number of sub- and superdiagonals is equal
nr_superdiag = nr_subdiag

-- Number of discretization points minus one (since we have q_0,..,q_n)
n = 8

-- Length
L = 10

-- Dissipative material constants
CGamd = { 5.0e1,
          5.0e1,
          5.0e1 }
CKd = { 5.0e1,
        5.0e1,
        5.0e1 }
-- Material properties
-- CGam = [ GA, GA, EA ]
CGam = { 1.0e4,
         1.0e4,
         1.0e4}
-- CK = [ EI, EI, GI ]
CK   = { 5.0e2,
         5.0e2,
         5.0e2}
-- Difference between two discretization points
ds = L/n
-- Mass of beam segment
-- m = A*rho * ds
m = 1.0 * ds
-- Inertial mass of a beam segment
-- mI = iner * rho * ds
mI = { 10.0 * ds,
       10.0 * ds,
       10.0 * ds }

-- -- -- Initial values -- -- --
-- Note that  0 <= s <= 1, independent of the beam length
-- Initial positions
function x0(s)
   return { 10.0*s,
            0.0,
            0.0 }
end
-- SIEHE UNTEN DEBUG
---- Initial velocities
function V0(s)
  return { 0, 0, 0} -- DEBUG
end

-- Helper functions
function cross(x,y)
   return { x[2]*y[3] - x[3]*y[2],
            x[3]*y[1] - x[1]*y[3],
            x[1]*y[2] - x[2]*y[1] }
end

function normalize(x)
   local norm = 0
   for i = 1, #x do
      norm = norm + x[i]*x[i]
   end
   norm = math.sqrt(norm)
   for i = 1, #x do
      x[i] = x[i]/norm
   end
   return x
end

-- Initial rotations
function p0(s)
   local v = normalize(cross({0, 0, 1},{1, 0, 0}))
   return {    1/math.sqrt(2),
            v[1]/math.sqrt(2),
            v[2]/math.sqrt(2),
            v[3]/math.sqrt(2) }
end
-- Initial angular velocities
function Om0(s)
   return {0,0,0}
end

-- External forces and moments

--external = 'flying_spaghetti'
--external_parameters = {
--   increasing_time = 2.5,
--   decreasing_time = 2.5,
--   maximum_height  = 200,
--   force_factors = {1/10, 0, 0},
--   moment_factors = {0, -1/2, -1}
--}

-- Roll-up
external = 'roll-up'
external_parameters = {
   factor = 2*math.pi*CK[1]/L,
}

-- Gravity (DEBUG)
--external = 'gravity'
--external_parameters = {
--   g = -100
--}

-- Fixing
fixed_x0 = 1
fixed_x0_position = x0(0)
fixed_xn = 0
fixed_xn_position = x0(1)

-- Fixing
fixed_p0 = 1
fixed_p0_orientation = p0(0)
fixed_pn = 0
fixed_pn_orientation = p0(1)

-- -- -- Output options -- -- --
output_t_at = 0
--t_output_at_multiples_of = 1/2^16
t_output_at_multiples_of = 1/2^7
output_s_at = 0
output_s = { 0/8, 1/8, 2/8, 3/8, 4/8, 5/8, 6/8, 7/8, 8/8 }
