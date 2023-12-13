-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Integrator options   -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- Algorithmic paramters for gena(for rho_infinity = 0.9)
rho_infinity = 0.3
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

-- Relative error bound for the Newton-Raphson method
rtol = [[((tol)) 1.0e-6 ]]
-- Absolute error bound for the Newton-Raphson method
atol = [[((tol)) 1.0e-8 ]]
-- Maximum unsuccessful iteration steps after which the method is considered not to converge
imax = 200

-- Integration interval and step size
t0 = 0
te = 20
--steps = (te - t0) * 2^ [--[ 16 || 15 || 14 || 13 || 12 || 11 || 10 || 9 ]]
steps = math.ceil((te - t0) * 2^11)

--[--[
-- Use stabilized index-2 formulation (only applies to the constrained case)
stab2 = 1
-- Omit stiffness matrix K and damping matrix D in the iteration matrix
no_K = 1
no_D = 1
--||
---- Use stabilized index-2 formulation (only applies to the constrained case)
--stab2 = 1
---- Omit stiffness matrix K and damping matrix D in the iteration matrix
--no_K = 1
--no_D = 1
--]]

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Problem options   -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
problem_name = 'vibrate'

-- Kirchhoff model
kirchhoff =  0
---- inextensible model
--inextensible = 0

-- Calculate number of subdiagonals and superdiagonals of the iteration matrix
additional_subdiag = 0
if kirchhoff then
   additional_subdiag = 2
end
--if inextensible then
--   additional_subdiag = additional_subdiag + 1
--end
if stab2 then
   additional_subdiag = 2*additional_subdiag
end
-- In the unconstrained case there are 11 subdiagonals
nr_subdiag = 11 + additional_subdiag
-- The number of sub- and superdiagonals is equal
nr_superdiag = nr_subdiag

-- Number of discretization points minus one (since we have q_0,..,q_n)
--n = 2^[--[ 8 || 7 || 6 || 5 || 4 || 3 || 2 ]]
--n = [--[ 16 || 8 ]]
n = 16
-- Length
L = 10

-- Dissipative material constants
CGamd = { 0.0e1,
          0.0e1,
          0.0e1 }
CKd = { 0.0e1,
        0.0e1,
        0.0e1 }
-- Material properties
-- CGam = [ GA, GA, EA ]
CGam = { 1.0e4,
         1.0e4,
         1.0e6}
-- CK = [ EI, EI, GI ]
CK   = { 5.0e4,
         5.0e4,
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
   --return { 10.0*s,
   --         0.3*s*(s-1),
   --         0.4*s*(s-0.5)*(s-1) }
   return { 11.0*s, 0, 0}
end
-- SIEHE UNTEN DEBUG
---- Initial velocities
--function V0(s)
--   return { 0, 0, 1} -- DEBUG
--end

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
   --local v = normalize(cross({0, 0, 1},{1, 0, 0}))
   --return {    1/math.sqrt(2),
   --         v[1]/math.sqrt(2),
   --         v[2]/math.sqrt(2),
   --         v[3]/math.sqrt(2) }
   return normalize({1,0,1+s,0})
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
--external = 'roll-up'
--external_parameters = {
--   factor = 2*math.pi*CK[1]/L,
--}

-- Gravity
external = 'gravity'
external_parameters = {
   g = 0
}

-- Fixing
fixed_x0 = 0
fixed_x0_position = x0(0)
fixed_xn = 0
fixed_xn_position = x0(1)

-- Fixing
fixed_p0 = 0
fixed_p0_orientation = p0(0)
fixed_pn = 0
fixed_pn_orientation = p0(1)

-- -- -- Output options -- -- --
output_t_at = 0
t_output_at_multiples_of = 1/2^7
output_s_at = 0
output_s = { 0/8, 1/8, 2/8, 3/8, 4/8, 5/8, 6/8, 7/8, 8/8 }





--
function qp(p,q)
   return {p[1]*q[1] - p[2]*q[2] - p[3]*q[3] - p[4]*q[4],
           p[2]*q[1] + p[1]*q[2] - p[4]*q[3] + p[3]*q[4],
           p[3]*q[1] + p[4]*q[2] + p[1]*q[3] - p[2]*q[4],
           p[4]*q[1] - p[3]*q[2] + p[2]*q[3] + p[1]*q[4]}
end
function conj_p(p)
   return {p[1], -p[2],-p[3],-p[4]}
end
function apply_p(p,v)
   tmp = {0, v[1], v[2], v[3]}
   tmp = qp(p,tmp)
   tmp = qp(tmp,conj_p(p))
   out = {tmp[2], tmp[3], tmp[4]}
   return out
end


-- Initial velocities
function V0(s)
   --return apply_p(conj_p(p0(s)), {0,  0 , 2*(s-0.5)})
   return {0,0,0}
end
