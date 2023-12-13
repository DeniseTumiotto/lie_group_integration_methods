-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Integrator options   -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- Algorithmic paramters (for rho_infinity = 0.9)
alpha_m =   8/ 19
alpha_f =   9/ 19
beta    = 100/361
gamma   =  21/ 38

-- Use constant mass matrix
const_mass_matrix = 1
-- Use diagonal mass matrix
diag_mass_matrix = 1
-- Use banded solvers for the iteration matrix
banded_iteration_matrix = 1
-- Recalculate the iteration matrix in ever Newton step
recalc_iteration_matrix = 0
-- Perturb initial values (only applies to the constrained case with direct integration of the index-3 formulation)
perturb   = 0
perturb_s = 1.0
-- Use numerical approximation for stiffness matrix K
use_num_K = 1
-- Use numerical approximation for damping matrix D
use_num_D = 1
-- Omit stiffness matrix K in the iteration matrix
no_K = 0
-- Omit damping matrix D in the iteration matrix
no_D = 0

-- Relative error bound for the Newton-Raphson method
rtol = 1.0e-6
-- Absolute error bound for the Newton-Raphson method
atol = 1.0e-8
-- Maximum unsuccessful iteration steps after which the method is considered not to converge
imax = 100

-- Integration interval and step size
t0 = 0
te = 1
steps = 2^12

-- Use stabilized index-2 formulation (only applies to the constrained case)
stab2 = 0

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Problem options   -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
problem_name = 'pendulum'

-- Kirchhoff model
kirchhoff = 0
-- inextensible model
inextensible = 0

-- Simplified model assumptions
use_SMA = 0

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
-- In the unconstrained case there are 9 subdiagonals
nr_subdiag = 9 + additional_subdiag
-- The number of sub- and superdiagonals is equal
nr_superdiag = nr_subdiag

-- Number of discretization points of the position
N = 2^4

-- Length
L = 1
-- Density
rho = 1.1e3
-- Elastic modulus
E = 5.0e6
-- Poisson number
nu = 0.5
-- Shear modulus
G = E/(2*(1+nu))
-- Radius
r = 5.0e-3
-- Area
A = math.pi* r^2
-- Inertia terms
iner = { r^4 * math.pi/4,
         r^4 * math.pi/4,
         r^4 * math.pi/2 }
-- Dissipative material constants
CGamd = { 1.0e-1,
          1.0e-1,
          2.0e2 }
CKd = { 2.0e-4,
        2.0e-4,
        8.0e-6 }
-- Timoshenko shear correction factors
kappa = { 6*(1+nu)/(7+6*nu),
          6*(1+nu)/(7+6*nu) }
-- Material properties
CGam = { A*G*kappa[1],
         A*G*kappa[2],
         A*E }
CK   = { E*iner[1],
         E*iner[2],
         G*iner[3] }
-- Difference between two discretization points
ds = L/(N-1)
-- Mass of beam segment
m = rho * A * ds
-- Inertial mass of a beam segment
mI = { rho*iner[1]*ds,
       rho*iner[2]*ds,
       rho*iner[3]*ds }

-- -- -- Initial values -- -- --
-- Note that  0 <= s <= 1, independent of the beam length
-- Initial positions
function x0(s)
   return { s*L, 0, 0 }
end
-- Initial velocities
function xd0(s)
   return { 0, 0, 0 }
end
-- Initial rotations
function p0(s)
   return { math.sqrt(0.5), 0, math.sqrt(0.5), 0 }
end
-- Initial angular velocities
function Om0(s)
   return { 0, 0, 0 }
end

-- External forces and moments
external = 'gravity'
external_parameters = { g = -9.81 }

-- Fixing
fixed_x1 = 1
fixed_x1_position = { 0, 0, 0 }
fixed_xN = 0
fixed_xN_position = nil

-- -- -- Output options -- -- --
output_t_at = 1
t_output_at_multiples_of = 1/2^7
output_xp_at = 0
s_output_x = { 0.0,        0.25,        0.5,        0.75,        1.0 }
s_output_p = {      0.125,       0.375,      0.625,       0.875      }
