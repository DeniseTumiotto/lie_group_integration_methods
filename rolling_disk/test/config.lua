-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Integrator options   -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


-- Use constant mass matrix
const_mass_matrix = 1
-- Use diagonal mass matrix
diag_mass_matrix = 1
-- Use banded solvers for the iteration matrix
banded_iteration_matrix = 0
-- Recalculate the iteration matrix in ever Newton step
recalc_iteration_matrix = 0

-- banded_iteration_matrix options (not used)
nr_superdiag = -1
nr_subdiag = -1

-- Relative error bound for the Newton-Raphson method
rtol = [[((tol)) 1.0e-6 ]]
-- Absolute error bound for the Newton-Raphson method
atol = [[((tol)) 1.0e-8 ]]
-- Maximum unsuccessful iteration steps after which the method is considered not to converge
imax = 200

-- Integration interval and step size
t0 = 0
--te = 2^- [--[ 12 || 21 || 20 || 19 || 18 || 17 || 16 || 15 || 14 || 13 || 11 || 10 || 9 || 8 ]]
te = 2
steps = math.ceil((te - t0)) * 2^ [[ 12 || 21 || 20 || 19 || 18 || 17 || 16 || 15 || 14 || 13 || 11 || 10 || 9 || 8 ]]
--steps = 2^ [--[ 12 || 21 || 20 || 19 || 18 || 17 || 16 || 15 || 14 || 13 || 11 || 10 || 9 || 8 ]]
--steps = math.ceil((te - t0) * 2^10)
--steps = 1

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- Problem options   -- -- -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
problem_name = 'rolling_disk'

-- Mass of the top
mass = 15

-- Gravity including direction
gravity = {0, 0, -9.81}
--gravity = {0,0,0}

-- Diagonal elements of inertial tensor wrt. the center of mass
inerJ = {1/4*mass*1^2, 1/4*mass*1^2, 1/2*mass*1^2}


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
   return {tmp[2], tmp[3], tmp[4]}
end

-- -- -- Initial values
-- Initial rotation (in S^3)
th = 0 /360 *2*math.pi
p0 = qp({math.cos(th/2),math.sin(th/2),0,0},{0, 0, 1/math.sqrt(2), 1/math.sqrt(2)})
-- Initial positions (in R^3)
x0 = {0, 0, 2*math.sqrt((p0[2]^2 + p0[3]^2)*(p0[1]^2 + p0[4]^2))}
-- Initial angular velocity
Om0 = {0, 15, 5}
-- Initial velocity is chosen to be consistent

-- -- -- Output options -- -- --
output_t_at = 1
t_output_at_multiples_of = 1/2^7
