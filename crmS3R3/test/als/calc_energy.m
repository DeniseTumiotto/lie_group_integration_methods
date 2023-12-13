function sol = calc_energy(sol)

% Define diagonals of mass and coefficient matrix (which are diagonal
% matrices)
diag_mass_matrix = [sol.mI, sol.m, sol.m, sol.m]';
diag_c_matrix = [sol.CK, sol.CGam]';

% Initialize energies
sol.rslt.kinetic_energy = zeros(size(sol.rslt.t));

sol.rslt.bending_energy = zeros(size(sol.rslt.t));
sol.rslt.extension_energy = zeros(size(sol.rslt.t));
sol.rslt.shearing_energy = zeros(size(sol.rslt.t));
sol.rslt.torsion_energy = zeros(size(sol.rslt.t));
sol.rslt.potential_energy = zeros(size(sol.rslt.t));

% Loop over all times
for i=1:length(sol.rslt.t)
   % Extract current q and v
   q = reshape(sol.rslt.q(:,i), [7, sol.n+1]);
   v = reshape(sol.rslt.v(:,i), [6, sol.n+1]);
   
   % % Kinetic energy
   % Treat ends differently
   sol.rslt.kinetic_energy(i) = ...
         1/2*v(:,  1)'*(diag_mass_matrix/2.*v(:,  1))...
       + 1/2*v(:,end)'*(diag_mass_matrix/2.*v(:,end));
   % Loop over all other discretization points
   for j=2:sol.n
      sol.rslt.kinetic_energy(i) = sol.rslt.kinetic_energy(i) ...
         + 1/2*v(:,j)'*(diag_mass_matrix.*v(:,j));
   end
   
   % % Potential energy
   % Calculate discrete spatial derivatives w
   w = zeros(6,sol.n);
   for j=1:sol.n
      w(:,j) = logt_S3R3(lp_S3R3(inv_S3R3(q(:,j)),q(:,j+1)))/sol.ds;
   end
   % Calculate potential energy
   for j=1:sol.n
      sol.rslt.shearing_energy(i) = sol.rslt.shearing_energy(i) ...
         + sol.ds/2 * w(1:2,j)' * (diag_c_matrix(1:2) .* w(1:2,j));
      sol.rslt.torsion_energy(i) = sol.rslt.torsion_energy(i) ...
         + sol.ds/2 * w(3,j)' * (diag_c_matrix(3) .* w(3,j));
      sol.rslt.bending_energy(i) = sol.rslt.bending_energy(i) ...
         + sol.ds/2 * w(4:5,j)' * (diag_c_matrix(4:5) .* w(4:5,j));
      sol.rslt.extension_energy(i) = sol.rslt.extension_energy(i) ...
         + sol.ds/2 * (w(6,j) - 1)' * (diag_c_matrix(6) .* (w(6,j) - 1));
   end
  
end

sol.rslt.potential_energy = sol.rslt.shearing_energy + sol.rslt.torsion_energy + sol.rslt.bending_energy + sol.rslt.extension_energy;

sol.rslt.energy = sol.rslt.kinetic_energy + sol.rslt.potential_energy;

function rslt = qp(p,q)

rslt(1)     = p(1)*q(1) - p(2:4)'*q(2:4);
rslt(2:4,1) = p(1)*q(2:4) + q(1)*p(2:4) + cross(p(2:4),q(2:4));

function rslt = apply_p(p,v)

pv = qp(p,[0;v]);
rslt(1:3,1) = -pv(1)*p(2:4) + p(1)*pv(2:4) - cross(pv(2:4),p(2:4));

function rslt = lp_S3R3(q1,q2)

rslt(1:4,1) = qp(q1(1:4), q2(1:4));
rslt(5:7)   = q1(5:7) + apply_p(q1(1:4), q2(5:7));

function rslt = inv_S3R3(q)

rslt(1:4,1) = [q(1); -q(2:4)];
rslt(5:7) = - apply_p(rslt(1:4), q(5:7));

function rslt = skw(v)

rslt = [ 0.0,  -v(3), v(2);
         v(3), 0.0,  -v(1);
        -v(2),  v(1), 0.0];

function rslt = Tan_inv_S3(Om)

EOm = skw(Om);
nOm = norm(Om);

if (nOm < 3.0e-3)
   f1 = 1.0/360 + nOm^2/7560 + nOm^4/201600;
else
   f1 = (2 - nOm / tan(nOm/2))/(2*nOm^2);
end
rslt = eye(3) + EOm/2 + f1*EOm*EOm;

function rslt = logt_S3R3(q)

if true
   if (abs(q(1)) > 1 - 1e-7)
      f1 = 2 - 2/3*(q(1)-1) + 4/15*(q(1)-1)^2 - 4/35*(q(1)-1)^3 + 16/315*(q(1)-1)^4;
   else
      f1 = 2*acos(q(1))/sqrt(1-q(1)^2);
   end

   rslt(1:3,1) = f1*q(2:4);
   rslt(4:6,1) = Tan_inv_S3(-rslt(1:3))*q(5:7);
else
   nrm = norm(q(2:4));
   if (nrm < 1.0e-12)
      rslt(1:3,1) = 0.0;
      rslt(4:6,1) = q(5:7);
   else
      rslt(1:3,1) = 2*acos(min(max(-1.0,q(1)),1.0))*q(2:4)/nrm;
      rslt(4:6,1) = Tan_inv_S3(-rslt(1:3))*q(5:7);
   end
end