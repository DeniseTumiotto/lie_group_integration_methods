function out = get_SO3xR3_u(q, v, sol)
% Calculates velocity u from configuration q, v given in
% a particular Lie group.

if sol.liegroup == 1 || sol.liegroup == 2
   refX = qpvkp(qconj(sol.p0), sol.x0')';
end

switch sol.liegroup
   case 1 % SO(3)
      out = reshape(q(1:9),[3 3]) * cross(v(1:3), refX);
   case 2 % S³
      out = qpvkp(q(1:4), cross(v(1:3), refX))';
   case {3, 4} % SO(3)xR³, S³xR³
      out = v(4:6);
   case 5 % SE(3)
      out = reshape(q(1:9),[3 3]) * v(4:6);
   case {6, 7} % S³|xR³, UDQ
      out = qpvkp(q(1:4), v(4:6))';
   otherwise
      error('get_SO3xR3_R:unknownLieGroup', 'Unknown Lie group');
end      
