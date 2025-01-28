function out = get_SO3xR3_R(q, sol)
% Calculates rotation matrix R (as vector) from configuration q given in
% a particular Lie group.

switch sol.liegroup
   case {1, 3, 5} % SO(3), SO(3)xR^3, SE(3)
      out = q(1:9);
   case {2, 4, 6, 7} % S^3, S^3xR^3, S^3|xR^3, UDQ
      out = reshape(qR(q(1:4)),[9 1]);
   otherwise
      error('get_SO3xR3_R:unknownLieGroup', 'Unknown Lie group');
end      
