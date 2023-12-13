function out = get_SO3xR3_Om(q, v, sol)
% Calculates angular velocity Omega from configuration q, v given in
% a particular Lie group.

switch sol.liegroup
   case {1, 2, 3, 4, 5, 6, 7}
      out = v(1:3);
   otherwise
      error('get_SO3xR3_R:unknownLieGroup', 'Unknown Lie group');
end      
