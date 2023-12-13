function out = get_SO3xR3_x(q, sol)
% Calculates position x from configuration q given in
% a particular Lie group.

if sol.liegroup == 1 || sol.liegroup == 2
   refX = qpvkp(qconj(sol.p0), sol.x0')';
end

switch sol.liegroup
   case 1 % SO(3)
      out = reshape(q(1:9),[3,3]) * refX;
   case 2 % S^3
      out = qpvkp(q(1:4), refX)';
   case {3, 4, 5, 6} % SO(3)xR^3, S^3xR^3, SE(3), S^3|xR^3
      out = q(end-2:end);
   case 7 % UDQ
      out = 2*qIm(qp(q(5:8),qconj(q(1:4))))';
   otherwise
      error('get_SO3xR3_R:unknownLieGroup', 'Unknown Lie group');
end      
