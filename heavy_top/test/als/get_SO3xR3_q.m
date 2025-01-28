function out = get_SO3xR3_q(q, sol)
% Calculates the configurations in the Lie group SO(3)xR³

out = zeros(12, size(q,2));
for i=1:size(q,2)
   out(:,i) = [get_SO3xR3_R(q(:,i), sol);
               get_SO3xR3_x(q(:,i), sol)];
end
