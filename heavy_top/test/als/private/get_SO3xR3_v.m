function out = get_SO3xR3_v(q, v, sol)
% Calculates the velocity of the configurations in the Lie group SO(3)xR³

out = zeros(6, size(q,2));
for i=1:size(q,2)
   out(:,i) = [get_SO3xR3_Om(q(:,i), v(:,i), sol);
               get_SO3xR3_u(q(:,i), v(:,i), sol)];
end
