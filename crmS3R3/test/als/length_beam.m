function sols = length_beam(sols)
    
    the_size = max(size(sols));
    
    for i = 1:the_size

        time_len = length(sols{i}.rslt.q(1,:));
        rslt = zeros(time_len, 1);
    
        for t = 1:time_len
            xx = q_to_xs(sols{i}.rslt.q(:,t),sols{i});
            for j = 2:sols{i}.n + 1
                rslt(t) = rslt(t) + sqrt(sum((xx(:,j)-xx(:,j-1)).^2));
            end
        end

        sols{i}.length = rslt;

    end

end

function xs = q_to_xs(q,sol)

    xs = zeros(3,sol.n+1);
    
    if sol.output_s_at == 0
    
       if sol.fixed_x0 == 0
          xs(:,1) = q(5:7);
          q_start = 8;
       else
          xs(:,1) = sol.fixed_x0_position;
          q_start = 5;
       end
       if sol.fixed_p0 == 1
          q_start = q_start - 4;
       end
    
       if sol.fixed_xn == 0
          xs(:,end) = q(end-2:end);
       else
          xs(:,end) = sol.fixed_xn_position;
       end
    
       tmp = reshape(q(q_start:q_start-1 + 7*(sol.n-1)),[7,sol.n-1]);
       xs(:,2:end-1) = tmp(5:7,:);
       
    else
       tmp = reshape(q,[7,length(sol.output_s)]);
       xs = tmp(5:7,:);
    end
   
end