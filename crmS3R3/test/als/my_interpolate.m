function s = my_interpolate(s,N)
% function to interpolate the results in more segments than the actual one
% it is valid only for multiple of the desired discretization
len_s = max(size(s));

for j = 1:len_s
    if mod(N,s{j}.n)~=0
        disp('Interpolation not possible')
        return
    else
        times = N/s{j}.n;
    end

    [ps,xs] = q_to_ps_xs(s{j}.rslt.q,s{j});
    
    interval = linspace(0,1,times+1);
    ps_interp = zeros(4,N+1);
    xs_interp = zeros(3,N+1);
    s{j}.rslt.q_interp = zeros(7*(N+1),1);
    
    % first point
    ps_interp(:,1) = ps(:,1);
    xs_interp(:,1) = xs(:,1);

    for i = 1:s{j}.n
    % segment between two consecutive nodes
    % ps_interp(:,times*(i-1)+2:times*(i-1)+times+1) = interval(2:end) .* ps(:,i+1) + (ones(1,times)-interval(2:end)) .* ps(:,i);
    ps_interp(:,times*(i-1)+2:times*(i-1)+times+1) = quatinterp(ones(times,1).*ps(:,i).', ones(times,1).*ps(:,i+1).', interval(2:end).').';
    xs_interp(:,times*(i-1)+2:times*(i-1)+times+1) = interval(2:end) .* xs(:,i+1) + (ones(1,times)-interval(2:end)) .* xs(:,i);
    end

    for i = 1:N+1
        s{j}.rslt.q_interp(7*(i-1)+1:7*(i-1)+7) = [ps_interp(:,i); xs_interp(:,i)];
    end

end
end

function [ps, xs] = q_to_ps_xs(q,sol)

    ps = zeros(4,sol.n+1);
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
            ps(:,1) = sol.fixed_p0_orientation;
            q_start = q_start - 4;
        else
            ps(:,1) = q(1:4);
        end
        if sol.fixed_xn == 0
            xs(:,end) = q(end-2:end);
            ps(:,end) = q(end-6:end-3);
            q_end = length(q) - 7;
        else
            xs(:,end) = sol.fixed_xn_position;
            ps(:,end) = q(end-3:end);
            q_end = length(q) - 4;
        end
        if sol.fixed_pn == 1
            ps(:,end) = sol.fixed_pn_orientation;
            q_end = q_end - 4;
        end

        tmp = reshape(q(q_start:q_start-1 + 7*(sol.n-1)),[7,sol.n-1]);
        ps(:,2:end-1) = tmp(1:4,:);
        xs(:,2:end-1) = tmp(5:7,:);
        else
        tmp = reshape(q,[7,length(sol.output_s)]);
        ps = tmp(1:4,:);
        xs = tmp(5:7,:);
        end
end