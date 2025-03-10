function sols = mixed_calc_end_errors(sols, refsol)

    if ~isfield(sols{1}, 'h')
        sols = create_h(sols);
    end

    the_size = max(size(sols));
    rll_index = {[]; []; []; []};
    rll_ref   = {0, 0; 0, 0; 0, 0; 0, 0};
    fly_index = {[]; []; []; []};
    fly_ref   = {0, 0; 0, 0; 0, 0; 0, 0};

    for i = 1:the_size
        if sols{i}.step_control
            warning('Step control not yet implemented!')
            continue
        end
        if strcmp(sols{i}.problem_name, 'roll-up')
            rll_index{sols{i}.order-1} = [rll_index{sols{i}.order-1}, i];
            if rll_ref{sols{i}.order-1, 1} == 0 || rll_ref{sols{i}.order-1, 2} > sols{i}.h
                rll_ref{sols{i}.order-1, 1} = i;
                rll_ref{sols{i}.order-1, 2} = sols{i}.h;
            end
        elseif strcmp(sols{i}.problem_name, 'flying_spaghetti')
            fly_index{sols{i}.order-1} = [fly_index{sols{i}.order-1}, i];
            if fly_ref{sols{i}.order-1, 1} == 0 || fly_ref{sols{i}.order-1, 2} > sols{i}.h
                fly_ref{sols{i}.order-1, 1} = i;
                fly_ref{sols{i}.order-1, 2} = sols{i}.h;
            end
        end
    end

    for i = 1:4
        if ~isempty(rll_index{i})
            if nargin < 2 || i == 4
                sols(rll_index{i}) = calc_end_errors(sols(rll_index{i}),sols{rll_ref{i,1}});
            else
                tmp_sols = calc_end_errors({sols{rll_index{i}},refsol},refsol);
                sols(rll_index{i}) = tmp_sols(1:end-1);
            end
        end
        if ~isempty(fly_index{i})
            if nargin < 2
                sols(fly_index{i}) = calc_end_errors(sols(fly_index{i}),sols{fly_ref{i,1}});
            else
                sols(fly_index{i}) = calc_end_errors(sols(fly_index{i}),refsol);
            end
        end
    end
end