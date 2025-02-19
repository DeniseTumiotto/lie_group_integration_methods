function sols = mixed_calc_end_errors(sols, refsol)

    if ~isfield(sols{1}, 'h')
        sols = create_h(sols);
    end

    the_size     = max(size(sols));
    so3r3_index  = {[]; []; []; []};
    so3r3_ref    = {0, 0; 0, 0; 0, 0; 0, 0};
    se3_index    = {[]; []; []; []};
    se3_ref      = {0, 0; 0, 0; 0, 0; 0, 0};
    s3r3_index   = {[]; []; []; []};
    s3r3_ref     = {0, 0; 0, 0; 0, 0; 0, 0};
    s3sdr3_index = {[]; []; []; []};
    s3sdr3_ref   = {0, 0; 0, 0; 0, 0; 0, 0};

    for i = 1:the_size
        if sols{i}.step_control
            warning('Step control not yet implemented!')
            continue
        end
        switch sols{i}.liegroup
        case 1 % SO(3)
            disp('function not implemented for SO(3)')
        case 2 % S³
            disp('function not implemented for S^3')
        case 3 % SO(3) x R³
            so3r3_index{sols{i}.order-1} = [so3r3_index{sols{i}.order-1}, i];
            if so3r3_ref{sols{i}.order-1, 1} == 0 || so3r3_ref{sols{i}.order-1, 2} > sols{i}.h
                so3r3_ref{sols{i}.order-1, 1} = i;
                so3r3_ref{sols{i}.order-1, 2} = sols{i}.h;
            end
        case 4 % S³ x R³    
            s3r3_index{sols{i}.order-1} = [s3r3_index{sols{i}.order-1}, i];
            if s3r3_ref{sols{i}.order-1, 1} == 0 || s3r3_ref{sols{i}.order-1, 2} > sols{i}.h
                s3r3_ref{sols{i}.order-1, 1} = i;
                s3r3_ref{sols{i}.order-1, 2} = sols{i}.h;
            end
        case 5 % SE(3)
            se3_index{sols{i}.order-1} = [se3_index{sols{i}.order-1}, i];
            if se3_ref{sols{i}.order-1, 1} == 0 || se3_ref{sols{i}.order-1, 2} > sols{i}.h
                se3_ref{sols{i}.order-1, 1} = i;
                se3_ref{sols{i}.order-1, 2} = sols{i}.h;
            end
        case 6 % S³ |x R³
            s3sdr3_index{sols{i}.order-1} = [s3sdr3_index{sols{i}.order-1}, i];
            if s3sdr3_ref{sols{i}.order-1, 1} == 0 || s3sdr3_ref{sols{i}.order-1, 2} > sols{i}.h
                s3sdr3_ref{sols{i}.order-1, 1} = i;
                s3sdr3_ref{sols{i}.order-1, 2} = sols{i}.h;
            end
        case 7 % UDQ
        end
    end

    for i = 1:4
        if ~isempty(so3r3_index{i})
            if nargin < 2 || i == 4
                sols(so3r3_index{i}) = calc_errors(sols(so3r3_index{i}),sols{so3r3_ref{i,1}});
            else
                tmp_sols = calc_errors({sols{so3r3_index{i}},refsol},refsol);
                sols(so3r3_index{i}) = tmp_sols(1:end-1);
            end
        end
        if ~isempty(s3r3_index{i})
            if nargin < 2 || i == 4
                sols(s3r3_index{i}) = calc_errors(sols(s3r3_index{i}),sols{s3r3_ref{i,1}});
            else
                tmp_sols = calc_errors({sols{s3r3_index{i}},refsol},refsol);
                sols(s3r3_index{i}) = tmp_sols(1:end-1);
            end
        end
        if ~isempty(se3_index{i})
            if nargin < 2 || i == 4
                sols(se3_index{i}) = calc_errors(sols(se3_index{i}),sols{se3_ref{i,1}});
            else
                tmp_sols = calc_errors({sols{se3_index{i}},refsol},refsol);
                sols(se3_index{i}) = tmp_sols(1:end-1);
            end
        end
        if ~isempty(s3sdr3_index{i})
            if nargin < 2 || i == 4
                sols(s3sdr3_index{i}) = calc_errors(sols(s3sdr3_index{i}),sols{s3sdr3_ref{i,1}});
            else
                tmp_sols = calc_errors({sols{s3sdr3_index{i}},refsol},refsol);
                sols(s3sdr3_index{i}) = tmp_sols(1:end-1);
            end
        end
    end
end