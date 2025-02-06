function solcell = calc_end_errors(solcell, refpattern, from)
% Calculates errors wrt a reference solution determined by a pattern


% Find reference solution
%refsols = load_all_config_and_bin(refpattern);
refsols = filter_solcell(solcell, refpattern);
if numel(refsols) > 1
   error('calc_errors:ambiguousReference',...
         'Reference pattern matches more than one solution');
end

if numel(refsols) == 0
   error('calc_errors:noReference',...
         'Reference pattern does not match any solution');
end
ref = refsols{1};

if ~ref.rslt.finished == 1
   error('calc_errors:refDidntFinish',...
         'Reference solution did not finish');
end



%
if nargin < 3
   % Error functions
   norm2 = @(x) sqrt(sum(x.^2,1)/size(x,1));
   abserr = @(x,xref) max(norm2(x(:,end)-xref(:,end)));
   relerr = @(x,xref) max(norm2(x(:,end)-xref(:,end))./norm2(xref(:,end)));
   errlog = @(x,xref) [log_s3(qp(inv_s3sdr3(xref(1:4)),x(1:4))); T_inv(log_s3(qp(inv_s3sdr3(xref(1:4)),x(1:4)))).' * (x(5:7)-xref(5:7))];
else
   % Error functions
   norm2 = @(x) sqrt(sum(x.^2,1)/size(x,1));
   abserr = @(x,xref) max(norm2(x(:,from:end)-xref(:,from:end)));
   relerr = @(x,xref) max(norm2(x(:,from:end)-xref(:,from:end))./norm2(xref(:,from:end)));
end

% Refconfig
refconfig = rmfield(ref, 'rslt');

% reference q=(p,x)
[refp, refx] = q_to_ps_xs(ref.rslt.q(:,end),ref);

% Calculate errata
for i=1:numel(solcell)
   has_vd = not(...
           strcmp(solcell{i}.integrator,'RATTLie') ...
        || strcmp(solcell{i}.integrator,'SHAKELie') ...
        || strcmp(ref.integrator,'RATTLie')...
        || strcmp(ref.integrator,'SHAKELie')...
        || strcmp(ref.integrator,'half_explicit'));
   if solcell{i}.rslt.finished == false
      solcell{i}.err.abs.q  = NaN;
      solcell{i}.err.abs.v  = NaN;
      if has_vd
         solcell{i}.err.abs.vd = NaN;
      end
      if isfield('l',solcell{i}.rslt)
         solcell{i}.err.abs.l = NaN;
         if solcell{i}.stab2 == 1 && ~strcmp(solcell{i}.integrator,'half_explicit')
            solcell{i}.err.abs.e = NaN;
         end
      end

      solcell{i}.err.rel.q  = NaN;
      solcell{i}.err.rel.v  = NaN;
      if has_vd
         solcell{i}.err.rel.vd = NaN;
      end
      if isfield('l',solcell{i}.rslt)
         solcell{i}.err.rel.l = NaN;
      end
   else
      solcell{i}.err.refconfig = refconfig;

      % error in q
      [ps, xs] = q_to_ps_xs(solcell{i}.rslt.q(:,end),solcell{i});
      tmp_err_q = [];
      for j = 1:solcell{i}.n+1
         tmp_err_q = [tmp_err_q; errlog([ps(:,j);xs(:,j)],[refp(:,j);refx(:,j)])];
      end
      solcell{i}.err.abs.q = max(norm2(tmp_err_q));

      % solcell{i}.err.abs.q  = abserr(solcell{i}.rslt.q, ref.rslt.q);
      solcell{i}.err.abs.v  = abserr(solcell{i}.rslt.v, ref.rslt.v);
      if has_vd
         solcell{i}.err.abs.vd = abserr(solcell{i}.rslt.vd,ref.rslt.vd);
      end
      if isfield(solcell{i}.rslt,'l')
         % Don't consider the error in the starting point, as it is not
         % calculated by all integrators
         solcell{i}.err.abs.l = abserr(solcell{i}.rslt.l(:,2:end),ref.rslt.l(:,2:end));
         if strcmp(solcell{i}.integrator,'RATTLie')
            solcell{i}.err.abs.lm = abserr(solcell{i}.rslt.e(:,2:end), ref.rslt.l(:,2:end));
            %solcell{i}.err.abs.lp = abserr(2*solcell{i}.rslt.l(:,2:end) - solcell{i}.rslt.e(:,2:end), ref.rslt.l(:,2:end));
         else
            if solcell{i}.stab2 == 1 && ~strcmp(solcell{i}.integrator,'half_explicit')
               solcell{i}.err.abs.e = abserr(solcell{i}.rslt.e, zeros(size(solcell{i}.rslt.e)));
            end
         end
      end

      % solcell{i}.err.rel.q  = relerr(solcell{i}.rslt.q, ref.rslt.q);
      solcell{i}.err.rel.v  = relerr(solcell{i}.rslt.v, ref.rslt.v);
      if has_vd
         solcell{i}.err.rel.vd = relerr(solcell{i}.rslt.vd,ref.rslt.vd);
      end
      if isfield(solcell{i}.rslt,'l')
         % Don't consider the error in the starting point, as it is not
         % calculated by all integrators
         solcell{i}.err.rel.l = relerr(solcell{i}.rslt.l(:,2:end),ref.rslt.l(:,2:end));
         if strcmp(solcell{i}.integrator,'RATTLie')
            solcell{i}.err.rel.lm = relerr(solcell{i}.rslt.e(:,2:end), ref.rslt.l(:,2:end));
            %solcell{i}.err.rel.lp = relerr(2*solcell{i}.rslt.l(:,2:end) - solcell{i}.rslt.e(:,2:end), ref.rslt.l(:,2:end));
         end
      end
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
