function solcell = calc_SO3xR3_errors(solcell, refcell)
% Calculates errors wrt a reference solution
% The errors are measured after converting into the Lie group SO3xR3

% Get reference
if numel(refcell) > 1
   error('calc_SO3xR3_errors:ambiguousReference',...
            'Reference pattern matches more than one solution');
elseif numel(refcell) == 0
   error('calc_SO3xR3_errors:noReference',...
      'Reference pattern does not match any solution');
end

ref = refcell{1};


if ~ref.rslt.finished == 1
   error('calc_SO3xR3_errors:refDidntFinish',...
      'Reference solution did not finish');
end
%if ~ref.liegroup == 3 % SO(3) x R³
%   error('calc_SO3xR3_errors:wrongReferenceLieGroup',...
%      'Lie group formulation of reference is not SO(3)xR³');
%end

% Error functions
norm2 = @(x) sqrt(sum(x.^2,1)/size(x,1));
abserr = @(x,xref) max(norm2(x-xref));
relerr = @(x,xref) max(norm2(x-xref)./norm2(xref));

id = @(x) x;

% Refconfig
refconfig = rmfield(ref, 'rslt');

% Calculate errata
for i=1:numel(solcell)
   
   if solcell{i}.rslt.finished == false
      solcell{i}.SO3xR3_err.abs.q  = NaN;
      solcell{i}.SO3xR3_err.abs.v  = NaN;
      solcell{i}.SO3xR3_err.abs.vd = NaN;
      if isfield('l',solcell{i}.rslt) && isfield('l',ref.rslt)
         solcell{i}.SO3xR3_err.abs.l = NaN;
         if solcell{i}.stab2 == 1
            solcell{i}.SO3xR3_err.abs.e = NaN;
         end
      end

      solcell{i}.SO3xR3_err.rel.q  = NaN;
      solcell{i}.SO3xR3_err.rel.v  = NaN;
      solcell{i}.SO3xR3_err.rel.vd = NaN;
      if isfield('l',solcell{i}.rslt) && isfield('l',ref.rslt)
         solcell{i}.SO3xR3_err.rel.l = NaN;
      end
   else
      solcell{i}.SO3xR3_err.refconfig = refconfig;

      solcell{i}.SO3xR3_err.abs.q  = abserr(get_SO3xR3_q(solcell{i}.rslt.q, solcell{i}), ...
                                            get_SO3xR3_q(ref.rslt.q, ref));
      solcell{i}.SO3xR3_err.abs.v  = abserr(get_SO3xR3_v(solcell{i}.rslt.q, solcell{i}.rslt.v, solcell{i}), ...
                                            get_SO3xR3_v(ref.rslt.q, ref.rslt.v, ref));
      if isfield(solcell{i}.rslt,'l') && isfield('l',ref.rslt)
         solcell{i}.SO3xR3_err.abs.l = abserr(solcell{i}.rslt.l,ref.rslt.l);
         if solcell{i}.stab2 == 1
            solcell{i}.SO3xR3_err.abs.e = abserr(solcell{i}.rslt.e, 0);
         end
      end

      solcell{i}.SO3xR3_err.rel.q  = relerr(get_SO3xR3_q(solcell{i}.rslt.q, solcell{i}), ...
                                            get_SO3xR3_q(ref.rslt.q, ref));
      solcell{i}.SO3xR3_err.rel.v  = relerr(get_SO3xR3_v(solcell{i}.rslt.q, solcell{i}.rslt.v, solcell{i}), ...
                                            get_SO3xR3_v(ref.rslt.q, ref.rslt.v, ref));   
      if isfield(solcell{i}.rslt,'l') && isfield('l',ref.rslt)
         solcell{i}.SO3xR3_err.rel.l = relerr(solcell{i}.rslt.l,ref.rslt.l);
      end
   end
end
   
