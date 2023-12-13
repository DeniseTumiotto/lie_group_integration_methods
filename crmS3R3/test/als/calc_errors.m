function solcell = calc_errors(solcell, refpattern)
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



% Error functions
norm2 = @(x) sqrt(sum(x.^2,1)/size(x,1));
abserr = @(x,xref) max(norm2(x-xref));
relerr = @(x,xref) max(norm2(x-xref)./norm2(xref));

% Refconfig
refconfig = rmfield(ref, 'rslt');

% Calculate errata
for i=1:numel(solcell)
   has_vd = not(...
           strcmp(solcell{i}.integrator,'RATTLie') ...
        || strcmp(solcell{i}.integrator,'SHAKELie') ...
        || strcmp(ref.integrator,'RATTLie')...
        || strcmp(ref.integrator,'SHAKELie'));
   if solcell{i}.rslt.finished == false
      solcell{i}.err.abs.q  = NaN;
      solcell{i}.err.abs.v  = NaN;
      if has_vd
         solcell{i}.err.abs.vd = NaN;
      end
      if isfield('l',solcell{i}.rslt)
         solcell{i}.err.abs.l = NaN;
         if solcell{i}.stab2 == 1
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

      solcell{i}.err.abs.q  = abserr(solcell{i}.rslt.q, ref.rslt.q);
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
            solcell{i}.err.abs.lmp = abserr((solcell{i}.rslt.l(:,2:end) + solcell{i}.rslt.e(:,2:end))/2, ref.rslt.l(:,2:end));
            %solcell{i}.err.abs.lp = abserr(2*solcell{i}.rslt.l(:,2:end) - solcell{i}.rslt.e(:,2:end), ref.rslt.l(:,2:end));
         else
            if solcell{i}.stab2 == 1
               solcell{i}.err.abs.e = abserr(solcell{i}.rslt.e, 0);
            end
         end
      end

      solcell{i}.err.rel.q  = relerr(solcell{i}.rslt.q, ref.rslt.q);
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
            solcell{i}.err.rel.lmp = relerr((solcell{i}.rslt.l(:,2:end) + solcell{i}.rslt.e(:,2:end))/2, ref.rslt.l(:,2:end));
            %solcell{i}.err.rel.lp = relerr(2*solcell{i}.rslt.l(:,2:end) - solcell{i}.rslt.e(:,2:end), ref.rslt.l(:,2:end));
         end
      end
   end
end
