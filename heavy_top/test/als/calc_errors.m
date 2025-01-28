function solcell = calc_errors(solcell, refpattern)
% Calculates errors wrt a reference solution determined by a pattern

% Find reference solution
foundref = false;
for i=1:numel(solcell)
   if structmatch(solcell{i}, refpattern)
      if foundref
         error('calc_errors:ambiguousReference',...
            'Reference pattern matches more than one solution');
      else
         foundref = true;
         ref = solcell{i};
      end
   end
end
if ~foundref
   error('calc_errors:noReference',...
      'Reference pattern does not match any solution');
end
if ~ref.rslt.finished == 1
   error('calc_errors:refDidntFinish',...
      'Reference solution did not finish');
end

% Error functions
norm2 = @(x) sqrt(sum(x.^2,1)/size(x,1));
abserr = @(x,xref) max(norm2(x(:,end)-xref(:,end)));
relerr = @(x,xref) max(norm2(x(:,end)-xref(:,end))./norm2(xref(:,end)));
abslog = @(x,xref) max(norm2([log_s3(qp(inv_s3sdr3(xref(1:4,end)),x(1:4,end))); T_inv(log_s3(qp(inv_s3sdr3(xref(1:4,end)),x(1:4,end)))).' * (x(5:7,end)-xref(5:7,end))]));
rellog = @(x,xref) max(norm2([log_s3(qp(inv_s3sdr3(xref(1:4,end)),x(1:4,end))); T_inv(log_s3(qp(inv_s3sdr3(xref(1:4,end)),x(1:4,end)))).' * (x(5:7,end)-xref(5:7,end))])./norm2(xref(:,end)));

% Refconfig
refconfig = rmfield(ref, 'rslt');

% Calculate errata
for i=1:numel(solcell)
   if solcell{i}.rslt.finished == false
      solcell{i}.err.abs.q  = NaN;
      solcell{i}.err.abs.v  = NaN;
      if isfield('vd',solcell{i}.rslt)
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
      if isfield('vd',solcell{i}.rslt)
          solcell{i}.err.rel.vd = NaN;
      end
      if isfield('l',solcell{i}.rslt)
         solcell{i}.err.rel.l = NaN;
      end
   else
      solcell{i}.err.refconfig = refconfig;

      % solcell{i}.err.abs.q  = abserr(solcell{i}.rslt.q, ref.rslt.q);
      solcell{i}.err.abs.q  = abslog(solcell{i}.rslt.q, ref.rslt.q);
      solcell{i}.err.abs.v  = abserr(solcell{i}.rslt.v, ref.rslt.v);
      if isfield('vd',solcell{i}.rslt)
          solcell{i}.err.abs.vd = abserr(solcell{i}.rslt.vd,ref.rslt.vd);
      end
      if isfield(solcell{i}.rslt,'l')
         solcell{i}.err.abs.l = abserr(solcell{i}.rslt.l,ref.rslt.l);
         if solcell{i}.stab2 == 1
            solcell{i}.err.abs.e = abserr(solcell{i}.rslt.e, 0);
         end
      end

      % solcell{i}.err.rel.q  = relerr(solcell{i}.rslt.q, ref.rslt.q);
      solcell{i}.err.rel.q  = rellog(solcell{i}.rslt.q, ref.rslt.q);
      solcell{i}.err.rel.v  = relerr(solcell{i}.rslt.v, ref.rslt.v);
      if isfield('vd',solcell{i}.rslt)
          solcell{i}.err.rel.vd = relerr(solcell{i}.rslt.vd,ref.rslt.vd);
      end
      if isfield(solcell{i}.rslt,'l')
         solcell{i}.err.rel.l = relerr(solcell{i}.rslt.l,ref.rslt.l);
      end
   end
end
   
