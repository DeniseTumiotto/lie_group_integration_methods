function solcell = calc_err_interp(solcell, refpattern)
% Calculates errors wrt a reference solution determined by a pattern
% ONLY FOR FLYING SPAGHETTI WITH DIFFERENT n

ref = refpattern;
if ~ref.rslt.finished == 1
   error('calc_errors:refDidntFinish',...
         'Reference solution did not finish');
end

% Error functions
norm2 = @(x) sqrt(sum(x.^2,1)/size(x,1));
errlog = @(x,xref) [log_s3(qp(inv_s3sdr3(xref(1:4)),x(1:4))); T_inv(log_s3(qp(inv_s3sdr3(xref(1:4)),x(1:4)))).' * (x(5:7)-xref(5:7))];

% reference q=(p,x)
my_n_ref = length(ref.rslt.q(:,end))/7;

% Calculate errata
for i=1:numel(solcell)
    my_n = length(solcell{i}.rslt.q(:,end))/7;
    solcell{i}.err_interp.abs.q = norm2(errlog(solcell{i}.rslt.q(7*(my_n-1)/2:7*(my_n-1)/2+7,end),ref.rslt.q(7*(my_n_ref-1)/2:7*(my_n_ref-1)/2+7,end)));
end
end