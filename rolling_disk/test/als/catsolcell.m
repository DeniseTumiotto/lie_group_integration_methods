function out = catsolcell(solcell, fname)

out = zeros(1,numel(solcell));
for i=1:numel(solcell)
   try
      out(i) = eval(['solcell{i}.' fname]);
   catch ME
      warning(['Error occured during call of ''solcell{' ...
         num2str(i) '}.' num2str(fname) 10 ...
         'Error message: ' ME.message 10 ...
         'giving back NaN instead']);
      out(i) = NaN;
   end
end