function sc = filter_solcell(solcell,pattern)
if ~isempty(pattern)
   sc = cell(1,0);
   for k=1:length(solcell)
      if structmatch(solcell{k}, pattern)
         sc{end+1} = solcell{k};
      end
   end
else
   sc = solcell;
end