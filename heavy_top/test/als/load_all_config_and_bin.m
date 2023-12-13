function out = load_all_config_and_bin(pattern)
% Load all config and bin that match a certain pattern

% If argument pattern is omitted, all files are loaded
if nargin == 0
   pattern = struct();
end

path = '../out/';

files = dir([path '*.lua']);

ind = 1;

for i=1:length(files)
   [sol, matched] = load_config_and_bin([path files(i).name(1:end-4)], pattern);
   if matched
      out{ind} = sol;
      ind = ind + 1;
   end
end