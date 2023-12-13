function delete_duplicate_config_and_bin(pattern)
% Load all config files and remove duplicates

% If argument pattern is omitted, all files are loaded
if nargin == 0
   pattern = struct();
end

path = '../out/';

files = dir([path '*.lua']);

[~, sorti] = sort(cat(1,files.name),'descend');

files = files(sorti);

conflist = [];

for i=1:length(files)
   [sol, matched] = load_config_and_bin([path files(i).name(1:end-4)], ...
                                        pattern, ...
                                        'only_conf');
   if matched
      for i=1:length(conflist)
         if structmatch(sol,conflist(i)) && structmatch(conflist(i),sol)
            % Structs sol and conflist(i) are equal in this case
            delete([path files(i).name(1:end-4) '.*']);
         end
      end
      conflist = [conflist sol];
   end
end