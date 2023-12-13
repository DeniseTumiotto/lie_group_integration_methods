function delete_all_config_and_bin(pattern, with_bin)
% Delete all config and bin that match a certain pattern

% If the with_bin argument is omitted, we assume we don't want to load the
% binary result files
if nargin < 2
   with_bin = false;
end

% If argument pattern is omitted, all files are loaded
if nargin == 0
   pattern = struct();
end

yn = input('Are you sure you want to delete these files? (y/n)','s');

if ~strcmp(yn,'y')
   display('Aborted');
   return;
end

path = '../out/';

files = dir([path '*.lua']);

for i=1:length(files)
   if strcmp(with_bin,'with_bin')
      [~, matched] = load_config_and_bin([path files(i).name(1:end-4)], pattern);
   else
      [~, matched] = load_config_and_bin([path files(i).name(1:end-4)], pattern, 'only_conf');
   end
   if matched
      delete([path files(i).name(1:end-4) '.*']);
   end
end