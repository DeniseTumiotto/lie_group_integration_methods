function out = load_latest_config_and_bin(which)
% Load (latest) config and bin
%
% load_latest_config_and_bin(0) is equivalent to
% load_latest_config_and_bin and loads the latest config and bin files.
%
% If which = -1 it would load the next to latest files and so on.
%
% For which = 1 it would load the oldest files, which = 2 the second
% oldest and so on.
%

if nargin == 0
   which = 0;
end

path = '../out/';

files = dir([path '*.lua']);

[~, sorti] = sort(cat(1,files.datenum));


if which <= 0
   out = load_config_and_bin([path files(sorti(end+which)).name(1:end-4)]);
else
   out = load_config_and_bin([path files(sorti(which)).name(1:end-4)]);
end