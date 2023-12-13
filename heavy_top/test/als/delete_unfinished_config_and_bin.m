function delete_unfinished_config_and_bin(pattern)
% Delete all unfinished config and bin files

pattern.rslt.finished = false;

delete_all_config_and_bin(pattern, 'with_bin');