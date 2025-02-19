function single_result_to_file(sol, values, filename, sample)

    if nargin < 4
        sample = 1;
    end

    my_matrix = [];
    
    for i = 1:length(values)
        tmp_matrix = eval(['sol.' values{i}]);
        my_matrix(:,i) = tmp_matrix(1:sample:end);
    end
    my_cell = [values; num2cell(squeeze(my_matrix))];
    writecell(my_cell, filename, 'Delimiter', ' ')
end