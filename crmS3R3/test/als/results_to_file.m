function results_to_file(sols, xvalue, yvalues, filename, my_zeros)

    the_size = max(size(sols));

    my_matrix = [];
    if nargin < 5
        my_zeros = NaN;
    end
    
    for i = 1:the_size
        for j = 1:length(yvalues)+1
            if j == 1
                my_matrix(i,j) = sols{i}.(xvalue);
            else
                if eval(['sols{i}.' yvalues{j-1}]) < eps
                    my_matrix(i,j) = my_zeros;
                else
                    my_matrix(i,j) = eval(['sols{i}.' yvalues{j-1}]);
                end
            end
        end
    end
    
    my_cell = [[{xvalue}, yvalues]; num2cell(squeeze(my_matrix))];
    writecell(my_cell, filename, 'Delimiter', ' ')
end