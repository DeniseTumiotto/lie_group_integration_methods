function result_to_file(sols, xvalue, yvalues, filename, upto)

the_size = max(size(sols));
if nargin < 5
    upto = the_size;
end

my_matrix = zeros(upto,length(yvalues)+1);

for i = 1:upto
    for j = 1:length(yvalues)+1
        if j == 1
            my_matrix(i,j) = sols{i}.(xvalue);
        else
            my_matrix(i,j) = eval(['sols{i}.' yvalues{j-1}]);
        end
    end
end

my_cell = [[{xvalue}, yvalues]; num2cell(my_matrix)];

writecell(my_cell,filename,'Delimiter',' ')

end