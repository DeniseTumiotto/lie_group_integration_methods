function results_to_file_by_order(sols, xvalue, yvalues, filename)

    the_size = max(size(sols));
    
    my_orders = zeros(the_size,1);
    cnt_orders = ones(4,1);
    my_matrix = [];
    
    for i = 1:the_size
        my_orders(i) = sols{i}.order-1;
        for j = 1:length(yvalues)+1
            if j == 1
                my_matrix(sols{i}.order-1,cnt_orders(my_orders(i)),j) = sols{i}.(xvalue);
            else
                if eval(['sols{i}.' yvalues{j-1}]) < eps
                    my_matrix(sols{i}.order-1,cnt_orders(my_orders(i)),j) = NaN;
                else
                    my_matrix(sols{i}.order-1,cnt_orders(my_orders(i)),j) = eval(['sols{i}.' yvalues{j-1}]);
                end
            end
        end
        cnt_orders(my_orders(i)) = cnt_orders(my_orders(i))+1;
    end
    
    for i = 1:4
        if any(my_orders==i)
            my_cell = [[{xvalue}, yvalues]; num2cell(squeeze(my_matrix(i,:,:)))];
            writecell(my_cell, [filename,'_order_',num2str(i+1)], 'Delimiter', ' ')
        end
    end
    
    
    
end