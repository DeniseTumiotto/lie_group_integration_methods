% useful_plots

txt_cell_array = {'snapshots';'step size control plots';'order of convergence plots';'animation';'CPU time'};

[indx,tf] = listdlg('ListString',txt_cell_array);

if tf == 0
    closing = questdlg('Do you want to leave?', ...
        'leaving',...
        'Yes, I want to leave', ...
        'No, it was a mistake',...
        'Yes, I want to leave');
    
    switch closing
        case 'Yes, I want to leave'
            disp('Thank you for using me! Have a nice day!')
            return
        case 'No, it was a mistake'
            useful_plots
    end
end

prompt = {'Which variable do I need to use?'};
dlgtitle = 'The solutions';
definput = {'s'};
dims = [1 40];
sols = inputdlg(prompt,dlgtitle,dims,definput);

for i = indx
    switch i
        case 1

        case 2
        case 3
        case 4
        case 5
    end
end


