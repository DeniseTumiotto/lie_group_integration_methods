% useful_plots

action_list = {'snapshots';'step size control plots';'order of convergence plots';'animation';'CPU time'};

[indx,tf] = listdlg('ListString',action_list);

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

if ~exist('s','var')
    prompt = {'Number of solutions'};
    dlgtitle = 'Sols';
    definput = {'all'};
    dims = [1 40];
    T = inputdlg(prompt,dlgtitle,dims,definput);
    if strcmp(T{1},'all')
        s = load_all_config_and_bin();
    else
        T = str2double(T{1});
        s = cell(1,T);
        for i = 0:-1:-T+1
            s{-i+1} = load_latest_config_and_bin(i);
        end
    end
end

for i = indx
    switch i
        case 1
            parameter_list = {'Plot extreme trajectory';'Plot in black&white';'Choose the number of shots';'none of the above'};
            [indx_parameter,tf] = listdlg('ListString',parameter_list,'InitialValue',4);
            traj_extr = 0; bw = 0; stride = 8;
            if tf ~= 0 && any(indx_parameter ~= 4)
                for k = indx_parameter
                    switch k
                        case 1
                            traj_extr = 1;
                        case 2
                            bw = 1;
                        case 3
                            prompt = {'Number of shots'};
                            dlgtitle = 'Choose a number';
                            definput = {'8'};
                            dims = [1 40];
                            stride = inputdlg(prompt,dlgtitle,dims,definput);
                            stride = uint8(str2num(stride{1}));
                        otherwise
                            continue
                    end
                end

            end
            for j = 1:T
                snapshotplot(s{j},'extreme',traj_extr,'stride',stride,'bw',bw)
            end
        case 2
            parameter_list = {'Choose the property';'Choose end time';'Error and time step size in different plot';...
                'Solutions in different plot';'none of the above'};
            [indx_parameter,tf] = listdlg('ListString',parameter_list,'InitialValue',5);
            the_name = 'problem_name'; t_end = 1e7; separate = 0; combine = 1;
            if tf ~= 0 && any(indx_parameter ~= 5)
                for k = indx_parameter
                    switch k
                        case 1
                            prompt = {'Name of property'};
                            dlgtitle = 'Choose a property';
                            definput = {'problem_name'};
                            dims = [1 40];
                            the_name = inputdlg(prompt,dlgtitle,dims,definput);
                        case 2
                            prompt = {'Select the end time'};
                            dlgtitle = 'End time';
                            definput = {'1'};
                            dims = [1 40];
                            t_end = inputdlg(prompt,dlgtitle,dims,definput);
                            t_end = str2num(t_end{1});
                        case 3
                            separate = 1;
                        case 4
                            combine = 0;
                        otherwise
                            continue
                    end
                end

            end           
            error_control_plot(s,'name',the_name,'end time',t_end,'separate',separate,'combine',combine)
        case 3
            if ~isfield(s{1},'err')
                s = mixed_calc_end_errors(s);
            end
        case 4
            parameter_list = {'Save a video';'Plot extreme trajectory';'Choose the step size';'none of the above'};
            [indx_parameter,tf] = listdlg('ListString',parameter_list,'InitialValue',4);
            make_video = 0; do_edges = 0; ts = [];
            if tf ~= 0 && any(indx_parameter ~= 4)
                for k = indx_parameter
                    switch k
                        case 1
                            make_video = 1;
                        case 2
                            do_edges = 1;
                        case 3
                            prompt = {'Indicate the step size'};
                            dlgtitle = 'Choose a number';
                            definput = {'2'};
                            dims = [1 40];
                            ts = inputdlg(prompt,dlgtitle,dims,definput);
                            ts = uint8(str2num(ts{1}));
                        otherwise
                            continue
                    end
                end

            end
            for j = 1:T            
                visualize_beam(s{j},'video',make_video,'edges',do_edges,'step',ts)
            end
        case 5
            parameter_list = {'Choose the error variable';'Combine over the problems';'none of the above'};
            [indx_parameter,tf] = listdlg('ListString',parameter_list,'InitialValue',3);
            my_var = 'q'; combine = 'no';
            if tf ~= 0 && any(indx_parameter ~= 5)
                for k = indx_parameter
                    switch k
                        case 1
                            parameter_list = {'q';'v';'l'};
                            [indx_parameter,tf] = listdlg('ListString',parameter_list,'InitialValue',1);
                            if tf == 0
                                warning('No selection, default value will be used!')
                            else
                                if length(indx_parameter) == 3
                                    my_var = 'all';
                                elseif length(indx_parameter) == 2
                                    my_var = string(parameter_list(indx_parameter(1)))+' '+string(parameter_list(indx_parameter(2)));
                                else
                                    my_var = string(parameter_list(indx_parameter));
                                end
                            end
                        case 2
                            combine = 'problem';
                        otherwise
                            continue
                    end
                end

            end
            if ~isfield(s{1},'err')
                s = mixed_calc_end_errors(s);
            end
            cpu_time_plot(s,'error',my_var,'combine',combine)
    end
end

clearvars -except s
