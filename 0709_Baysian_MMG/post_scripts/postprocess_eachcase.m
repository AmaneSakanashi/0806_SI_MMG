function [] = postprocess_eachcase(setting_switch, usr_color,CmaSettings, case_attribute, CMAsetting_path, CMAprob_path,outputcmaes_path, pp)
    close all
    
    disp(case_attribute.name)
    fid_CMAsetting = fopen(CMAsetting_path);
    fid_CMAprob = fopen(CMAprob_path);
    fid_outputcmaes=fopen(outputcmaes_path); 
    directory_path = strcat(case_attribute.folder, '/', case_attribute.name,'/');
    
    %% make save directory for processed data
    save_directory_path = strcat(directory_path,'/','processed','/');
    status_mkdir = mkdir(save_directory_path);
    
    %% read outputCmaes.txt
    if setting_switch.read_outputcames == 1
        [CmaSettings, itr, minj, minarffeas, arffeas,arf,df,dx,sigma,vecD,xmean,...
            sigma_sqrt_eigC,sigma_sqrt_daigC,dim,itrplot,Xopt] ...
            = read_outputcmaes(fid_CMAsetting, fid_CMAprob, fid_outputcmaes, CmaSettings);
            savename_cmaes_itr = strcat(save_directory_path, 'cmaes_itr_result.mat');
            save(savename_cmaes_itr)
    else
        loadname_cmaes_itr = strcat(save_directory_path, 'cmaes_itr_result.mat');
        load(loadname_cmaes_itr,'CmaSettings', 'itr', 'minj', 'minarffeas', 'arffeas','arf',...
                                'df','dx','sigma','vecD','xmean',...
                                'sigma_sqrt_eigC', 'sigma_sqrt_daigC', 'dim','itrplot','Xopt')
    end
    %% make figure of cma-es iteration
    if setting_switch.draw_itr == true
        mkfigure_CMAES_itr(itr,arffeas,sigma_sqrt_eigC,sigma_sqrt_daigC,dim,save_directory_path)
    end
     %% update coefficients by cma-es results
    initial_params_csv_path = strcat(directory_path, '/MMG_params.csv');
    [hydro_derivative, mmgparams, mmgparams_ori, update_index, update_params]=update_coeff(dim, Xopt, CmaSettings, initial_params_csv_path);
    
    % save mmg params to csv 
    save_mmgparams_tocsv(save_directory_path, case_attribute.name, update_params, update_index, hydro_derivative,Xopt)
    
    % show the parameters which Xopt are 1<<Xopt(i)
    not_conv_index = Xopt > 1.5;
    if not_conv_index == false (size(Xopt))
        formatSpec = 'all parameter converged';
        fprintf(formatSpec)
    else
        formatSpec = 'not converved parameter:';
        fprintf( formatSpec)
        disp(hydro_derivative.parameter(not_conv_index))
    end

   %% compute obj function for train file
%     CmaSettings.inputdata = 1; % 1: train data, 2: test data
%     CmaSettings.validation_type = 1; %1: compute obj func, 2: compare trajectory and stateval, 3;compare trajectory and stateval(z and stop) 
%     start_time_step = 1; % timestep number to start MMG sim; dont change for train data
%     
%     [ Obj_train, Obj_ori_train,~,~, ~, ~,~, ~, ~,~,~,~,~, ~, ~,~  ] = ...
%     compute_state(directory_path, save_directory_path, case_attribute.name, CmaSetting, pp, mmgparams, mmgparams_ori, arffeas,minj,itr,start_time_step);
% 
%     varNames = {'Obj_train_Fortran','Obj_train_matlab','Obj_test_EFD_Coeff'};
%     minobjectives = minarffeas*ones(3,1);
%     savetable=table(minobjectives, Obj_train, Obj_ori_train,'variablenames',varNames); 
%     savename=strcat(save_directory_path, 'train_Objfunc.csv');
%     writetable(savetable,savename)
    %% compute obj function for test data
     CmaSettings.validation_type = 2; %1: compute obj func, 2: compare trajectory and stateval, 3;compare trajectory and stateval(z and stop)
     CmaSettings.inputdata = 2; % 1: train data, 2: test data
     start_time_step = 1; % timestep number to start MMG sim;

    compute_state(setting_switch, directory_path, save_directory_path, case_attribute.name,...
                  usr_color, CmaSettings, pp, mmgparams, mmgparams_ori, arffeas,minj,itr,start_time_step);

end