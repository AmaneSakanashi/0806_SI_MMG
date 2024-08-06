function save_mmgparams_tocsv(save_dir_path, case_name, update_params, update_index, hydro_derivative,Xopt)
    varNames ={'parameter_name','EFD','optimized', 'optimized_EFD_ratio','Xpot_by_cma'};
    savename = strcat(save_dir_path, '/',case_name,'_optimized_coefficient.csv');
    ratio = zeros(length(update_params),1);
    for m = 1:length(update_params)
        ratio(m) =update_params(m) / hydro_derivative.params_init(update_index(m)) ;
    end
   savetable = table(hydro_derivative.parameter(update_index), hydro_derivative.params_init(update_index),...
   update_params, ratio, transpose(Xopt),'variablenames',varNames); 
    writetable(savetable,savename)
end