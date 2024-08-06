function compute_state(setting_switch, directory_path, save_directory_path, condition_name,...
                    usr_color, CmaSettings, pp, mmgparams, mmgparams_ori, arffeas,minj,itr, start_time_step)
%% read input fils


if(CmaSettings.inputdata ==1)
    %%% use train data
    train_list_path = dir(strcat(directory_path,'/','*list.txt'));
%     files.
elseif(CmaSettings.inputdata ==2)
    %%%% use test data
    files = dir("../testdata/*.csv");
end

step_max = zeros(length(files),1);
for i=1:length(files)
    if(CmaSettings.inputdata ==1)
        %%% use train data
        xfileplace = ["../traindata/", files(i).name];
    elseif(CmaSettings.inputdata ==2)
        %%%% use test data
        xfileplace = ["../testdata/", files(i).name];
    end
       
    fileplace = xfileplace(1)+xfileplace(2);
    data = readtable(fileplace);
    step_max(i) = height(data);
end

time_input=zeros(max(step_max)-start_time_step+1,length(files));
x_input=zeros(max(step_max)-start_time_step+1,length(files));
u_input=zeros(max(step_max)-start_time_step+1,length(files));
y_input=zeros(max(step_max)-start_time_step+1,length(files));
vm_input=zeros(max(step_max)-start_time_step+1,length(files));
psi_input=zeros(max(step_max)-start_time_step+1,length(files));
r_input=zeros(max(step_max)-start_time_step+1,length(files));
apparent_wind_velo =zeros(max(step_max)-start_time_step+1, 1);
apparent_wind_dir=zeros(max(step_max)-start_time_step+1, 1);
Wind_Velocity =zeros(max(step_max)-start_time_step+1,length(files));
Wind_Direction =zeros(max(step_max)-start_time_step+1,length(files));
n_input =zeros(max(step_max)-start_time_step+1,length(files));
delta_input =zeros(max(step_max)-start_time_step+1,length(files));

for i=1:length(files)
    if(CmaSettings.inputdata ==1)
        %%% use train data
        xfileplace = ["../traindata/", files(i).name];
    elseif(CmaSettings.inputdata ==2)
        %%%% use test data
        xfileplace = ["../testdata/", files(i).name];
    end
    fileplace = xfileplace(1)+xfileplace(2);
    data = readtable(fileplace);
    for j =start_time_step:step_max(i)
        time_input(j-start_time_step+1,i) = data.t_s_(j);
        x_input(j-start_time_step+1,i) = data.x_position_mid_m_(j);
        u_input(j-start_time_step+1,i) =data.u_velo_m_s_(j);
        y_input(j-start_time_step+1,i) =data.y_position_mid_m_(j);
        vm_input(j-start_time_step+1,i) =data.vm_velo_m_s_(j);
        psi_input(j-start_time_step+1,i) =data.psi_hat_rad_(j);
        r_input(j-start_time_step+1,i) =data.r_angvelo_rad_s_(j);
        apparent_wind_velo(j-start_time_step+1) = data.wind_velo_relative_mid_m_s_(j);
        apparent_wind_dir(j-start_time_step+1) = data.wind_dir_relative_mid_rad_(j);
        Wind_Velocity(j-start_time_step+1,i) =data.wind_velo_true_m_s_(j);
        Wind_Direction(j-start_time_step+1,i) =data.wind_dir_true_rad_(j);
        n_input(j-start_time_step+1,i) =data.n_prop_rps_(j);
        delta_input(j-start_time_step+1,i) =data.delta_rudder_rad_(j);
    end
end

Obj=zeros(3,1);
Obj_ori = zeros(3,1);
%%%%%%%%%%%%%%%%%%%%%%%
%%%% compute Object function part %%%
%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(files)
    % reset obj for each test data
    if(CmaSettings.inputdata ==2)
        Obj=zeros(3,1);
        Obj_ori = zeros(3,1);
    end
    % initialize
    average_input_sin = 0.0d0;
    average_sim_sin = 0.0d0;
    average_sim_ori_sin = 0.0d0;
    
    average_input_cos = 0.0d0;
    average_sim_cos = 0.0d0;
    average_sim_ori_cos = 0.0d0;
    
    stdev_input_sin = 0.0d0;
    stdev_sim_sin = 0.0d0;
    stdev_sim_ori_sin = 0.0d0;
    
    stdev_input_cos = 0.0d0;
    stdev_sim_cos = 0.0d0;
    stdev_sim_ori_cos = 0.0d0;
    
    average_input = zeros(6,1);
    average_sim = zeros(6,1);
    average_sim_ori = zeros(6,1);
    stdev_input  = zeros(6,1);
    stdev_sim = zeros(6,1);
    stdev_sim_ori =zeros(6,1);
    timestep_len = length(time_input);
    standardrized_input = zeros(timestep_len,6);
    standardrized_sim  = zeros(timestep_len,6);
    standardrized_sim_ori = zeros(timestep_len,6);
    standardrized_input_sin = zeros(timestep_len,1);
    standardrized_sim_sin = zeros(timestep_len,1);
    standardrized_sim_ori_sin = zeros(timestep_len,1);
    standardrized_input_cos = zeros(timestep_len,1);
    standardrized_sim_cos = zeros(timestep_len,1);
    standardrized_sim_ori_cos = zeros(timestep_len,1);
    state_sim = zeros(timestep_len,6);
    state_sim_ori = zeros(timestep_len,6);
    state_instant  = zeros(6,1) ;
    
    % make state variables of train data    
        state_input(:,1) = x_input(:,i);
        state_input(:,2) = u_input(:,i);
        state_input(:,3) = y_input(:,i);
        state_input(:,4) = vm_input(:,i);
        state_input(:,5) = psi_input(:,i);
        state_input(:,6) = r_input(:,i);
   % define time step size
%         time_step_size = round((time_input(3,i)-time_input(2,i)),3);
        time_step_size = round((time_input(3,i)-time_input(2,i)),3);
    % define reset frequency
    if(CmaSettings.validation_type==1)
        reset_freq = round(CmaSettings.integration_period/time_step_size); 
    elseif(CmaSettings.validation_type==2)
        reset_freq = round(CmaSettings.validation_period/time_step_size);  
    elseif(CmaSettings.validation_type==3)
        reset_freq = 1;
    end
     for j = 1:step_max(i)-start_time_step
        if(CmaSettings.validation_type==1 || CmaSettings.validation_type==2 )
            if(mod(j, reset_freq)==1)
                state_sim(j,:) = state_input(j, :);
                state_sim_ori(j,:) = state_input(j, :);
            end
        elseif(CmaSettings.validation_type==3 )
            if j==1
                state_sim(j,:) = state_input(j, :);
                state_sim_ori(j,:) = state_input(j, :);
            end
        end
         % time stepping by EFD MMG params
         for k=1:6
            state_instant(k,1) = state_sim_ori(j,k);
         end
         
         k1 = Right_MMGLS(time_input(j,i),  state_instant,                                                  delta_input(j,i), n_input(j,i), pp, mmgparams_ori, CmaSettings.switch_wind, Wind_Direction(j,i), Wind_Velocity(j,i));
         k2 = Right_MMGLS(time_input(j,i),  state_instant + 0.5d0 * k1 * time_step_size, delta_input(j,i), n_input(j,i), pp, mmgparams_ori, CmaSettings.switch_wind, Wind_Direction(j,i), Wind_Velocity(j,i));
         k3 = Right_MMGLS(time_input(j,i),  state_instant + 0.5d0 * k2 * time_step_size, delta_input(j,i), n_input(j,i), pp, mmgparams_ori, CmaSettings.switch_wind, Wind_Direction(j,i), Wind_Velocity(j,i));
         k4 = Right_MMGLS(time_input(j,i),  state_instant +              k3 * time_step_size, delta_input(j,i), n_input(j,i), pp, mmgparams_ori, CmaSettings.switch_wind, Wind_Direction(j,i), Wind_Velocity(j,i));
                 
         right = ( k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) / 6.0d0;
         for k=1:6
            state_sim_ori(j+1,k) = state_instant(k) + right(k) * time_step_size; 
         end
         
         % time stepping by CMA MMG params
         for k=1:6
            state_instant(k,1) = state_sim(j,k);
         end
         
         k1 = Right_MMGLS(time_input(j,i),  state_instant,                                                  delta_input(j,i), n_input(j,i), pp, mmgparams, CmaSettings.switch_wind, Wind_Direction(j,i), Wind_Velocity(j,i));
         k2 = Right_MMGLS(time_input(j,i),  state_instant + 0.5d0 * k1 * time_step_size, delta_input(j,i), n_input(j,i), pp, mmgparams, CmaSettings.switch_wind, Wind_Direction(j,i), Wind_Velocity(j,i));
         k3 = Right_MMGLS(time_input(j,i),  state_instant + 0.5d0 * k2 * time_step_size, delta_input(j,i), n_input(j,i), pp, mmgparams, CmaSettings.switch_wind, Wind_Direction(j,i), Wind_Velocity(j,i));
         k4 = Right_MMGLS(time_input(j,i),  state_instant +              k3 * time_step_size, delta_input(j,i), n_input(j,i), pp, mmgparams, CmaSettings.switch_wind, Wind_Direction(j,i), Wind_Velocity(j,i));
                 
         right = ( k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) / 6.0d0;
%                     !!! limit max velocity and acceleration 
            grav = 9.80665;
            acc_limit = grav*1.0d10;
            rdot_limit = acc_limit/(pp.lpp/2.0d0);
            velo_limit= 1.0d10;
            angvelo_limit= velo_limit/(pp.lpp/2.0d0);
            for k =1:6
                if (k==1.||k==3)
                    if(abs(right(k))>velo_limit)
                        comment = "limit of velocity"
                        pause
                    end
                elseif (k==2 || k==4) 
                    if (abs(right(k)) >acc_limit)
                        comment = "limit of acc"
                        pause
                    end 
               elseif(k==5)
                    if(abs(right(k))>angvelo_limit)
                        comment = "limit of ang velocity"
                        pause
                    end
               elseif(k == 6) 
                    if (abs(right(k)) > rdot_limit)
                        comment = "limit of ang acc"
                        pause
                    end
                end 
            end
            
         for k=1:6
            state_sim(j+1,k) = state_instant(k) + right(k) * time_step_size; 
         end
         

         for k =1:6
             if(k==5)
                average_input_sin = average_input_sin +sin(state_input(j,k));
                average_sim_sin = average_sim_sin +sin(state_sim(j,k));
                average_sim_ori_sin = average_sim_ori_sin +sin(state_sim_ori(j,k));
                
                average_input_cos = average_input_cos +cos(state_input(j,k));
                average_sim_cos = average_sim_cos +cos(state_sim(j,k));
                average_sim_ori_cos = average_sim_ori_cos +cos(state_sim_ori(j,k));
             else
                average_input(k)= average_input(k)+state_input(j,k);
                average_sim(k) = average_sim(k) + state_sim(j,k);
                average_sim_ori(k) = average_sim_ori(k) + state_sim_ori(j,k);
             end
         end   
     end
     
     %%% compute object func
     for k =1:6
         if k==5
           average_input_sin = average_input_sin/(step_max(i)-start_time_step+1);
           average_sim_sin = average_sim_sin/(step_max(i)-start_time_step+1);
           average_sim_ori_sin = average_sim_ori_sin/(step_max(i)-start_time_step+1);
           
           average_input_cos = average_input_cos/(step_max(i)-start_time_step+1);
           average_sim_cos = average_sim_cos/(step_max(i)-start_time_step+1);
           average_sim_ori_cos = average_sim_ori_cos/(step_max(i)-start_time_step+1);
           
           %! standard deviation of val
           for j = 1:step_max(i)-start_time_step+1
                    stdev_input_sin = stdev_input_sin + (sin(state_input(j,k))-average_input_sin)^2.0d0;
                    stdev_sim_sin   = stdev_sim_sin   + (sin(state_sim(j,k))  -average_sim_sin)^2.0d0;
                    stdev_sim_ori_sin   = stdev_sim_ori_sin   + (sin(state_sim_ori(j,k))  -average_sim_ori_sin)^2.0d0;
                    
                    stdev_input_cos = stdev_input_cos + (cos(state_input(j,k))-average_input_cos)^2.0d0;
                    stdev_sim_cos   = stdev_sim_cos   + (cos(state_sim(j,k))  -average_sim_cos)^2.0d0;
                    stdev_sim_ori_cos   = stdev_sim_ori_cos   + (cos(state_sim_ori(j,k))  -average_sim_ori_cos)^2.0d0;
           end 
           stdev_input_sin = sqrt(stdev_input_sin/(step_max(i)-start_time_step+1));
           stdev_sim_sin = sqrt(stdev_sim_sin/(step_max(i)-start_time_step+1));
           stdev_sim_ori_sin = sqrt(stdev_sim_ori_sin/(step_max(i)-start_time_step+1));
           
           stdev_input_cos = sqrt(stdev_input_cos/(step_max(i)-start_time_step+1));
           stdev_sim_cos = sqrt(stdev_sim_cos/(step_max(i)-start_time_step+1));
           stdev_sim_ori_cos = sqrt(stdev_sim_ori_cos/(step_max(i)-start_time_step+1));
           
           %!!! standarization
           for j = 1:step_max(i)-start_time_step+1
                    standardrized_input_sin(j) = (sin(state_input(j,k))-average_input_sin)/stdev_input_sin;
                    standardrized_sim_sin(j) = (sin(state_sim(j,k))-average_sim_sin)/stdev_sim_sin;
                    standardrized_sim_ori_sin(j) = (sin(state_sim_ori(j,k))-average_sim_ori_sin)/stdev_sim_ori_sin;
                    
                    standardrized_input_cos(j) = (cos(state_input(j,k))-average_input_cos)/stdev_input_cos;
                    standardrized_sim_cos(j) = (cos(state_sim(j,k))-average_sim_cos)/stdev_sim_cos;
                    standardrized_sim_ori_cos(j) = (cos(state_sim_ori(j,k))-average_sim_ori_cos)/stdev_sim_ori_cos;
           end  
         else
            % average of state varibale
            average_input(k)= average_input(k)/(step_max(i)-start_time_step+1);
            average_sim(k) = average_sim(k)/(step_max(i)-start_time_step+1);
            average_sim_ori(k) = average_sim_ori(k)/(step_max(i)-start_time_step+1);
            
            % standard deviation of variable
            for j = 1:step_max(i)-start_time_step+1
                stdev_input(k) = stdev_input(k) + (state_input(j,k)-average_input(k))^2.0d0;
                stdev_sim(k) = stdev_sim(k) + (state_sim(j,k)-average_sim(k))^2.0d0;
                stdev_sim_ori(k) = stdev_sim_ori(k) + (state_sim_ori(j,k)-average_sim_ori(k))^2.0d0;
            end
            stdev_input(k) = sqrt(stdev_input(k)/(step_max(i)-start_time_step+1));
            stdev_sim(k) = sqrt(stdev_sim(k)/(step_max(i)-start_time_step+1));
            stdev_sim_ori(k) = sqrt(stdev_sim_ori(k)/(step_max(i)-start_time_step+1));
             
            %%% standarization
            for j = 1:step_max(i)-start_time_step+1
                standardrized_input(j,k) = (state_input(j,k)-average_input(k))/stdev_input(k);
                standardrized_sim(j,k) = (state_sim(j,k)-average_sim(k))/stdev_sim(k);
                standardrized_sim_ori(j,k) = (state_sim_ori(j,k)-average_sim_ori(k))/stdev_sim_ori(k);
            end
         end
     end
      

     for j = 1:step_max(i)-start_time_step+1
              % compute obj func form standardrized u, vm and r  
            Obj(1) = Obj(1) + ((standardrized_input(j,2) - standardrized_sim(j,2)) ^2.0d0 +...
                        (standardrized_input(j,4) - standardrized_sim(j,4)) ^2.0d0 +...
                        (standardrized_input(j,6) - standardrized_sim(j,6)) ^2.0d0 ) * time_step_size;
            Obj_ori(1) = Obj_ori(1) + ((standardrized_input(j,2) - standardrized_sim_ori(j,2)) ^2.0d0 +...
                        (standardrized_input(j,4) - standardrized_sim_ori(j,4)) ^2.0d0 +...
                        (standardrized_input(j,6) - standardrized_sim_ori(j,6)) ^2.0d0 ) * time_step_size;
             % compute obj func form standardrized X, Y, psi               
            Obj(3) = Obj(3) + ((standardrized_input(j,1) - standardrized_sim(j,1)) ^2.0d0 +...
                        (standardrized_input(j,3) - standardrized_sim(j,3)) ^2.0d0 +...
                        (standardrized_input_sin(j) - standardrized_sim_sin(j)) ^2.0d0 +...
                        (standardrized_input_cos(j) - standardrized_sim_cos(j)) ^2.0d0 ) * time_step_size;
            Obj_ori(3) = Obj_ori(3) + ((standardrized_input(j,1) - standardrized_sim_ori(j,1)) ^2.0d0 +...
                        (standardrized_input(j,3) - standardrized_sim_ori(j,3)) ^2.0d0 +...
                        (standardrized_input_sin(j) - standardrized_sim_ori_sin(j)) ^2.0d0 +...
                        (standardrized_input_cos(j) - standardrized_sim_ori_cos(j)) ^2.0d0) * time_step_size;    
     end
          % compute obj func form standardrized all six componets 
     for k =1:6
         for j = 1:step_max(i)-start_time_step+1
             if(k==5)
                  Obj(2) = Obj(2) +( (standardrized_input_sin(j) - standardrized_sim_sin(j)) ^2.0d0  ...
                            +  (standardrized_input_cos(j) - standardrized_sim_cos(j)) ^2.0d0)* time_step_size;
                  Obj_ori(2) = Obj_ori(2) +( (standardrized_input_sin(j) - standardrized_sim_ori_sin(j)) ^2.0d0  ....
                            +  (standardrized_input_cos(j) - standardrized_sim_ori_cos(j)) ^2.0d0)* time_step_size;          
             else
                  Obj(2) = Obj(2) + (standardrized_input(j,k) - standardrized_sim(j,k)) ^2.0d0  * time_step_size;
                  Obj_ori(2) = Obj_ori(2) + (standardrized_input(j,k) - standardrized_sim_ori(j,k)) ^2.0d0  * time_step_size;                  
             end
         end
     end
    %%%% compute Object function for each mini-batch 
%     if CmaSettings.validation_type ==2
%         number_batch = fix(step_max( CmaSettings.number_testfiles)/reset_freq);
%     
%         Obj_batch=zeros(3, number_batch);
%         Obj_ori_batch = zeros(3, number_batch);
% 
%         for batch = 1:number_batch
%              lot = [1+reset_freq*(batch-1):reset_freq*batch];
%              for j = lot(1):lot(end)
%                 % compute obj func form standardrized u, vm and r  
%                 Obj_batch(1,batch) = Obj_batch(1,batch) + ((standardrized_input(j,2) - standardrized_sim(j,2)) ^2.0d0 +...
%                             (standardrized_input(j,4) - standardrized_sim(j,4)) ^2.0d0 +...
%                             (standardrized_input(j,6) - standardrized_sim(j,6)) ^2.0d0 ) * time_step_size;
%                 Obj_ori_batch(1,batch) = Obj_ori_batch(1,batch) + ((standardrized_input(j,2) - standardrized_sim_ori(j,2)) ^2.0d0 +...
%                             (standardrized_input(j,4) - standardrized_sim_ori(j,4)) ^2.0d0 +...
%                             (standardrized_input(j,6) - standardrized_sim_ori(j,6)) ^2.0d0 ) * time_step_size;
% 
%                 % compute obj func form standardrized X, Y, psi               
%                 Obj_batch(3,batch) = Obj_batch(3,batch) + ((standardrized_input(j,1) - standardrized_sim(j,1)) ^2.0d0 +...
%                             (standardrized_input(j,3) - standardrized_sim(j,3)) ^2.0d0 +...
%                             (standardrized_input_sin(j) - standardrized_sim_sin(j)) ^2.0d0 +...
%                             (standardrized_input_cos(j) - standardrized_sim_cos(j)) ^2.0d0 ) * time_step_size;
%                 Obj_ori_batch(3,batch) = Obj_ori_batch(3,batch) + ((standardrized_input(j,1) - standardrized_sim_ori(j,1)) ^2.0d0 +...
%                             (standardrized_input(j,3) - standardrized_sim_ori(j,3)) ^2.0d0 +...
%                             (standardrized_input_sin(j) - standardrized_sim_ori_sin(j)) ^2.0d0 +...
%                             (standardrized_input_cos(j) - standardrized_sim_ori_cos(j)) ^2.0d0) * time_step_size;    
%             end
%               % compute obj func form standardrized all six componets 
%             for k =1:6
%              for j = lot(1):lot(end)
%                  if(k==5)
%                       Obj_batch(2,batch) = Obj_batch(2,batch) +( (standardrized_input_sin(j) - standardrized_sim_sin(j)) ^2.0d0  ...
%                                 +  (standardrized_input_cos(j) - standardrized_sim_cos(j)) ^2.0d0)* time_step_size;
%                       Obj_ori_batch(2,batch) = Obj_ori_batch(2,batch) +( (standardrized_input_sin(j) - standardrized_sim_ori_sin(j)) ^2.0d0  ....
%                                 +  (standardrized_input_cos(j) - standardrized_sim_ori_cos(j)) ^2.0d0)* time_step_size;          
%                  else
%                       Obj_batch(2,batch) = Obj_batch(2,batch) + (standardrized_input(j,k) - standardrized_sim(j,k)) ^2.0d0  * time_step_size;
%                       Obj_ori_batch(2,batch) = Obj_ori_batch(2,batch)  + (standardrized_input(j,k) - standardrized_sim_ori(j,k)) ^2.0d0  * time_step_size;                  
%                  end
%              end
%             end
%         end
%     else
        Obj_batch=0.0;
        Obj_ori_batch=0.0;
%     end
    if (CmaSettings.validation_type ==2 || CmaSettings.validaton_type ==3)
        formatSpec = 'Obj (Matlab)  is %e %e %e at %i itr\n';
        fprintf( formatSpec, Obj(1),Obj(2),Obj(3), itr(minj));
%         formatSpec = 'Obj (Matlab, EFD Coefs.)  is %e %e %e';
%         fprintf( formatSpec, Obj_ori(1), Obj_ori(2), Obj_ori(3));
        fprintf(newline)
        
        varNames = {'Obj_test_CMAES','Obj_test_EFD'};
        savetable=table(Obj, Obj_ori,'variablenames',varNames); 
        savename=strcat(save_directory_path,condition_name,'_',files(i).name,'_DTeval=',num2str(CmaSettings.validation_period), '_Objfunc.csv');
 
        writetable(savetable,savename)
    end
    if setting_switch.draw_traj_fig == true
         draw_traj_figure(save_directory_path, condition_name, files(i).name, start_time_step, usr_color, CmaSettings,...
                          state_input(1:step_max(i),:), state_sim(1:step_max(i),:), state_sim_ori(1:step_max(i),:),...
                          n_input(1:step_max(i),i), delta_input(1:step_max(i),i), ...
                          step_max(i), apparent_wind_velo(1:step_max(i)), apparent_wind_dir(1:step_max(i)),...
                          time_input(1:step_max(i),i),pp )
    end
end



if(CmaSettings.validation_type==1)
    formatSpec = 'Obj (Fortran) is %e at %i itr\n';
    fprintf( formatSpec, arffeas( minj ), itr(minj))
    formatSpec = 'Obj (Matlab)  is %e %e %e at %i itr\n';
    fprintf( formatSpec, Obj(1),Obj(2),Obj(3), itr(minj));
%     formatSpec = 'Obj (Matlab, EFD Coefs.)  is %e %e %e';
%     fprintf( formatSpec, Obj_ori(1), Obj_ori(2), Obj_ori(3));
    fprintf(newline)

end


end