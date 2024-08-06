function [hydro_derivative, mmgparams, mmgparams_ori, update_index, update_params ] = update_coeff(dim, Xopt, CmaSetting,initial_param_csv_path)
% Numerical simulation and display 

% Numerical simulation setting.
%
hydro_derivative = readtable(initial_param_csv_path);

    mmgparams.massx_nd = hydro_derivative.params_init(1);
    mmgparams.massy_nd= hydro_derivative.params_init(2);
    mmgparams.IzzJzz_nd= hydro_derivative.params_init(3);
    mmgparams.xuu_nd = hydro_derivative.params_init(4);
    mmgparams.Xvr_nd = hydro_derivative.params_init(5);
    mmgparams.Yv_nd = hydro_derivative.params_init(6);
    mmgparams.Yr_nd = hydro_derivative.params_init(7);
    mmgparams.Nv_nd = hydro_derivative.params_init(8);
    mmgparams.Nr_nd= hydro_derivative.params_init(9);
    mmgparams.CD = hydro_derivative.params_init(10);
    mmgparams.C_rY = hydro_derivative.params_init(11);
    mmgparams.C_rN= hydro_derivative.params_init(12);
    mmgparams.X_0F_nd = hydro_derivative.params_init(13);
    mmgparams.X_0A_nd= hydro_derivative.params_init(14);
    mmgparams.kt_coeff0 = hydro_derivative.params_init(15);
    mmgparams.kt_coeff1 = hydro_derivative.params_init(16);
    mmgparams.kt_coeff2 = hydro_derivative.params_init(17);
    mmgparams.t_prop = hydro_derivative.params_init(18);
    mmgparams.wP0 = hydro_derivative.params_init(19);
    mmgparams.tau_prop = hydro_derivative.params_init(20);
    mmgparams.CP_nd = hydro_derivative.params_init(21);
    mmgparams.xP_nd = hydro_derivative.params_init(22);
    mmgparams.A1 = hydro_derivative.params_init(23);
    mmgparams.A2 = hydro_derivative.params_init(24);
    mmgparams.A3 = hydro_derivative.params_init(25);
    mmgparams.A4 = hydro_derivative.params_init(26);
    mmgparams.A5 = hydro_derivative.params_init(27);
    mmgparams.A6 = hydro_derivative.params_init(28);
    mmgparams.A7 = hydro_derivative.params_init(29);
    mmgparams.A8= hydro_derivative.params_init(30);
    mmgparams.B1 = hydro_derivative.params_init(31);
    mmgparams.B2 = hydro_derivative.params_init(32);
    mmgparams.B3 = hydro_derivative.params_init(33);
    mmgparams.B4 = hydro_derivative.params_init(34);
    mmgparams.B5 = hydro_derivative.params_init(35);
    mmgparams.B6 = hydro_derivative.params_init(36);
    mmgparams.B7 = hydro_derivative.params_init(37);
    mmgparams.B8= hydro_derivative.params_init(38);
    mmgparams.C3 = hydro_derivative.params_init(39);
    mmgparams.C6 = hydro_derivative.params_init(40);
    mmgparams.C7 = hydro_derivative.params_init(41);
    mmgparams.C10= hydro_derivative.params_init(42);
    mmgparams.Jmin = hydro_derivative.params_init(43);
    mmgparams.alpha_prop= hydro_derivative.params_init(44);
    mmgparams.tR = hydro_derivative.params_init(45);
    mmgparams.aH_rudder = hydro_derivative.params_init(46);
    mmgparams.xh_rudder_nd = hydro_derivative.params_init(47);
    mmgparams.lR_nd = hydro_derivative.params_init(48);
    mmgparams.kx_rudder = hydro_derivative.params_init(49);
    mmgparams.kx_rudder_reverse = hydro_derivative.params_init(50);
    mmgparams.epsilon_rudder = hydro_derivative.params_init(51);
    mmgparams.cpr_rudder = hydro_derivative.params_init(52);
    mmgparams.gammaN = hydro_derivative.params_init(53);
    mmgparams.gammaP= hydro_derivative.params_init(54);
    mmgparams.XX0 = hydro_derivative.params_init(55);
    mmgparams.XX1 = hydro_derivative.params_init(56);
    mmgparams.XX3 = hydro_derivative.params_init(57);
    mmgparams.XX5 = hydro_derivative.params_init(58);
    mmgparams.YY1 = hydro_derivative.params_init(59);
    mmgparams.YY3 = hydro_derivative.params_init(60);
    mmgparams.YY5 = hydro_derivative.params_init(61);
    mmgparams.NN1 = hydro_derivative.params_init(62);
    mmgparams.NN2 = hydro_derivative.params_init(63);
    mmgparams.NN3 = hydro_derivative.params_init(64);
    mmgparams_ori = mmgparams;
    
    if(CmaSetting.type_si==1)
    % update all valiable
    % Ndim == 44
    % choose values to be updated
    update_index =[5, 6, 7, 8, 9, 10, 11, 12, ...
                   14, ...
                   18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, ...
                   45, 46, 47, 48, 49, 50, 51, 52, 53, 54] ;   
    update_params =zeros(dim,1);
    % update mmg params(hydro derivatives) by CMA-ES
    for i = 1:dim
        alt = update_index(i);
        update_params(i) = hydro_derivative.params_min(alt) +(hydro_derivative.params_max(alt) - hydro_derivative.params_min(alt)) * ( Xopt(i)+ 1.0d0 ) * 0.5d0;
    end 
    
    mmgparams.Xvr_nd = update_params(1);
    mmgparams.Yv_nd = update_params(2);
    mmgparams.Yr_nd = update_params(3);
    mmgparams.Nv_nd = update_params(4);
    mmgparams.Nr_nd= update_params(5);
    mmgparams.CD = update_params(6);
    mmgparams.C_rY = update_params(7);
    mmgparams.C_rN= update_params(8);
    mmgparams.X_0A_nd= update_params(9);
    mmgparams.t_prop = update_params(10);
    mmgparams.wP0 = update_params(11);
    mmgparams.tau_prop = update_params(12);
    mmgparams.CP_nd = update_params(13);
    mmgparams.xP_nd = update_params(14);
    mmgparams.A1 = update_params(15);
    mmgparams.A2 = update_params(16);
    mmgparams.A3 = update_params(17);
    mmgparams.A4 = update_params(18);
    mmgparams.A5 = update_params(19);
    mmgparams.A6 = update_params(20);
    mmgparams.A7 = update_params(21);
    mmgparams.A8= update_params(22);
    mmgparams.B1 = update_params(23);
    mmgparams.B2 = update_params(24);
    mmgparams.B3 = update_params(25);
    mmgparams.B4 = update_params(26);
    mmgparams.B5 = update_params(27);
    mmgparams.B6 = update_params(28);
    mmgparams.B7 = update_params(29);
    mmgparams.B8= update_params(30);
    mmgparams.C3 = update_params(31);
    mmgparams.C6 = update_params(32);
    mmgparams.C7 = update_params(33);
    mmgparams.C10= update_params(34);
    mmgparams.tR = update_params(35);
    mmgparams.aH_rudder = update_params(36);
    mmgparams.xh_rudder_nd = update_params(37);
    mmgparams.lR_nd = update_params(38);
    mmgparams.kx_rudder = update_params(39);
    mmgparams.kx_rudder_reverse = update_params(40);
    mmgparams.epsilon_rudder = update_params(41);
    mmgparams.cpr_rudder = update_params(42);
    mmgparams.gammaN = update_params(43);
    mmgparams.gammaP= update_params(44);

    elseif(CmaSetting.type_si==2)
    % update variables related to lowspeed and backward 
    % Ndim == 31
    update_index =[5, 6, 7, 8, 9, 10, 11, 12, ...
                   14, ...
                   23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, ...
                   50, 52];    
    
    % update mmg params(hydro derivatives) by CMA-ES
     update_params =zeros(dim,1);
    for i = 1:dim
        alt = update_index(i);
        update_params(i) = hydro_derivative.params_min(alt) +(hydro_derivative.params_max(alt) - hydro_derivative.params_min(alt)) * ( Xopt(i)+ 1.0d0 ) * 0.5d0;
    end 
        mmgparams.Xvr_nd = update_params(1);
        mmgparams.Yv_nd = update_params(2);
        mmgparams.Yr_nd = update_params(3);
        mmgparams.Nv_nd = update_params(4);
        mmgparams.Nr_nd= update_params(5);
        mmgparams.CD = update_params(6);
        mmgparams.C_rY = update_params(7);
        mmgparams.C_rN= update_params(8);
        mmgparams.X_0A_nd= update_params(9);
        mmgparams.A1 = update_params(10);
        mmgparams.A2 = update_params(11);
        mmgparams.A3 = update_params(12);
        mmgparams.A4 = update_params(13);
        mmgparams.A5 = update_params(14);
        mmgparams.A6 = update_params(15);
        mmgparams.A7 = update_params(16);
        mmgparams.A8= update_params(17);
        mmgparams.B1 = update_params(18);
        mmgparams.B2 = update_params(19);
        mmgparams.B3 = update_params(20);
        mmgparams.B4 = update_params(21);
        mmgparams.B5 = update_params(22);
        mmgparams.B6 = update_params(23);
        mmgparams.B7 = update_params(24);
        mmgparams.B8= update_params(25);
        mmgparams.C3 = update_params(26);
        mmgparams.C6 = update_params(27);
        mmgparams.C7 = update_params(28);
        mmgparams.C10= update_params(29);
        mmgparams.kx_rudder_reverse = update_params(30);
        mmgparams.cpr_rudder = update_params(31);

    elseif(CmaSetting.type_si==3)
        % update all valiable
        % Ndim == 47
        %choose values to be updated
        update_index =[1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, ...
                   14, ...
                   18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, ...
                   45, 46, 47, 48, 49, 50, 51, 52, 53, 54];    

    % update mmg params(hydro derivatives) by CMA-ES
     update_params =zeros(dim,1);
    for i = 1:dim
        alt = update_index(i);
        update_params(i) = hydro_derivative.params_min(alt) +(hydro_derivative.params_max(alt) - hydro_derivative.params_min(alt)) * ( Xopt(i)+ 1.0d0 ) * 0.5d0;
    end 
    
    mmgparams.massx_nd = update_params(1);
    mmgparams.massy_nd= update_params(2);
    mmgparams.IzzJzz_nd= update_params(3);
    mmgparams.Xvr_nd = update_params(4);
    mmgparams.Yv_nd = update_params(5);
    mmgparams.Yr_nd = update_params(6);
    mmgparams.Nv_nd = update_params(7);
    mmgparams.Nr_nd= update_params(8);
    mmgparams.CD = update_params(9);
    mmgparams.C_rY = update_params(10);
    mmgparams.C_rN= update_params(11);
    mmgparams.X_0A_nd= update_params(12);
    mmgparams.t_prop = update_params(13);
    mmgparams.wP0 = update_params(14);
    mmgparams.tau_prop = update_params(15);
    mmgparams.CP_nd = update_params(16);
    mmgparams.xP_nd = update_params(17);
    mmgparams.A1 = update_params(18);
    mmgparams.A2 = update_params(19);
    mmgparams.A3 = update_params(20);
    mmgparams.A4 = update_params(21);
    mmgparams.A5 = update_params(22);
    mmgparams.A6 = update_params(23);
    mmgparams.A7 = update_params(24);
    mmgparams.A8= update_params(25);
    mmgparams.B1 = update_params(26);
    mmgparams.B2 = update_params(27);
    mmgparams.B3 = update_params(28);
    mmgparams.B4 = update_params(29);
    mmgparams.B5 = update_params(30);
    mmgparams.B6 = update_params(31);
    mmgparams.B7 = update_params(32);
    mmgparams.B8= update_params(33);
    mmgparams.C3 = update_params(34);
    mmgparams.C6 = update_params(35);
    mmgparams.C7 = update_params(36);
    mmgparams.C10= update_params(37);
    mmgparams.tR = update_params(38);
    mmgparams.aH_rudder = update_params(39);
    mmgparams.xh_rudder_nd = update_params(40);
    mmgparams.lR_nd = update_params(41);
    mmgparams.kx_rudder = update_params(42);
    mmgparams.kx_rudder_reverse = update_params(43);
    mmgparams.epsilon_rudder = update_params(44);
    mmgparams.cpr_rudder = update_params(45);
    mmgparams.gammaN = update_params(46);
    mmgparams.gammaP= update_params(47);
    
    elseif(CmaSetting.type_si==4)
%  ! update forward valiable
%  ! Ndim == 22
%  ! choose values to be updated
        update_index =[5, 6, 7, 8, 9, 10, 11, 12, ...
                   14, ...
                   18, 19, 20, 21, 22, ...
                   45, 46, 47, 48, 49, 51, 53, 54];    

    % update mmg params(hydro derivatives) by CMA-ES
     update_params =zeros(dim,1);
    for i = 1:dim
        alt = update_index(i);
        update_params(i) = hydro_derivative.params_min(alt) +(hydro_derivative.params_max(alt) - hydro_derivative.params_min(alt)) * ( Xopt(i)+ 1.0d0 ) * 0.5d0;
    end 
    
    mmgparams.Xvr_nd = update_params(1);
    mmgparams.Yv_nd = update_params(2);
    mmgparams.Yr_nd = update_params(3);
    mmgparams.Nv_nd = update_params(4);
    mmgparams.Nr_nd= update_params(5);
    mmgparams.CD = update_params(6);
    mmgparams.C_rY = update_params(7);
    mmgparams.C_rN= update_params(8);
    mmgparams.X_0A_nd= update_params(9);
    mmgparams.t_prop = update_params(10);
    mmgparams.wP0 = update_params(11);
    mmgparams.tau_prop = update_params(12);
    mmgparams.CP_nd = update_params(13);
    mmgparams.xP_nd = update_params(14);
    mmgparams.tR = update_params(15);
    mmgparams.aH_rudder = update_params(16);
    mmgparams.xh_rudder_nd = update_params(17);
    mmgparams.lR_nd = update_params(18);
    mmgparams.kx_rudder = update_params(19);
    mmgparams.epsilon_rudder = update_params(20);
    mmgparams.gammaN = update_params(21);
    mmgparams.gammaP= update_params(22);
    
        elseif(CmaSetting.type_si==5)
% ! update all valiable include wind
%     ! Ndim == 54
%     ! choose values to be updated
        update_index =[5, 6, 7, 8, 9, 10, 11, 12, ...
                   14, ...
                   18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, ...
                   45, 46, 47, 48, 49, 50, 51, 52, 53, 54, ...
                   55, 56, 57, 58, 59, 60, 61, 62, 63, 64];    

    % update mmg params(hydro derivatives) by CMA-ES
     update_params =zeros(dim,1);
    for i = 1:dim
        alt = update_index(i);
        update_params(i) = hydro_derivative.params_min(alt) +(hydro_derivative.params_max(alt) - hydro_derivative.params_min(alt)) * ( Xopt(i)+ 1.0d0 ) * 0.5d0;
    end 
    
    mmgparams.Xvr_nd = update_params(1);
    mmgparams.Yv_nd = update_params(2);
    mmgparams.Yr_nd = update_params(3);
    mmgparams.Nv_nd = update_params(4);
    mmgparams.Nr_nd= update_params(5);
    mmgparams.CD = update_params(6);
    mmgparams.C_rY = update_params(7);
    mmgparams.C_rN= update_params(8);
    mmgparams.X_0A_nd= update_params(9);
    mmgparams.t_prop = update_params(10);
    mmgparams.wP0 = update_params(11);
    mmgparams.tau_prop = update_params(12);
    mmgparams.CP_nd = update_params(13);
    mmgparams.xP_nd = update_params(14);
    mmgparams.A1 = update_params(15);
    mmgparams.A2 = update_params(16);
    mmgparams.A3 = update_params(17);
    mmgparams.A4 = update_params(18);
    mmgparams.A5 = update_params(19);
    mmgparams.A6 = update_params(20);
    mmgparams.A7 = update_params(21);
    mmgparams.A8= update_params(22);
    mmgparams.B1 = update_params(23);
    mmgparams.B2 = update_params(24);
    mmgparams.B3 = update_params(25);
    mmgparams.B4 = update_params(26);
    mmgparams.B5 = update_params(27);
    mmgparams.B6 = update_params(28);
    mmgparams.B7 = update_params(29);
    mmgparams.B8= update_params(30);
    mmgparams.C3 = update_params(31);
    mmgparams.C6 = update_params(32);
    mmgparams.C7 = update_params(33);
    mmgparams.C10= update_params(34);
    mmgparams.tR = update_params(35);
    mmgparams.aH_rudder = update_params(36);
    mmgparams.xh_rudder_nd = update_params(37);
    mmgparams.lR_nd = update_params(38);
    mmgparams.kx_rudder = update_params(39);
    mmgparams.kx_rudder_reverse = update_params(40);
    mmgparams.epsilon_rudder = update_params(41);
    mmgparams.cpr_rudder = update_params(42);
    mmgparams.gammaN = update_params(43);
    mmgparams.gammaP= update_params(44);
    mmgparams.XX0 = update_params(45);
    mmgparams.XX1 = update_params(46);
    mmgparams.XX3 = update_params(47);
    mmgparams.XX5 = update_params(48);
    mmgparams.YY1 = update_params(49);
    mmgparams.YY3 = update_params(50);
    mmgparams.YY5 = update_params(51);
    mmgparams.NN1 = update_params(52);
    mmgparams.NN2 = update_params(53);
    mmgparams.NN3 = update_params(54);
    
    elseif(CmaSetting.type_si==6)
%  ! update variables related to lowspeed, backward and wind
%     ! Ndim == 41
%     ! choose values to be updated
        update_index =[5, 6, 7, 8, 9, 10, 11, 12, ...
                   14, ...
                   23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, ...
                   50, 52, ...
                   55, 56, 57, 58, 59, 60, 61, 62, 63, 64];    

    % update mmg params(hydro derivatives) by CMA-ES
     update_params =zeros(dim,1);
    for i = 1:dim
        alt = update_index(i);
        update_params(i) = hydro_derivative.params_min(alt) +(hydro_derivative.params_max(alt) - hydro_derivative.params_min(alt)) * ( Xopt(i)+ 1.0d0 ) * 0.5d0;
    end 
        mmgparams.Xvr_nd = update_params(1);
        mmgparams.Yv_nd = update_params(2);
        mmgparams.Yr_nd = update_params(3);
        mmgparams.Nv_nd = update_params(4);
        mmgparams.Nr_nd= update_params(5);
        mmgparams.CD = update_params(6);
        mmgparams.C_rY = update_params(7);
        mmgparams.C_rN= update_params(8);
        mmgparams.X_0A_nd= update_params(9);
        mmgparams.A1 = update_params(10);
        mmgparams.A2 = update_params(11);
        mmgparams.A3 = update_params(12);
        mmgparams.A4 = update_params(13);
        mmgparams.A5 = update_params(14);
        mmgparams.A6 = update_params(15);
        mmgparams.A7 = update_params(16);
        mmgparams.A8= update_params(17);
        mmgparams.B1 = update_params(18);
        mmgparams.B2 = update_params(19);
        mmgparams.B3 = update_params(20);
        mmgparams.B4 = update_params(21);
        mmgparams.B5 = update_params(22);
        mmgparams.B6 = update_params(23);
        mmgparams.B7 = update_params(24);
        mmgparams.B8= update_params(25);
        mmgparams.C3 = update_params(26);
        mmgparams.C6 = update_params(27);
        mmgparams.C7 = update_params(28);
        mmgparams.C10= update_params(29);
        mmgparams.kx_rudder_reverse = update_params(30);
        mmgparams.cpr_rudder = update_params(31);
        mmgparams.XX0 = update_params(32);
        mmgparams.XX1 = update_params(33);
        mmgparams.XX3 = update_params(34);
        mmgparams.XX5 = update_params(35);
        mmgparams.YY1 = update_params(36);
        mmgparams.YY3 = update_params(37);
        mmgparams.YY5 = update_params(38);
        mmgparams.NN1 = update_params(39);
        mmgparams.NN2 = update_params(40);
        mmgparams.NN3 = update_params(41);
        
    elseif(CmaSetting.type_si==7)
%         ! update all valiable including wind
%         ! Ndim == 57
%     ! choose values to be updated
        update_index =[1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, ...
                   14, ...
                   18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, ...
                   45, 46, 47, 48, 49, 50, 51, 52, 53, 54, ...
                   55, 56, 57, 58, 59, 60, 61, 62, 63, 64];    

    % update mmg params(hydro derivatives) by CMA-ES
     update_params =zeros(dim,1);
    for i = 1:dim
        alt = update_index(i);
        update_params(i) = hydro_derivative.params_min(alt) +(hydro_derivative.params_max(alt) - hydro_derivative.params_min(alt)) * ( Xopt(i)+ 1.0d0 ) * 0.5d0;
    end 
     
    mmgparams.massx_nd = update_params(1);
    mmgparams.massy_nd= update_params(2);
    mmgparams.IzzJzz_nd= update_params(3);
    mmgparams.Xvr_nd = update_params(4);
    mmgparams.Yv_nd = update_params(5);
    mmgparams.Yr_nd = update_params(6);
    mmgparams.Nv_nd = update_params(7);
    mmgparams.Nr_nd= update_params(8);
    mmgparams.CD = update_params(9);
    mmgparams.C_rY = update_params(10);
    mmgparams.C_rN= update_params(11);
    mmgparams.X_0A_nd= update_params(12);
    mmgparams.t_prop = update_params(13);
    mmgparams.wP0 = update_params(14);
    mmgparams.tau_prop = update_params(15);
    mmgparams.CP_nd = update_params(16);
    mmgparams.xP_nd = update_params(17);
    mmgparams.A1 = update_params(18);
    mmgparams.A2 = update_params(19);
    mmgparams.A3 = update_params(20);
    mmgparams.A4 = update_params(21);
    mmgparams.A5 = update_params(22);
    mmgparams.A6 = update_params(23);
    mmgparams.A7 = update_params(24);
    mmgparams.A8= update_params(25);
    mmgparams.B1 = update_params(26);
    mmgparams.B2 = update_params(27);
    mmgparams.B3 = update_params(28);
    mmgparams.B4 = update_params(29);
    mmgparams.B5 = update_params(30);
    mmgparams.B6 = update_params(31);
    mmgparams.B7 = update_params(32);
    mmgparams.B8= update_params(33);
    mmgparams.C3 = update_params(34);
    mmgparams.C6 = update_params(35);
    mmgparams.C7 = update_params(36);
    mmgparams.C10= update_params(37);
    mmgparams.tR = update_params(38);
    mmgparams.aH_rudder = update_params(39);
    mmgparams.xh_rudder_nd = update_params(40);
    mmgparams.lR_nd = update_params(41);
    mmgparams.kx_rudder = update_params(42);
    mmgparams.kx_rudder_reverse = update_params(43);
    mmgparams.epsilon_rudder = update_params(44);
    mmgparams.cpr_rudder = update_params(45);
    mmgparams.gammaN = update_params(46);
    mmgparams.gammaP= update_params(47);
    mmgparams.XX0 = update_params(48);
    mmgparams.XX1 = update_params(49);
    mmgparams.XX3 = update_params(50);
    mmgparams.XX5 = update_params(51);
    mmgparams.YY1 = update_params(52);
    mmgparams.YY3 = update_params(53);
    mmgparams.YY5 = update_params(54);
    mmgparams.NN1 = update_params(55);
    mmgparams.NN2 = update_params(56);
    mmgparams.NN3 = update_params(57);
    
    elseif(CmaSetting.type_si==8)
%     ! update forward valiable include wind
%     ! Ndim == 32
%     ! choose values to be updated
        update_index =[5, 6, 7, 8, 9, 10, 11, 12, ...
                   14, ...
                   18, 19, 20, 21, 22, ...
                   45, 46, 47, 48, 49, 51, 53, 54, ...
                    55, 56, 57, 58, 59, 60, 61, 62, 63, 64];    

    % update mmg params(hydro derivatives) by CMA-ES
     update_params =zeros(dim,1);
    for i = 1:dim
        alt = update_index(i);
        update_params(i) = hydro_derivative.params_min(alt) +(hydro_derivative.params_max(alt) - hydro_derivative.params_min(alt)) * ( Xopt(i)+ 1.0d0 ) * 0.5d0;
    end 
    
    mmgparams.Xvr_nd = update_params(1);
    mmgparams.Yv_nd = update_params(2);
    mmgparams.Yr_nd = update_params(3);
    mmgparams.Nv_nd = update_params(4);
    mmgparams.Nr_nd= update_params(5);
    mmgparams.CD = update_params(6);
    mmgparams.C_rY = update_params(7);
    mmgparams.C_rN= update_params(8);
    mmgparams.X_0A_nd= update_params(9);
    mmgparams.t_prop = update_params(10);
    mmgparams.wP0 = update_params(11);
    mmgparams.tau_prop = update_params(12);
    mmgparams.CP_nd = update_params(13);
    mmgparams.xP_nd = update_params(14);
    mmgparams.tR = update_params(15);
    mmgparams.aH_rudder = update_params(16);
    mmgparams.xh_rudder_nd = update_params(17);
    mmgparams.lR_nd = update_params(18);
    mmgparams.kx_rudder = update_params(19);
    mmgparams.epsilon_rudder = update_params(20);
    mmgparams.gammaN = update_params(21);
    mmgparams.gammaP= update_params(22);
    mmgparams.XX0 = update_params(23);
    mmgparams.XX1 = update_params(24);
    mmgparams.XX3 = update_params(25);
    mmgparams.XX5 = update_params(26);
    mmgparams.YY1 = update_params(27);
    mmgparams.YY3 = update_params(28);
    mmgparams.YY5 = update_params(29);
    mmgparams.NN1 = update_params(30);
    mmgparams.NN2 = update_params(31);
    mmgparams.NN3 = update_params(32);
    
    end
end