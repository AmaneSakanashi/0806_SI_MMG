function [F_right] = Right_MMGLS(time, x, delta_rudder, n_prop, pp, mmgparams,...
                                 switch_wind, Wind_Direction, AVG_Wind_Velocity)

%  % read State Variables and control
Xpos  = x(1);
u_velo     = x(2);
Ypos  = x(3);
vm_velo    = x(4);
psi   = x(5);
r_angvelo     = x(6);
U_ship = sqrt(u_velo ^ 2.0d0 + vm_velo ^ 2.0d0); 
beta  = atan2(vm_velo, u_velo); 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%  Principal Particulars and parameters
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% principal particular

    %%% tentative! save as structure
    lpp  = pp.lpp;
    breadth = pp.lpp;
    
    
rho_fresh  = 1000.0/9.80665;
lpp = 3.0d0;
breadth = 0.48925d0;
draft = 0.20114d0;
Mass_nd  = 0.27078d0;
xG_nd    = 0.03119d0;
%% prop
Dp    =  0.084d0;
pitchr =  0.7151d0;

% Rudder
AREA    =  0.01063d0;
lambda_rudder  =  1.539d0;
xRR     = -1.5d0;
% wind force
   
    AT    = 0.13544d0;  %    AT --- Transverse Projected area (m^2)
    AL    = 0.52038d0;  %    AL --- Lateral Projected area (m^2)
    AOD   = 0.14230d0;  %    AOD -- Lateral Projected area of superstructure Ass and LNG tanks,
                        %           container etc. on the deck (m^2)
                        %           Here, Ass is defined as lateral Projected area of superstructure (m^2)
    LCW   = -0.01299d0; %    C ---- Distance from midship section to center of lateral projected area (m)
    LCBR  = -0.06209d0; %    CBR -- Distance from midship section to center of the Ass (m)
    HBR   = 0.294d0;    %    HBR -- Height from free surface to top of the superstructure (bridge) (m) 
    HC    = 0.10359d0;  %    HC --- Height to center of lateral projected area (m)
    SBW   = breadth;         %    SBW -- breadth for wind area (in usual, it coincides with breadth, but it is defined for Trimaran vessel)
    Lz    = 0.0d0;      %    Lz --- Acting position of sway force from center of gravity (it is necessaty to calculate roll moment)


%%
  
% variables to add dimenstion
D1_ad2 = 0.5 * rho_fresh * lpp ^ 2 * draft;      %
D2_ad2 = 0.5 * rho_fresh * lpp ^ 4 * draft;      %
Mass  = Mass_nd * D1_ad2;
MassX =mmgparams.massx_nd * D1_ad2;
MassY = mmgparams.massy_nd * D1_ad2;
IJzz  = mmgparams.IzzJzz_nd * D2_ad2;
xG    = xG_nd * lpp;

%	
if (abs(U_ship) < 1.0d-5)
    v_nd = 0.0000001;
    r_nd = 0.0000001;
    U_ship    = 0.0000001;
else
    v_nd = vm_velo / U_ship; 
    r_nd = r_angvelo * lpp / U_ship;
end 

    %%%%%%%%%%%%%%%%%%%%%
    %%% Force of Hull
    %%%%%%%%%%%%%%%%%%%%%


%integration parameters
i_max = 200;
Y_ad = 0.0;
N_ad = 0.0;

for i = 1:i_max
    x0 = -0.5+(i-1)/i_max;
    x1 = -0.5+(i)  /i_max;
    %   In order to apply this model to super low speed, erase u_nd, v_nd and r_nd.
    comp0 = vm_velo + mmgparams.C_rY * r_angvelo * lpp * x0; % v -->> vm_velo (definition change 2018/4/5)
    comp1 = vm_velo + mmgparams.C_rY * r_angvelo * lpp * x1; % v -->> vm_velo (definition change 2018/4/5)
    comp2 = vm_velo + mmgparams.C_rN * r_angvelo * lpp * x0; % v -->> vm_velo (definition change 2018/4/5)
    comp3 = vm_velo + mmgparams.C_rN * r_angvelo * lpp * x1; % v -->> vm_velo (definition change 2018/4/5)
% %   EFD description of Hull force
%     comp0=v_nd+mmgparams.C_rY*r_nd*x0;
%     comp1=v_nd+mmgparams.C_rY*r_nd*x1;
%     comp2=v_nd+mmgparams.C_rN*r_nd*x0;
%     comp3=v_nd+mmgparams.C_rN*r_nd*x1;
    %
    Y_ad = Y_ad + 0.5 * ( abs( comp0 ) * comp0 +abs( comp1 ) * comp1 ) / i_max;
    N_ad = N_ad + 0.5 * ( abs( comp2 ) * comp2 * x0 + abs( comp3 ) * comp3 * x1 )/i_max;
end

YHN_nd = -mmgparams.CD*Y_ad;
NHN_nd = -mmgparams.CD*N_ad;

XH = 0.5d0 * rho_fresh * lpp * draft * ( ( mmgparams.X_0F_nd + ( mmgparams.X_0A_nd - mmgparams.X_0F_nd ) * ( abs( beta ) / pi ) ) * u_velo * U_ship + mmgparams.Xvr_nd * vm_velo*r_angvelo*lpp);
YH = 0.5d0 * rho_fresh * lpp * draft * ( mmgparams.Yv_nd * vm_velo * abs( u_velo ) + mmgparams.Yr_nd * r_angvelo * lpp * u_velo + YHN_nd );
NH = 0.5d0 * rho_fresh * lpp ^ 2 *draft * ( mmgparams.Nv_nd * vm_velo * u_velo + mmgparams.Nr_nd * r_angvelo * lpp * abs( u_velo ) + NHN_nd );


%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Force of Propeller %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%

Jsyn  = -0.35d0;
Jsyn0 = -0.06d0;
pitch = Dp * pitchr;
%
wP     = mmgparams.wP0 - mmgparams.tau_prop * abs( v_nd + mmgparams.xP_nd * r_nd ) - mmgparams.CP_nd * ( v_nd + mmgparams.xP_nd * r_nd ) ^ 2;
if (u_velo>0.0d0) 
        one_minus_wprop = 1.0d0 - wP;
else
        one_minus_wprop = 1.0d0; %u = up at 2nd and 4th quad.
end

if( one_minus_wprop > 1.0d0 ) 
        one_minus_wprop = 1.0d0;
elseif( one_minus_wprop <= 0.0d0 ) 
        one_minus_wprop = eps(1.0d0);
end

uprop = u_velo * one_minus_wprop;           
if(abs(n_prop) < eps(1.0d0))
       Js = 1.0d+4;
elseif(abs(u_velo) < eps(1.0d0))
       Js = eps(1.0d0) ;
else
       Js = u_velo / (Dp * n_prop);
end
%
 if  n_prop>=0.0d0  &&  u_velo>=0.0d0   % 1st quad
       J_prop  = Js * one_minus_wprop;
       KT = mmgparams.kt_coeff0 + mmgparams.kt_coeff1 * J_prop + mmgparams.kt_coeff2 * J_prop ^ 2.0d0;
       XP = rho_fresh * Dp ^ 4.0d0 * n_prop ^ 2.0d0 * (1.0d0 - mmgparams.t_prop) * KT;
       YP = 0.0d0;
       NP = 0.0d0;
 elseif n_prop>=0.0d0 && u_velo<0.0d0  % 2nd quad
       J_prop = Js;
       KT = mmgparams.kt_coeff0 + mmgparams.kt_coeff1 * J_prop + mmgparams.kt_coeff2 * J_prop ^ 2.0d0;
       XP = rho_fresh * Dp ^ 4.0d0 * n_prop ^ 2.0d0 * (1.0d0 - mmgparams.t_prop) * KT;
       YP = 0.5d0 * rho_fresh * lpp    * draft * (n_prop * pitch)^2.0d0 * (mmgparams.A6 * Js^2.0d0 + mmgparams.A7 * Js + mmgparams.A8);
       NP = 0.5d0 * rho_fresh * lpp^2.0d0 * draft * (n_prop * pitch)^2.0d0 * (mmgparams.B6 * Js^2.0d0 + mmgparams.B7 * Js + mmgparams.B8);
 else % 3rd & 4th quad(n <0)
         if(u_velo>=0.0d0)
            J_prop  = Js * one_minus_wprop;
        else
            J_prop  = Js ;
         end
        if Js >= mmgparams.C10
          XP = rho_fresh * n_prop ^ 2.0d0 * Dp ^ 4.0d0 * (mmgparams.C6 + mmgparams.C7 * Js);
       else
          XP = rho_fresh * n_prop ^ 2.0d0 * Dp ^ 4.0d0 * mmgparams.C3;
       end
       %
       
       if(Jsyn<=Js && Js<=Jsyn0)
          YP = 0.5d0 * rho_fresh * lpp          * draft * (n_prop * Dp) ^ 2.0d0 * (mmgparams.A1 + mmgparams.A2 * Js);
          NP = 0.5d0 * rho_fresh * lpp ^ 2.0d0 * draft * (n_prop * Dp) ^ 2.0d0 * (mmgparams.B1 + mmgparams.B2 * Js);  
       elseif(Js<Jsyn)
          YP = 0.5d0 * rho_fresh * lpp          * draft * (n_prop * Dp) ^ 2.0d0 * (mmgparams.A3 + mmgparams.A4 * Js);
          NP = 0.5d0 * rho_fresh * lpp ^ 2.0d0 * draft * (n_prop * Dp) ^ 2.0d0 * (mmgparams.B3 + mmgparams.B4 * Js);
       elseif(Jsyn0<Js)
          YP = 0.5d0 * rho_fresh * lpp          * draft * (n_prop * Dp) ^ 2.0d0 * mmgparams.A5;
          NP = 0.5d0 * rho_fresh * lpp ^ 2.0d0 * draft * (n_prop * Dp) ^ 2.0d0 * mmgparams.B5;
       else
           print('stio')
       end
 end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Force of Rudder    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %  parameters 
    hight_rudder = sqrt(AREA * lambda_rudder);
    eta_rudder = Dp/hight_rudder;

    xh_rudder = mmgparams.xh_rudder_nd * lpp;
    lR  = mmgparams.lR_nd * lpp;
    %    Fujii's formula
    fa = 6.13d0 * lambda_rudder / (2.25d0 + lambda_rudder);
    
    % Calculate effective velocity to rudder uR
    kappa_rudder = mmgparams.kx_rudder / mmgparams.epsilon_rudder ;     
    one_minus_wrudder = mmgparams.epsilon_rudder * one_minus_wprop;
    
    % compute KT from XP to consider about 2nd,3rd,4th Quad
    if(n_prop>eps(1.0d0))
        KT = XP / (rho_fresh * Dp^4.0d0 * n_prop^2.0d0 * (1-mmgparams.t_prop));
    elseif(abs(n_prop)>eps(1.0d0))
        KT = XP / (rho_fresh * Dp^4.0d0 * n_prop^2.0d0 );
    else
        KT =0.0d0;
    end
    
    if( n_prop >= 0.0d0 && KT>0.0d0)
        uR = mmgparams.epsilon_rudder * sqrt(eta_rudder * (uprop+ kappa_rudder *  ...
                (sqrt(uprop^2 + 8.0d0 * KT * n_prop^2 * Dp ^2/ (pi )) - uprop))^2 + (1- eta_rudder) * uprop^2); %!!!<= normarl mmg model for low speed (Yoshimura's)
    else % n<0 
    % Kitagawa's model for uR in n<0
        if(u_velo < 0.0d0)  %4th quad
            uR = u_velo;
        else                   %3rd quad
            urpr1 = u_velo * one_minus_wrudder + n_prop * Dp * mmgparams.kx_rudder_reverse * sqrt(8.0d0 * abs(KT)/ pi);
            urpr2 = u_velo * one_minus_wrudder;
            ursq  = eta_rudder * sign(urpr1) * urpr1^2 + (1- eta_rudder) * urpr2^2.0d0 + mmgparams.cpr_rudder * u_velo^2.0d0;
            uR =  sign(ursq)*sqrt(abs(ursq)) ;   
        end
    end
    if(vm_velo+xRR*r_angvelo>=0.0d0) 
        vR = -1.0d0 * mmgparams.gammaP * (vm_velo + lR * r_angvelo);
    else
        vR = -1.0d0 *mmgparams.gammaN * (vm_velo + lR * r_angvelo);
    end
    UUR = sqrt(uR ^ 2.0d0 + vR ^ 2.0d0);
    aR = delta_rudder - atan2(vR, uR);
    FN = 0.5d0 * rho_fresh * AREA * fa * UUR ^ 2.0d0 * sin(aR);
    %    Rudder forces and moments
    XR = - (1.0d0 - mmgparams.tR) * FN * sin(delta_rudder);
    YR = - (1.0d0 + mmgparams.aH_rudder) * FN * cos(delta_rudder);
    NR = - (xRR + mmgparams.aH_rudder * xh_rudder) * FN * cos(delta_rudder);

%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% Force of Bow and stern thruster    %%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Bow and stern thruster (Super simple modelling)
%     % Here, KT characteristics is assumed as symmetry in forward and back side.
% 
% Dthrust =  0.084d0 * 0.25d0; % Assumption (qurter size%)
% tthrust =  0.0d0;
% wthrust =  0.0d0; % It was changed at 01/31/2018
% lB      =  lpp * 0.45d0;
% lS      = - lpp * 0.45d0;
% %   Bow thruster
% if Brps == 0.0d0
%    AJB = 0.0d0;
% else
%    AJB = ( vm_velo + lB * r_angvelo ) / ( Dthrust * abs( Brps ) ) * ( 1.0d0 - wthrust );
% end
% %
% KTB = mmgparams.kt_coeff0 + mmgparams.kt_coeff1 * AJB + mmgparams.kt_coeff2 * AJB ^ 2.0d0;
% %
% YBT = rho_fresh * Dthrust ^ 4.0d0 * Brps ^ 2.0d0 * ( 1.0d0 - tthrust ) * KTB * sign( Brps );
% NBT = YBT * lB;
% %   Stern thruster
% if Srps == 0.0d0
%    AJS = 0.0d0;
% else
%    AJS = ( vm_velo + lS * r_angvelo ) / ( Dthrust * abs( Srps ) ) * ( 1.0d0 - wthrust );
% end
% %
% KTS= mmgparams.kt_coeff0 + mmgparams.kt_coeff1 * AJS + mmgparams.kt_coeff2 * AJS ^ 2.0d0;
% %
% YST= rho_fresh * Dthrust ^ 4.0d0 * Srps ^ 2.0d0 * ( 1.0d0 - tthrust ) * KTS * sign( Srps );
% NST= YST * lS;
% %   Restriction of force and moment due to forward speed
% if abs( u_velo ) > Thruster_speed_max
%     YBT = 0.0d0;
%     NBT = 0.0d0;
%     YST = 0.0d0;
%     NST = 0.0d0;
% end
% %
% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% if IBowThruster == 0
%     YBT = 0.0d0;
%     NBT = 0.0d0;
% end
% if ISternThruster == 0
%     YST = 0.0d0;
%     NST = 0.0d0;
% end
% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Force of wind      %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Save as Windcoef
    Windcoef(1) = AT;   Windcoef(2) = AL;  Windcoef(3) = AOD; Windcoef(4) = LCW;
    Windcoef(5) = LCBR; Windcoef(6) = HBR; Windcoef(7) = HC;  Windcoef(8) = SBW;%
    Wind_Velocity = AVG_Wind_Velocity;
    % if switch_wind == 2
    %     for i = 1: DivOme
    %         %
    %         Wind_Velocity = Wind_Velocity + ampwind (i) * sin(Omegaswind(i) * time + Randomwind(i));
    %     end
    % end

% compute apparent wind
[ Relative_Wind_Velocity, Rerative_Wind_Direction] = true2apparent(Wind_Velocity, U_ship, Wind_Direction, beta, psi);
Angle_of_attack = 2.0d0 * pi - Rerative_Wind_Direction;

%   Calculation of wind force and moment by Fujiwara's method
% [CXwind, CYwind, CNwind, CKwind, FXwind, FYwind, FNwind, FKwind] = ...
% WindForce(Relative_Wind_Velocity, Angle_of_attack, lpp, Lz, Windcoef);

%  Calculation of wind force and moment by Fujiwara's method: coefficient
%  optimized by CMA-ES
    KK1 = 0.0d0;
    KK2 = 0.0d0;
    KK3 = 0.0d0;
    KK5 = 0.0d0;
    
[~, ~, ~, ~, FXwind, FYwind, FNwind, ~] = ...
WindForceOnly(Relative_Wind_Velocity, Angle_of_attack, lpp, Lz, Windcoef, mmgparams.XX0, mmgparams.XX1, mmgparams.XX3, mmgparams.XX5, ...
mmgparams.YY1, mmgparams.YY3, mmgparams.YY5, mmgparams.NN1, mmgparams.NN2, mmgparams.NN3,...
KK1, KK2, KK3, KK5);

XA = FXwind; % EFD positive definition of FXwind is backward. Therefore, it should be converted. 
YA = FYwind;
NA = FNwind;

%
if switch_wind == 0
    XA = 0.0d0;
    YA = 0.0d0;
    NA = 0.0d0;
end
%
%   Summation of every force and moment
X = XH + XP + XA + XR;
Y = YH + YP + YA + YR;
N = NH + NP + NA + NR;

AA1 = Mass+MassY;
AA2 = xG*Mass;
AA3 = Y-(Mass+MassX)*u_velo*r_angvelo;
BB1 = IJzz+xG^2*Mass;
BB2 = xG*Mass;
BB3 = N-xG*Mass*u_velo*r_angvelo;
u_dot  = (X+(Mass+MassY)*vm_velo*r_angvelo+xG*Mass*r_angvelo^2)/(Mass+MassX);
vm_dot = (AA3*BB1-AA2*BB3)/(AA1*BB1-AA2*BB2);
r_dot  = (AA3*AA2-BB3*AA1)/(AA2*BB2-AA1*BB1);
%
phi_dot   = r_angvelo;
X_dot     = u_velo*cos(psi)-vm_velo*sin(psi); % v -->> vm_velo (2018/4/5)
Y_dot     = u_velo*sin(psi)+vm_velo*cos(psi); % v -->> vm_velo (2018/4/5)

F_right(1, 1) = X_dot;
F_right(2, 1) = u_dot;
F_right(3, 1) = Y_dot;
F_right(4, 1) = vm_dot;
F_right(5, 1) = phi_dot;
F_right(6, 1) = r_dot;
 
end
%%
function[UA, gamma_a] = true2apparent(UT, U_ship, zeta_wind, beta, psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute apparent(relative) wind (or current) speed and direction from true wind measured on ship
% !!! angle are RADIANS !!!
%
% [input]
% UT      : norm of true wind speed (speed to ground).
% U_ship  : norm of speed of arbitary point on ship [m/s], given by GNSS
% zeta_wind:true wind direction (direction to ground). 0 at x direction of
%           earth fix coordinate. (depend on earth fix coor. system. NOT north = zero)
% beta    : angle between ship heading and U_ship (drift angle) [rad] [0, 2pi].
% psi     : heading of ship related to earth fix coordinate [rad] [-pi, pi].
%
% [output]
% UA      : apparent wind speed (norm) on midship. [m/s]
% gamma_a : apparent wind direction on midship.[rad] [0, 2pi]. 0 at ship
% heading direction. clock wise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute chi(agle between ship heading and true wind direction)
    chi = mod(3*pi - zeta_wind + psi, 2*pi);
    
    if beta < 0
        beta = beta + 2*pi;
    end
    
    % compute chi_beta(angle between ship traveling direction and true wind)
    chi_beta = chi+beta ;
    
    % convert chi_beta tp [0, 2pi]
    if chi_beta < 0
        chi_beta = chi_beta + 2*pi;
    elseif chi_beta > 2*pi
        chi_beta = chi_beta - 2*pi;
    end
    
    % compute apparent wind speed
    UA = sqrt(UT^2 + U_ship^2 -2*UT*U_ship*cos(chi_beta));
    
    % avoid numerical error by acos (acos input are [-1 1])
    RV = (UA^2 + U_ship^2 -UT^2)/(2*UA*U_ship);
    if abs(RV) >1
       RV = sign(RV);
    end
    
    % compute gamma_beta(angle between ship traveling direction and apparent wind)
    % use if sentence because acos range is (0 pi) and gamma_beta (0 2pi)
    if chi_beta <= pi
        gamma_beta = acos(RV);
    else
        gamma_beta = 2*pi-acos(RV);
    end
    
    % compute apparent wind direction gamma_a 
    gamma_a = gamma_beta + beta; 
    
    % avoid NaN at wind or shipspeed are zero
    if UT < 1e-5
        gamma_a = beta;
    elseif U_ship < 1e-5
        gamma_a = mod(3*pi-chi, 2*pi);
    end
    
    % convert to [0, 2pi]
    if gamma_a > 2*pi
        gamma_a = gamma_a-2*pi;
    elseif gamma_a < 0
        gamma_a = gamma_a + 2*pi;
    end   
end


function [CXwind, CYwind, CNwind, CKwind,FXwind, FYwind, FNwind, FKwind] = ...
WindForce(V, psi, lpp, Lz, Windcoef)
%
%	Difinition of physical values
grav = 9.80665d0; 
rhoA = 1.220d0 / grav;

%
%   Wind coefficients
AT   = Windcoef(1); AL  = Windcoef(2); AOD = Windcoef(3); LC  = Windcoef(4);
LCBR = Windcoef(5); HBR = Windcoef(6); HC  = Windcoef(7); SBW = Windcoef(8);
%   Coefficients of Fujiwara's regression
%   X-directional Coefficients
X00 = -0.330D0;  X01 =  0.293D0;  X02 =  0.0193D0; X03 =  0.6820D0;
X10 = -1.353D0;  X11 =  1.700D0;  X12 =  2.8700D0; X13 = -0.4630D0;
X14 = -0.570D0;  X15 = -6.640D0;  X16 = -0.0123D0; X17 =  0.0202D0;
X30 =  0.830D0;  X31 = -0.413D0;  X32 = -0.0827D0; X33 = -0.5630D0;
X34 =  0.804D0;  X35 = -5.670D0;  X36 =  0.0401D0; X37 = -0.1320D0;
X50 =  0.0372D0; X51 = -0.0075D0; X52 = -0.1030D0; X53 =  0.0921D0;
%   Y-directional Coefficients
Y10 =  0.684d0;  Y11 =   0.717d0;  Y12 = -3.2200d0; Y13 =  0.0281d0; Y14 =  0.0661d0; Y15 =  0.2980d0;
Y30 = -0.400d0;  Y31 =   0.282d0;  Y32 =  0.3070d0; Y33 =  0.0519d0; Y34 =  0.0526d0; Y35 = -0.0814d0; Y36 =  0.0582d0;
Y50 =  0.122d0;  Y51 =  -0.166d0;  Y52 = -0.0054d0; Y53 = -0.0481d0; Y54 = -0.0136d0; Y55 =  0.0864d0; Y56 = -0.0297d0;
%   N-directional Coefficients
N10 =  0.2990d0; N11 =   1.710d0;  N12 =  0.183d0;  N13 = -1.09d0;   N14 = -0.0442d0; N15 = -0.289d0;  N16 =  4.24d0;  N17 = -0.0646d0; N18 =  0.0306d0;
N20 =  0.1170d0; N21 =   0.123d0;  N22 = -0.323d0;  N23 =  0.0041d0; N24 = -0.166d0;  N25 = -0.0109d0; N26 =  0.174d0; N27 =  0.214d0;  N28 = -1.06d0; 
N30 =  0.0230d0; N31 =   0.0385d0; N32 = -0.0339d0; N33 =  0.0023d0; 
%   K-directional Coefficients
K10 =  3.63d0;   K11 = -30.7d0;    K12 = 16.8d0;    K13 =  3.270d0;  K14 = -3.03d0;   K15 =  0.552d0;  K16 = -3.03d0;   K17 = 1.82d0;   K18 = -0.224d0; 
K20 = -0.480d0;  K21 =   0.166d0;  K22 =  0.318d0;  K23 =  0.132d0;  K24 = -0.148d0;  K25 =  0.408d0;  K26 = -0.0394d0; K27 = 0.0041d0; 
K30 =  0.164d0;  K31 =  -0.170d0;  K32 =  0.0803d0; K33 =  4.920d0;  K34 = -1.780d0;  K35 =  0.0404d0; K36 = -0.739d0; 
K50 =  0.449d0;  K51 =  -0.148d0;  K52 = -0.0049d0; K53 = -0.396d0;  K54 = -0.0109d0; K55 = -0.0726d0;
%  
XX0 = X00 + X01 * (SBW * HBR / AT) + X02 * (LC / HC) + X03 * (AOD / lpp / lpp);
XX1 = X10 + X11 * (AL / lpp / SBW) + X12 * (lpp * HC / AL) + X13 * (lpp * HBR / AL) + X14 * (AOD / AL) + X15 * (AT / lpp / SBW) +...
      X16 * (lpp * lpp / AT) + X17 * (lpp / HC);
XX3 = X30 + X31 * (AL / lpp / HBR) + X32 * (AL / AT) + X33 * (lpp * HC / AL) + X34 * (AOD / AL) + X35 * (AOD / lpp / lpp) +...
      X36 * (LC / HC) + X37 * (LCBR / lpp);
XX5 = X50 + X51 * (AL / AOD) + X52 * (LCBR / lpp) + X53 * (AL / lpp / SBW);
%
YY1 = Y10 + Y11 * (LCBR / lpp) + Y12 * (LC / lpp) + Y13 * (AL / AOD) + Y14 * (LC / HC) + Y15 * (AT / (SBW * HBR));
YY3 = Y30 + Y31 * (AL / (lpp * SBW)) + Y32 * (lpp * HC / AL) + Y33 * (LCBR / lpp) + Y34 * (SBW / HBR) + Y35 * (AOD / AL) + Y36 * (AT / (SBW * HBR));
YY5 = Y50 + Y51 * (AL / (lpp * SBW)) + Y52 * (lpp / HBR) + Y53 * (LCBR / lpp) + Y54 * (SBW ^ 2 / AT) + Y55 * (LC / lpp) + Y56 * (LC * HC / AL);
%
NN1 = N10 + N11 * (LC / lpp) + N12 * (lpp * HC / AL) + N13 * (AT / AL) + N14 * (LC / HC) + N15 * (AL / (lpp * SBW)) + N16 * (AT / lpp ^ 2) + N17 * (SBW ^ 2 / AT) + N18 * (LCBR / lpp);
NN2 = N20 + N21 * (LCBR / lpp) + N22 * (LC / lpp) + N23 * (AL / AOD) + N24 * (AT / SBW ^ 2) + N25 * (lpp / HBR) + N26 * (AT / (SBW * HBR)) + N27 * (AL / (lpp * SBW)) + N28 * (AL / lpp ^ 2);
NN3 = N30 + N31  *(LCBR / lpp) + N32 * (AT / (SBW * HBR)) + N33 * (AL / AT);
%
KK1 =K10 + K11 * (HBR / lpp) + K12 * (AT / (lpp * SBW)) + K13 * (lpp * HC / AL) + K14 * (LC / lpp) + K15 * (LCBR / lpp) + K16 * (SBW / HBR) + K17 * (SBW ^ 2 / AT) + K18 * (lpp / SBW);
KK2 =K20 + K21 * (SBW / HBR) + K22 * (AT / SBW ^ 2) + K23 * (AL / (lpp * HC)) + K24 * (LCBR / lpp) + K25 * (HBR * LC / AL) + K26 * (lpp / SBW) + K27 * (lpp ^ 2 / AL);
KK3 =K30 + K31 * (SBW ^ 2 / AT) + K32 * (LCBR / lpp) + K33 * (HC / lpp) + K34 * (AT / (lpp * SBW)) + K35 * (lpp * SBW / AL) + K36 * (AOD / lpp ^ 2);
KK5 =K50 + K51 * (AL / (lpp * HC)) + K52 * (AL / AOD) + K53 * (AT / AL) + K54 * (lpp / SBW) + K55 * (AL / (lpp * SBW));
%
%  Cal of non-dimentionalized coefficients
CXwind = XX0 + XX1 * cos(psi) + XX3 * cos(3.0D0 * psi) + XX5 * cos(5.0D0 * psi);
CYwind = YY1 * sin(psi) + YY3 * sin(3.0D0 * psi) + YY5 * sin(5.0D0 * psi); % YY5 is corrected (XX5 -->> YY5) 05/13/2019
CNwind = NN1 * sin(psi) + NN2 * sin(2.0D0 * psi) + NN3 * sin(3.0D0 * psi);
CKwind = KK1 * sin(psi) + KK2 * sin(2.0D0 * psi) + KK3 * sin(3.0D0 * psi) + KK5 * sin(5.0D0 * psi);
%  Dimentionalization
FXwind = CXwind * (0.5D0 * rhoA * V * V) * AT;
FYwind = CYwind * (0.5D0 * rhoA * V * V) * AL;
FNwind = CNwind * (0.5D0 * rhoA * V * V) * lpp *AL;
FKwind = CKwind * (0.5D0 * rhoA * V * V) * AL * (AL / lpp);
%  Convert K morment around G
FKwind = FKwind + FYwind * Lz;
CKwind = FKwind / ((0.5D0 * rhoA * V * V) * AL * (AL / lpp));
%
end

function [CXwind, CYwind, CNwind, CKwind,FXwind, FYwind, FNwind, FKwind] = ...
WindForceOnly(V, psi, lpp, Lz, Windcoef, XX0, XX1, XX3, XX5, YY1, YY3, YY5, NN1, NN2, NN3, ...
KK1, KK2, KK3, KK5)
%
%	Difinition of physical values
grav = 9.80665d0; 
rhoA = 1.220d0 / grav;

%   Wind coefficients
AT   = Windcoef(1); AL  = Windcoef(2); AOD = Windcoef(3); LC  = Windcoef(4);
LCBR = Windcoef(5); HBR = Windcoef(6); HC  = Windcoef(7); SBW = Windcoef(8);
%
%  Cal of non-dimentionalized coefficients
CXwind = XX0 + XX1 * cos(psi) + XX3 * cos(3.0D0 * psi) + XX5 * cos(5.0D0 * psi);
CYwind = YY1 * sin(psi) + YY3 * sin(3.0D0 * psi) + YY5 * sin(5.0D0 * psi); % YY5 is corrected (XX5 -->> YY5) 05/13/2019
CNwind = NN1 * sin(psi) + NN2 * sin(2.0D0 * psi) + NN3 * sin(3.0D0 * psi);
CKwind = KK1 * sin(psi) + KK2 * sin(2.0D0 * psi) + KK3 * sin(3.0D0 * psi) + KK5 * sin(5.0D0 * psi);
%  Dimentionalization
FXwind = CXwind * (0.5D0 * rhoA * V * V) * AT;
FYwind = CYwind * (0.5D0 * rhoA * V * V) * AL;
FNwind = CNwind * (0.5D0 * rhoA * V * V) * lpp *AL;
FKwind = CKwind * (0.5D0 * rhoA * V * V) * AL * (AL / lpp);
%  Convert K morment around G
FKwind = FKwind + FYwind * Lz;
CKwind = FKwind / ((0.5D0 * rhoA * V * V) * AL * (AL / lpp));
%
end