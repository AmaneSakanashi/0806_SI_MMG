close all
clear 
clc

%%% set default for plot
set(0,'defaultAxesFontSize',20);
set(0,'defaultAxesFontName','times new roman');
set(0,'defaultAxesXcolor','k');
set(0,'defaultAxesYcolor','k');
set(0,'defaultAxesZcolor','k');
set(0,'defaultTextFontSize',20);
set(0,'defaultTextFontName','times new roman');
set(0,'defaultTextColor','k');
set(0,'defaultLineLineWidth',1.2);
set(0,'defaultLegendFontSize',16);
set(0,'defaultLegendFontName','times new roman');
set(0,'DefaultAxesXGrid','on');
set(0,'DefaultAxesYGrid','on');

%%% set user color set %%%
usr_color.paleblue = [204,242,255]/255;
usr_color.lightblue = [0 0.4470 0.7410];
usr_color.orange = [0.8500 0.3250 0.0980];


%% set parameters
setting_switch.read_outputcames = 1;  % 1: read from Outputcmaes.txt, 0: load from presaved .mat file
setting_switch.draw_itr = true;      % draw figure of CMA iteration process
setting_switch.compute_train = false; % compute obj. function of train data
setting_switch.compute_test = true;   % cpmpute obj. function of test data
setting_switch.draw_traj_fig = true;  % draw figure
setting_switch.draw_traj_video = false; % draw video 
CmaSettings.validation_period = 100.0;


%% principal perticulars (for drawing)
pp.lpp = 3.0d0;
pp.breadth = 0.48925d0;

% search outfiles directiory
outputfiles = dir('../case*');



% open output files
for cases = 1:length(outputfiles)
% for cases=2 % use this for sentense for debug
    CMAsetting_path = strcat(outputfiles(cases).folder,'/', outputfiles(cases).name, '/', 'settings.txt');
    CMAprob_path    = strcat(outputfiles(cases).folder,'/', outputfiles(cases).name, '/', 'InitProblem31.txt');
    outputcmaes_path    = strcat(outputfiles(cases).folder,'/', outputfiles(cases).name, '/', 'OutputCmaes.txt');
    
    postprocess_eachcase(setting_switch,usr_color, CmaSettings, outputfiles(cases), CMAsetting_path, CMAprob_path, outputcmaes_path, pp);
end
%


%% read output file of CMEA-ES: OutputCmaes.txt

% for i=1:5
%     tline = fgetl(fid_CMAsetting);
% 
%     tline = extractBefore( tline, "!" );
%     Data = str2double( tline );
%     if i==2
%         CmaSetting.switch_wind = Data;
%     elseif i==3
%         CmaSetting.type_si = Data;
%     elseif i==4
%         CmaSetting.number_files = Data;
%     elseif i==5
%         CmaSetting.integration_period = Data;
%     end
% end

% initprob = textscan(fid_CMAprob,'%c %c %c %d');
% dim  = cell2mat(initprob(4)); % Number of variables

% Data input and sort

% j=0;
% 
% %
% 
% tline = fgetl(fid_outputcmaes);
% while ischar(tline)
%     disp(tline)
%     tline = fgetl(fid_outputcmaes);
%     if ischar(tline)==0
%         break;
%     end
%     Data = str2num(tline);
%     % Find Infinity and NaN output. 
%     % In this case, utilize the previous values
%     TF1 = contains(tline,"Infinity");
%     TF2 = contains(tline,"NaN");
%     %
%     j=j+1;
%     %
%     if j==1
%         [rows,cols] = size(Data);
%         vecD=zeros(rows,dim);
%         xmean=zeros(rows,dim);
%         diagSqrtC=zeros(rows,dim);
%     end
%     %
%     if TF1==1 || TF2==1
%         itr    (j)=itr    (j-1);
%         arffeas(j)=arffeas(j-1);
%         arf    (j)=arf    (j-1);
%         df     (j)=df     (j-1);
%         dx     (j)=dx     (j-1);
%         sigma  (j)=sigma  (j-1);
%         %
%         for i=1:dim
%             vecD     (j,i)=vecD     (j-1,i);
%             xmean    (j,i)=xmean    (j-1,i);
%             diagSqrtC(j,i)=diagSqrtC(j-1,i);
%         end
%     else
%         itr    (j)=Data(1,1);
%         arffeas(j)=Data(1,2);
%         arf    (j)=Data(1,3);
%         df     (j)=Data(1,4);
%         dx     (j)=Data(1,5);
%         sigma  (j)=Data(1,6);
%         %
%         for i=1:dim
%             vecD     (j,i)=Data(1,6      +i);
%             xmean    (j,i)=Data(1,6+  dim+i);
%             diagSqrtC(j,i)=Data(1,6+2*dim+i);
%         end
%     end
%  end
% %
% % choose minimum J
% 
% [ minarffeas, minj]=min(arffeas);
% [rows,cols] = size(itr);
% 
% sigma_sqrt_eigC = zeros(cols, dim);
% sigma_sqrt_daigC= zeros(cols, dim);
% 
% for i=1:cols
%     sigma_sqrt_eigC (i,:)=sigma(i)*vecD(i,:);
%     sigma_sqrt_daigC(i,:)=sigma(i)*diagSqrtC(i,:);
% end
% %
% %   Set itrplot (horizontal axis of the graph)
% Npower=fix(log10(itr(cols)));
% itrplot=round(itr(cols),-Npower);
% if itrplot<itr(cols)
%     itrplot=itrplot+10^(Npower);
% end
% %
% fclose(fid_outputcmaes);
% Xopt=zeros(1,dim);
% for i=1:dim
%     Xopt(i)=xmean(minj,i);
% end

% save optimized_coeffs.mat
% 
% %% CMAES convergence output
% load optimized_coeffs.mat 
% % CmaSetting.validation_period = 100.0;
%%
 %% draw figures
%  % make minibatch
%  if CmaSettings.validation_type ==2
%     number_batch = fix(step_max( CmaSettings.number_testfiles)/reset_freq);
%  elseif CmaSettings.validation_type == 3
%      number_batch = 360;
%      reset_freq = length(state_input);
%  else
%      message = 'not sufficient setting'
%      pause
%  end
%  
% filetype = 'png';
% plot_time =  [0: reset_freq-1] * time_step_size;
  %% plot Obj
%   figure(20)
%   barplot=horzcat(transpose(Obj_test_batch(1,:)),transpose(Obj_ori_test_batch(1,:)));
%   bar(barplot)
%   hold on
% %   bar(transpose(Obj_ori_test_batch),'facealpha',0.5)
%   xlabel('batch no.')
%   ylabel('\itJ_1')
%   legend('CMA-ES Coeff.','EFD Coef.')
%   figure(21)
%   barplot=horzcat(transpose(Obj_test_batch(3,:)),transpose(Obj_ori_test_batch(3,:)));
%   bar(barplot)
%   hold on
% %   bar(transpose(Obj_ori_test_batch),'facealpha',0.5)
%   xlabel('batch no.')
%   ylabel('\itJ_3')
%   legend('CMA-ES Coeff.','EFD Coef.')
%   
%   figure(22)
%   barplot=horzcat(transpose(Obj_test_batch(2,:)),transpose(Obj_ori_test_batch(2,:)));
%   bar(barplot)
%   hold on
% %   bar(transpose(Obj_ori_test_batch),'facealpha',0.5)
%   xlabel('batch no.')
%   ylabel('\itJ_2')
%   legend('CMA-ES Coeff.','EFD Coef.')
%   fig = gcf;
%   figurename = strcat(condition_name,'_', inputfilename,'_DTeval=',num2str(CmaSetting.validation_period),'_J2');
%   SaveFig(fig,figurename)
 
    %%
%     %% animation of trajectory
% %     %   Prepare the figure block
%     time_video = time_input(:)-time_input(1);
%     S = get(0,'screensize');
%     figvideo = figure(11);
%     figvideo.Color='white';
%     figvideo.Position=[ S(1) S(2)  S(3) S(4)];
%     
%     videoname=strcat('./',condition_name,'_', inputfilename,'_video.mp4');
%     v = VideoWriter(videoname,'MPEG-4');
%     v.FrameRate = 30;
%     open(v);
%     
%     
%     %%%
%     minxaxis_t =min(min(state_input(:,1)),min(state_sim(:,1)))/pp.lpp-1;
%     minyaxis_t =min(min(state_input(:,3)),min(state_sim(:,3)))/pp.lpp-1;
%     maxxaxis_t =max(max(state_input(:,1)),max(state_sim(:,1)))/pp.lpp+1;
%     maxyaxis_t =max(max(state_input(:,3)),max(state_sim(:,3)))/pp.lpp+1;
%     %%%
%      % u velocity plot
%      
%      minyaxis =min(min(state_input(:,2)),min(state_sim_ori(:,2)))*1.1;
%      maxyaxis =max(max(state_input(:,2)),max(state_sim_ori(:,2)))*1.1;
%     
%      fig_u=subplot(4,2,5);
%      plot(time_video,state_input(:,2),'k-')
%      hold(fig_u,'on')
%      plot(time_video,state_sim(:,2),'-','color',lightblue)
%      plot(time_video,state_sim_ori(:,2),'-','color',orange)
%      p_u_input = plot(time_video(1),state_input(1,2),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
%      p_u_sim = plot(time_video(1),state_sim(1,2),'LineStyle','none','Marker','o','MarkerFaceColor',lightblue,'MarkerEdgeColor',lightblue,'MarkerSize',6);
%      p_u_sim_ori = plot(time_video(1),state_sim_ori(1,2),'LineStyle','none','Marker','o','MarkerFaceColor',orange,'MarkerEdgeColor',orange,'MarkerSize',6);
%      l_u_input = plot([time_video(1) time_video(1)], [-36 36],'k-');
%      
%      patch_time_u = patch(fig_u, [time_video(1) time_video(1) time_video(1) time_video(1)],[-36 -36 36 36],paleblue,'facealpha',0.2);
%      axis(fig_u,[time_video(1) time_video(end) minyaxis maxyaxis]);
%      set(gca,'Xgrid','on','Ygrid','on');
%      xtickformat('%.0f')
%      ytickformat('%.1f')
%      set(gca,'XMinorTick','on','YMinorTick','on');
%      xlabel('Time \rm[s\rm]' );
%      ylabel('\itu \rm[m/s]')
%      hold(fig_u, 'off')
%      
%       % vm velocity plot
%       
%      minyaxis =min(min(state_input(:,4)),min(state_sim_ori(:,4)))*1.1;
%      maxyaxis =max(max(state_input(:,4)),max(state_sim_ori(:,4)))*1.1;
%     
%      fig_vm=subplot(4,2,6);
%      plot(time_video,state_input(:,4),'k-')
%      hold(fig_vm,'on')
%      plot(time_video,state_sim(:,4),'-','color',lightblue)
%      plot(time_video,state_sim_ori(:,4),'-','color',orange)
%      p_vm_input = plot(time_video(1),state_input(1,4),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
%      p_vm_sim = plot(time_video(1),state_sim(1,4),'LineStyle','none','Marker','o','MarkerFaceColor',lightblue,'MarkerEdgeColor',lightblue,'MarkerSize',6);
%      p_vm_sim_ori = plot(time_video(1),state_sim_ori(1,4),'LineStyle','none','Marker','o','MarkerFaceColor',orange,'MarkerEdgeColor',orange,'MarkerSize',6);
%      l_vm_input = plot([time_video(1) time_video(1)], [-36 36],'k-');
%      
%      patch_time_vm = patch(fig_vm, [time_video(1) time_video(1) time_video(1) time_video(1)],[-36 -36 36 36],paleblue,'facealpha',0.2);
%      axis(fig_vm,[time_video(1) time_video(end) minyaxis maxyaxis]);
%      set(gca,'Xgrid','on','Ygrid','on');
%      xtickformat('%.0f')
%      ytickformat('%.1f')
%      set(gca,'XMinorTick','on','YMinorTick','on');
%      xlabel('Time \rm[s\rm]' );
%      ylabel('\itv_m \rm[m/s]')
%      hold(fig_vm, 'off')
%      
%      % r plot
%      
%      minyaxis = rad2deg(min(min(state_input(:,6)),min(state_sim(:,6))))*1.1;
%      maxyaxis = rad2deg(max(max(state_input(:,6)),max(state_sim(:,6))))*1.1;
%      
%      fig_r=subplot(4,2,7);
%      plot(time_video,rad2deg(state_input(:,6)),'k-')
%      hold(fig_r,'on')
%      plot(time_video,rad2deg(state_sim(:,6)),'-','color',lightblue)
%      plot(time_video,rad2deg(state_sim_ori(:,6)),'-','color',orange)
%      p_r_input = plot(time_video(1),rad2deg(state_input(1,6)),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
%      p_r_sim = plot(time_video(1),rad2deg(state_sim(1,6)),'o','LineStyle','none','Marker','o','MarkerFaceColor',lightblue,'MarkerEdgeColor',lightblue,'MarkerSize',6);
%      p_r_sim_ori = plot(time_video(1),rad2deg(state_sim_ori(1,6)),'LineStyle','none','Marker','o','MarkerFaceColor',orange,'MarkerEdgeColor',orange,'MarkerSize',6);
%      l_r_input = plot([time_video(1) time_video(1)], [-36 36],'k-');
%      
%      patch_time_r = patch(fig_r, [time_video(1) time_video(1) time_video(1) time_video(1)],[-36 -36 36 36],paleblue,'facealpha',0.2);
%      axis(fig_r,[time_video(1) time_video(end) minyaxis maxyaxis]);
%      set(gca,'Xgrid','on','Ygrid','on');
%      xtickformat('%.0f')
%      ytickformat('%.1f')
%      set(gca,'XMinorTick','on','YMinorTick','on');
%      xlabel('Time \rm[s\rm]' );
%      ylabel('\itr \rm[degree/s]')
%      hold(fig_r, 'off')
%      
%     
%      % control plot
% 
%      
%      fig_control=subplot(4,2,8);
%      plot(time_video,rad2deg(delta_input(:)),'k-')
%      hold(fig_control,'on')
%      plot(time_video,n_input,'-','color',lightblue)
%      p_delta = plot(time_video(1),rad2deg(delta_input(1)),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
%      p_n = plot(time_video(1),n_input(1),'o','LineStyle','none','Marker','o','MarkerFaceColor',lightblue,'MarkerEdgeColor',lightblue,'MarkerSize',6);
%      
%      patch_time_control = patch(fig_control, [time_video(1) time_video(1) time_video(1) time_video(1)],[-36 -36 36 36],paleblue,'facealpha',0.2);
%      l_control = plot([time_video(1) time_video(1)], [-36 36],'k-');
%      axis(fig_control,[time_video(1) time_video(end) -36 36]);
%      set(gca,'Xgrid','on','Ygrid','on');
%      xtickformat('%.0f')
%      ytickformat('%.1f')
%      set(gca,'XMinorTick','on','YMinorTick','on');
%      xlabel('Time \rm[s\rm]' );
%      ylabel('\delta \rm[degree], \itn \rm[rps]')
%      legend('\delta', '\itn','location','northeastoutside')
%      hold(fig_control, 'off')
%      
%      % trajctory anime 
%         for i=1:step_max( CmaSetting.number_testfiles)-start_time_step+1
%             if i==1 || mod(i, 5)==1 || i==step_max( CmaSetting.number_testfiles)
%             
%             
%             fig_traj=subplot(4,2,[1, 2, 3, 4]);
%             plot(state_input(1:i,1)/pp.lpp,state_input(1:i,3)/pp.lpp,'k--' )
%             hold(fig_traj,'on')
%             plot(fig_traj, state_sim(1:i,1)/pp.lpp,state_sim(1:i,3)/pp.lpp,'linestyle','--','color',lightblue )
%              plot(fig_traj, state_sim_ori(1:i,1)/pp.lpp,state_sim_ori(1:i,3)/pp.lpp,'linestyle','--','color',orange )
%              
%             % plot input ship
%             Y = [state_input( i, 3 ) - pp.breadth / 2.0d0, state_input( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_input( i, 3 ) + pp.breadth / 2.0d0, state_input( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_input( i, 3 ) + pp.breadth / 2.0d0, state_input( i, 3 );...
%                  state_input( i, 3 ), state_input( i, 3 ) - pp.breadth / 2.0d0;...
%                  state_input( i, 3 ) - pp.breadth / 2.0d0, state_input( i, 3 ) - pp.breadth / 2.0d0 ]/pp.lpp;
%             X = [state_input( i, 1 ) - pp.lpp / 2.0d0, state_input( i, 1 ) - pp.lpp / 2.0d0;...
%                  state_input( i, 1 ) - pp.lpp / 2.0d0, state_input( i, 1 ) + 1.0d0;...
%                  state_input( i, 1 ) + 1.0d0, state_input( i, 1 ) + pp.lpp / 2.0d0;...
%                  state_input( i, 1 ) + pp.lpp / 2.0d0, state_input( i, 1 ) + 1.0d0;...
%                  state_input( i, 1 ) + 1.0d0, state_input( i, 1 ) - pp.lpp / 2.0d0]/pp.lpp;
%             Cy = state_input( i, 3 )/pp.lpp;
%             Cx = state_input( i, 1 )/pp.lpp;
%             rot =  state_input(i,5);
%             %   Rotation of the coordinate
%             x2 = cos(rot) * X - sin(rot) * Y + Cx - Cx * cos(rot) + Cy * sin(rot);
%             y2 = sin(rot) * X + cos(rot) * Y + Cy - Cx * sin(rot) - Cy * cos(rot);           
% %             plot( y2 , x2 , 'k-', 'Linewidth', 0.5 );
%             plot( fig_traj, x2 , y2 , 'k-', 'Linewidth', 0.5 );
%             
%             % apparent wind of input data
% %             windend = [state_input(i, 3)/pp.lpp, state_input(i, 1)/pp.lpp];
% 
%             windend = [state_input(i, 1)/pp.lpp, state_input(i, 3)/pp.lpp];
%             windstart = [1 1];
%             arrowratio = 0.2;
%             windstart(1) = windend(1) + arrowratio* apparent_wind_velo(i) * sin(apparent_wind_dir(i) +state_input(i, 5));
%             windstart(2) = windend(2) + arrowratio* apparent_wind_velo(i) * cos(apparent_wind_dir(i) +state_input(i, 5));
%             % plot wind
%             quiver(windstart(1), windstart(2), windend(1)-windstart(1), windend(2)-windstart(2), 0, 'color','r' ,'linewidth',2.0,...
%                 'MaxHeadSize',1.5)
% 
%             %plot sim ship
%             Y = [state_sim( i, 3 ) - pp.breadth / 2.0d0, state_sim( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_sim( i, 3 ) + pp.breadth / 2.0d0, state_sim( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_sim( i, 3 ) + pp.breadth / 2.0d0, state_sim( i, 3 );...
%                  state_sim( i, 3 ), state_sim( i, 3 ) - pp.breadth / 2.0d0;...
%                  state_sim( i, 3 ) - pp.breadth / 2.0d0, state_sim( i, 3 ) - pp.breadth / 2.0d0 ]/pp.lpp;
%             X = [state_sim( i, 1 ) - pp.lpp / 2.0d0, state_sim( i, 1 ) - pp.lpp / 2.0d0;...
%                  state_sim( i, 1 ) - pp.lpp / 2.0d0, state_sim( i, 1 ) + 1.0d0;...
%                  state_sim( i, 1 ) + 1.0d0, state_sim( i, 1 ) + pp.lpp / 2.0d0;...
%                  state_sim( i, 1 ) + pp.lpp / 2.0d0, state_sim( i, 1 ) + 1.0d0;...
%                  state_sim( i, 1 ) + 1.0d0, state_sim( i, 1 ) - pp.lpp / 2.0d0]/pp.lpp;
%             Cy = state_sim( i, 3 )/pp.lpp;
%             Cx = state_sim( i, 1 )/pp.lpp;
%             rot =  state_sim(i,5);
%             %   Rotation of the coordinate
%             x2 = cos(rot) * X - sin(rot) * Y + Cx - Cx * cos(rot) + Cy * sin(rot);
%             y2 = sin(rot) * X + cos(rot) * Y + Cy - Cx * sin(rot) - Cy * cos(rot); 
%             
% %             plot( y2 , x2 ,'color',lightblue,'linestyle','-','linewidth',0.5)
%              plot( x2 , y2 ,'color',lightblue,'linestyle','-','linewidth',0.5)
%             %plot EFD Coef. ship
%              Y = [state_sim_ori( i, 3 ) - pp.breadth / 2.0d0, state_sim_ori( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_sim_ori( i, 3 ) + pp.breadth / 2.0d0, state_sim_ori( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_sim_ori( i, 3 ) + pp.breadth / 2.0d0, state_sim_ori( i, 3 );...
%                  state_sim_ori( i, 3 ), state_sim_ori( i, 3 ) - pp.breadth / 2.0d0;...
%                  state_sim_ori( i, 3 ) - pp.breadth / 2.0d0, state_sim_ori( i, 3 ) - pp.breadth / 2.0d0 ]/pp.lpp;
%             X = [state_sim_ori( i, 1 ) - pp.lpp / 2.0d0, state_sim_ori( i, 1 ) - pp.lpp / 2.0d0;...
%                  state_sim_ori( i, 1 ) - pp.lpp / 2.0d0, state_sim_ori( i, 1 ) + 1.0d0;...
%                  state_sim_ori( i, 1 ) + 1.0d0, state_sim_ori( i, 1 ) + pp.lpp / 2.0d0;...
%                  state_sim_ori( i, 1 ) + pp.lpp / 2.0d0, state_sim_ori( i, 1 ) + 1.0d0;...
%                  state_sim_ori( i, 1 ) + 1.0d0, state_sim_ori( i, 1 ) - pp.lpp / 2.0d0]/pp.lpp;
%             Cy = state_sim_ori( i, 3 )/pp.lpp;
%             Cx = state_sim_ori( i, 1 )/pp.lpp;
%             rot =  state_sim_ori(i,5);
%             %   Rotation of the coordinate
%             x2 = cos(rot) * X - sin(rot) * Y + Cx - Cx * cos(rot) + Cy * sin(rot);
%             y2 = sin(rot) * X + cos(rot) * Y + Cy - Cx * sin(rot) - Cy * cos(rot);           
% %             plot( y2 , x2 ,'color',orange,'linestyle','-','linewidth',0.5)
%             plot( x2 , y2 ,'color',orange,'linestyle','-','linewidth',0.5)
%             
%             axis equal
%             axis ij
%             axis([minxaxis_t, maxxaxis_t, minyaxis_t, maxyaxis_t])
%             ylabel ('$$\hat{y}/L_{pp}$$ ','interpreter','latex')
%             xlabel ('$$\hat{x}/L_{pp}$$ ','interpreter','latex')
%             legend('Input','CMA-ES Coef.', 'EFD Coef.','location','northeastoutside')
%             
%             % indecate time
%             refpointx = 2.0;
%             refpointy = -2.0;
% 
%             formatSpec = 'time = %.3f s';
%             str=sprintf( formatSpec, time_video(i));
%             text(refpointx, refpointy ,str, 'FontName','Times','fontsize',16);
%             
%             % plot current control value
%             p_u_input.XData = time_video(i);
%             p_u_input.YData = state_input(i,2);
%             
%             p_u_sim.XData = time_video(i);
%             p_u_sim.YData = state_sim(i,2);
%             
%             p_u_sim_ori.XData = time_video(i);
%             p_u_sim_ori.YData = state_sim_ori(i,2);
%             
%             p_vm_input.XData = time_video(i);
%             p_vm_input.YData = state_input(i,4);
%             
%             p_vm_sim.XData = time_video(i);
%             p_vm_sim.YData = state_sim(i,4);
%             
%             p_vm_sim_ori.XData = time_video(i);
%             p_vm_sim_ori.YData = state_sim_ori(i,4);
%             
%             p_r_input.XData = time_video(i);
%             p_r_input.YData = rad2deg(state_input(i,6));
%             
%             p_r_sim.XData = time_video(i);
%             p_r_sim.YData = rad2deg(state_sim(i,6));
%             
%             p_r_sim_ori.XData = time_video(i);
%             p_r_sim_ori.YData = rad2deg(state_sim_ori(i,6));
%             
%             p_n.XData = time_video(i);
%             p_n.YData = n_input(i);
% 
%             p_delta.XData = time_video(i);
%             p_delta.YData = rad2deg(delta_input(i));
%             
%             % plot time bar
%             l_u_input.XData = [time_video(i) time_video(i)];
%             l_u_input.YData = [-36 36];
% 
%             l_vm_input.XData = [time_video(i) time_video(i)];
%             l_vm_input.YData = [-36 36];
% 
%             l_r_input.XData = [time_video(i) time_video(i)];
%             l_r_input.YData = [-36 36];
%             
%             l_control.XData = [time_video(i) time_video(i)];
%             l_control.YData = [-36 36];
%             
%             % plot time area
%             patch_time_u.XData = [time_video(1) time_video(i) time_video(i) time_video(1)];
%             patch_time_u.YData = [-36 -36 36 36];
%             
%             patch_time_vm.XData = [time_video(1) time_video(i) time_video(i) time_video(1)];
%             patch_time_vm.YData = [-36 -36 36 36];
%             
%             patch_time_r.XData = [time_video(1) time_video(i) time_video(i) time_video(1)];
%             patch_time_r.YData = [-36 -36 36 36];
%             
%             patch_time_control.XData = [time_video(1) time_video(i) time_video(i) time_video(1)];
%             patch_time_control.YData = [-36 -36 36 36];
%             drawnow
%             
%             % save video
%             movmov = getframe( gcf );
%             writeVideo( v, movmov );
% 
%             hold ( fig_traj, 'off' );
%             plot( 0 )
% 
%             end
%         end
%      
%      close(v);
    %%
%%
%  for batch = 1:number_batch
% %  for batch = 6
%     close all
%     %%%  state variables
%     figure(2)
% 
%     lot = [1+reset_freq*(batch-1):reset_freq*batch];
%     
%     plot(plot_time, state_input(lot,2 ),'k-')
%     hold on
%     plot(plot_time, state_sim(lot, 2),'color',lightblue,'linestyle','--')
%     plot(plot_time, state_sim_ori(lot, 2),'color',orange,'linestyle','--')
%     xlabel('\it t \rm(\its\rm)')
%     ylabel('\itu \rm(m/s)')
%     legend('Input','CMA-ES Coeff.', 'EFD Coef.','location','best')
%     xlim([0 plot_time(end)+1])
%     hold off
%     %%%
% %     ax = gca;
% %     outerpos = ax.OuterPosition;
% %     ti = ax.TightInset; 
% %     left = outerpos(1) + ti(1);
% %     bottom = outerpos(2) + ti(2);
% %     ax_width = outerpos(3) - ti(1) - ti(3);
% %     ax_height = outerpos(4) - ti(2) - ti(4);
% %     ax.Position = [left bottom ax_width ax_height];
%     
%     fig = gcf;
%      figurename=strcat('./',inputfilename,'_',num2str(batch),'_u.',filetype);
%     saveas(gca,figurename,filetype)
%     figurename_pdf=strcat('./',inputfilename,'_',num2str(batch),'_u.pdf');
%     print(fig,figurename_pdf,'-dpdf')
%     
%     figure(3)
%     plot(plot_time, state_input(lot,4 ),'k-')
%     hold on
%     plot(plot_time, state_sim(lot, 4),'color',lightblue,'linestyle','--')
%     plot(plot_time, state_sim_ori(lot, 4),'color',orange,'linestyle','--')
%     xlabel('\it t \rm(\its\rm)')
%     ylabel('\itv_{m} \rm(m/s)')
% %     legend('Input','CMA-ES Coeff.', 'EFD Coef.','location','best')
%     xlim([0 plot_time(end)+1])
%     hold off
%          %%%
% %     ax = gca;
% %     outerpos = ax.OuterPosition;
% %     ti = ax.TightInset; 
% %     left = outerpos(1) + ti(1);
% %     bottom = outerpos(2) + ti(2);
% %     ax_width = outerpos(3) - ti(1) - ti(3);
% %     ax_height = outerpos(4) - ti(2) - ti(4);
% %     ax.Position = [left bottom ax_width ax_height];
% %     
%     fig = gcf;
%      figurename=strcat('./',inputfilename,'_',num2str(batch),'_vm.',filetype);
%     saveas(gca,figurename,filetype)
%     figurename_pdf=strcat('./',inputfilename,'_',num2str(batch),'_vm.pdf');
%     print(fig,figurename_pdf,'-dpdf')
%     
%     
%     figure(4)
%     lot = [1+reset_freq*(batch-1):reset_freq*batch];
%     plot(plot_time, state_input(lot,6 ),'k-')
%     hold on
%     plot(plot_time, state_sim(lot, 6),'color',lightblue,'linestyle','--')
%     plot(plot_time, state_sim_ori(lot, 6),'color',orange,'linestyle','--')
%     xlabel('\it t \rm(\its\rm)')
%     ylabel('\itr \rm(rad/s)')
% %     legend('Input','CMA-ES Coeff.', 'EFD Coef.','location','best')
%     xlim([0 plot_time(end)+1])
%     hold off
%     %%%
% %     ax = gca;
% %     outerpos = ax.OuterPosition;
% %     ti = ax.TightInset; 
% %     left = outerpos(1) + ti(1);
% %     bottom = outerpos(2) + ti(2);
% %     ax_width = outerpos(3) - ti(1) - ti(3);
% %     ax_height = outerpos(4) - ti(2) - ti(4);
% %     ax.Position = [left bottom ax_width ax_height];
%     
%     fig = gcf;
%      figurename=strcat('./',inputfilename,'_',num2str(batch),'_r.',filetype);
%     saveas(gca,figurename,filetype)
%     figurename_pdf=strcat('./',inputfilename,'_',num2str(batch),'_r.pdf');
%     print(fig,figurename_pdf,'-dpdf')
%     
%     figure(5)
%     plot(plot_time, n_input(lot),'k-')
%     hold on
%     plot(plot_time, rad2deg(delta_input(lot)),'k--')
%     xlabel('\it t \rm(\its\rm)')
%     ylabel('\delta (deg), \it n_p \rm(rps)')
%     legend('\it n_p','\delta','location','best')
%     xlim([0 plot_time(end)+1])
%     ylim([-36 36])
%     hold off
% 
%          %%%
% %     ax = gca;
% %     outerpos = ax.OuterPosition;
% %     ti = ax.TightInset; 
% %     left = outerpos(1) + ti(1);
% %     bottom = outerpos(2) + ti(2);
% %     ax_width = outerpos(3) - ti(1) - ti(3);
% %     ax_height = outerpos(4) - ti(2) - ti(4);
% %     ax.Position = [left bottom ax_width ax_height];
% %     
%     fig = gcf;
%      figurename=strcat('./',inputfilename,'_',num2str(batch),'_control.',filetype);
%     saveas(gca,figurename,filetype)
%     figurename_pdf=strcat('./',inputfilename,'_',num2str(batch),'_control.pdf');
%     print(fig,figurename_pdf,'-dpdf')
%     
%     
%     % set t=0 x,y = 0
%     state_input_traj = zeros(length(state_input),6);
%     state_sim_traj = zeros(length(state_input),6);
%     state_sim_ori_traj = zeros(length(state_input),6);
%     
%     for i= 1:6
%         if (i==1 || i==3)
%             state_input_traj(:,i) = state_input(:,i)-state_input(lot(1),i);
%             state_sim_traj(:,i) = state_sim(:,i)-state_sim(lot(1),i);
%             state_sim_ori_traj(:,i) = state_sim_ori(:,i)-state_sim_ori(lot(1),i);
%         else
%             state_input_traj(:,i) = state_input(:,i);
%             state_sim_traj(:,i) = state_sim(:,i);
%             state_sim_ori_traj(:,i) = state_sim_ori(:,i);
%         end
%     end
%     
%         % psi to [-pi pi]
%     for j=1:length(state_sim_traj)
%         if state_sim_traj(j,5) > pi
%             state_sim_traj(j,5) = state_sim_traj(j,5)-2*pi; 
%         elseif state_sim_traj(j,5) < -pi
%             state_sim_traj(j,5) = state_sim_traj(j,5)+2*pi;
%         end
%         
%         if state_sim_ori_traj(j,5) > pi
%             state_sim_ori_traj(j,5) = state_sim_ori_traj(j,5)-2*pi; 
%         elseif state_sim_ori_traj(j,5) < -pi
%             state_sim_ori_traj(j,5) = state_sim_ori_traj(j,5)+2*pi;
%         end
%         
%     end
%     
%     %%%  trajectory time history
%     figure('Color','white','Position',[2*S(4)/12 2*S(4)/12 5*S(3)/12 9*S(4)/12])
%     subplot(2,2,1);
%     lot = [1+reset_freq*(batch-1):reset_freq*batch];
%     plot(plot_time, state_input_traj(lot,1 )/pp.lpp,'k-')
%     hold on
%     plot(plot_time, state_sim_traj(lot, 1)/pp.lpp,'color',lightblue,'linestyle','--')
%     plot(plot_time, state_sim_ori_traj(lot, 1)/pp.lpp,'color',orange,'linestyle','--')
%     xlabel('\it t \rm(\its\rm)')
%     ylabel ('$$\hat{x}/L_{pp}$$ ','interpreter','latex')
%     legend('Input','CMA-ES Coeff.', 'EFD Coef.','location','best')
%     xlim([0 plot_time(end)+1])
%     hold off
%         
%     subplot(2,2,2);
%     lot = [1+reset_freq*(batch-1):reset_freq*batch];
%     plot(plot_time, state_input_traj(lot,3 )/pp.lpp,'k-')
%     hold on
%     plot(plot_time, state_sim_traj(lot, 3)/pp.lpp,'color',lightblue,'linestyle','--')
%     plot(plot_time, state_sim_ori_traj(lot, 3)/pp.lpp,'color',orange,'linestyle','--')
% 
%     xlabel('\it t \rm(\its\rm)')
%     ylabel ('$$\hat{y}/L_{pp}$$ ','interpreter','latex')
%     legend('Input','CMA-ES Coeff.', 'EFD Coef.','location','best')
%     xlim([0 plot_time(end)+1])
%     hold off
%       
%     subplot(2,2,3);
%     lot = [1+reset_freq*(batch-1):reset_freq*batch];
%     plot(plot_time, state_input_traj(lot,5 ),'k-')
%     hold on
%     plot(plot_time, state_sim_traj(lot, 5),'color',lightblue,'linestyle','--')
%     plot(plot_time, state_sim_ori_traj(lot, 5),'color',orange,'linestyle','--')
%     xlabel('\it t \rm(\its\rm)')
%     ylabel ('$$\hat{\psi}$$ ','interpreter','latex')
%     legend('Input','CMA-ES Coeff.', 'EFD Coef.','location','best')
%     xlim([0 plot_time(end)+1])
%     hold off
%     
% 
%     subplot(2,2,4)
%     yyaxis left
%     plot(plot_time, apparent_wind_velo(lot),'k-')
%     hold on
%     yyaxis right
%     plot(plot_time, rad2deg(apparent_wind_dir(lot)),'k.')
%     xlabel('\it t \rm(\its\rm)')
%     yyaxis left
%     ylabel('\it U_A \rm(m/s)')
%     ylim([0 6])
%     yyaxis right
%     ylabel('\gamma_A (degree)')
%     legend('\it U_A','\gamma_A','location','best')
%     xlim([0 plot_time(end)+1])
%     ax = gca;
%     ax.YAxis(1).Color = 'k';
%     ax.YAxis(2).Color = 'k';
%     hold off
%     figurename=strcat('./',inputfilename,'_',num2str(batch),'_trajhist.',filetype);
%     saveas(gca,figurename,filetype)
%     fig = gcf;
%     figurename_pdf=strcat('./',inputfilename,'_',num2str(batch),'_trajhist.pdf');
%     print(fig,figurename_pdf,'-dpdf')
% 
%     figure(8)
%     yyaxis left
%     plot(plot_time, apparent_wind_velo(lot),'k-')
%     hold on
%     yyaxis right
%     plot(plot_time, rad2deg(apparent_wind_dir(lot)),'k.')
%     xlabel('\it t \rm(\its\rm)')
%     yyaxis left
%     ylabel('\it U_A \rm(m/s)')
%     ylim([0 6])
%     yyaxis right
%     ylabel('\gamma_A (degree)')
%     legend('\it U_A','\gamma_A','location','best')
%     xlim([0 plot_time(end)+1])
%     ax = gca;
%     ax.YAxis(1).Color = 'k';
%     ax.YAxis(2).Color = 'k';
%     hold off
%     figurename=strcat('./',inputfilename,'_',num2str(batch),'_wind.',filetype);
%     saveas(gca,figurename,filetype)
%     fig = gcf;
%     figurename_pdf=strcat('./',inputfilename,'_',num2str(batch),'_wind.pdf');
%     print(fig,figurename_pdf,'-dpdf')
%     
%     %%% trajectory 
%     figure(7)
%     plot(state_input_traj(lot, 3)/pp.lpp, state_input_traj(lot, 1)/pp.lpp,'k:','linewidth',1.5)
%     hold on
%     plot(state_sim_traj(lot, 3)/pp.lpp, state_sim_traj(lot, 1)/pp.lpp,'color',lightblue,'linestyle',':','linewidth',1.5)
%     plot(state_sim_ori_traj(lot, 3)/pp.lpp, state_sim_ori_traj(lot, 1)/pp.lpp,'color',orange,'linestyle',':','linewidth',1.5)
%         for i=lot(1):lot(end)
%             if i==lot(1) || mod(i, 200)==0 || i==lot(end)
%             % plot input ship
%             Y = [state_input_traj( i, 3 ) - pp.breadth / 2.0d0, state_input_traj( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_input_traj( i, 3 ) + pp.breadth / 2.0d0, state_input_traj( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_input_traj( i, 3 ) + pp.breadth / 2.0d0, state_input_traj( i, 3 );...
%                  state_input_traj( i, 3 ), state_input_traj( i, 3 ) - pp.breadth / 2.0d0;...
%                  state_input_traj( i, 3 ) - pp.breadth / 2.0d0, state_input_traj( i, 3 ) - pp.breadth / 2.0d0 ]/pp.lpp;
%             X = [state_input_traj( i, 1 ) - pp.lpp / 2.0d0, state_input_traj( i, 1 ) - pp.lpp / 2.0d0;...
%                  state_input_traj( i, 1 ) - pp.lpp / 2.0d0, state_input_traj( i, 1 ) + 1.0d0;...
%                  state_input_traj( i, 1 ) + 1.0d0, state_input_traj( i, 1 ) + pp.lpp / 2.0d0;...
%                  state_input_traj( i, 1 ) + pp.lpp / 2.0d0, state_input_traj( i, 1 ) + 1.0d0;...
%                  state_input_traj( i, 1 ) + 1.0d0, state_input_traj( i, 1 ) - pp.lpp / 2.0d0]/pp.lpp;
%             Cy = state_input_traj( i, 3 )/pp.lpp;
%             Cx = state_input_traj( i, 1 )/pp.lpp;
%             rot =  state_input_traj(i,5);
%             %   Rotation of the coordinate
%             x2 = cos(rot) * X - sin(rot) * Y + Cx - Cx * cos(rot) + Cy * sin(rot);
%             y2 = sin(rot) * X + cos(rot) * Y + Cy - Cx * sin(rot) - Cy * cos(rot);           
%             plot( y2 , x2 , 'k-', 'Linewidth', 0.5 );
%             
%             % apparent wind of input data
%             windend = [state_input_traj(i, 3)/pp.lpp, state_input_traj(i, 1)/pp.lpp];
%             windstart = [1 1];
%             arrowratio = 0.2;
%             windstart(1) = windend(1) + arrowratio* apparent_wind_velo(i) * sin(apparent_wind_dir(i) +state_input_traj(i, 5));
%             windstart(2) = windend(2) + arrowratio* apparent_wind_velo(i) * cos(apparent_wind_dir(i) +state_input_traj(i, 5));
%             % plot wind
%             quiver(windstart(1), windstart(2), windend(1)-windstart(1), windend(2)-windstart(2), 0, 'color','r' ,'linewidth',2.0,...
%                 'MaxHeadSize',1.5)
% 
%             %plot sim ship
%             Y = [state_sim_traj( i, 3 ) - pp.breadth / 2.0d0, state_sim_traj( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_sim_traj( i, 3 ) + pp.breadth / 2.0d0, state_sim_traj( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_sim_traj( i, 3 ) + pp.breadth / 2.0d0, state_sim_traj( i, 3 );...
%                  state_sim_traj( i, 3 ), state_sim_traj( i, 3 ) - pp.breadth / 2.0d0;...
%                  state_sim_traj( i, 3 ) - pp.breadth / 2.0d0, state_sim_traj( i, 3 ) - pp.breadth / 2.0d0 ]/pp.lpp;
%             X = [state_sim_traj( i, 1 ) - pp.lpp / 2.0d0, state_sim_traj( i, 1 ) - pp.lpp / 2.0d0;...
%                  state_sim_traj( i, 1 ) - pp.lpp / 2.0d0, state_sim_traj( i, 1 ) + 1.0d0;...
%                  state_sim_traj( i, 1 ) + 1.0d0, state_sim_traj( i, 1 ) + pp.lpp / 2.0d0;...
%                  state_sim_traj( i, 1 ) + pp.lpp / 2.0d0, state_sim_traj( i, 1 ) + 1.0d0;...
%                  state_sim_traj( i, 1 ) + 1.0d0, state_sim_traj( i, 1 ) - pp.lpp / 2.0d0]/pp.lpp;
%             Cy = state_sim_traj( i, 3 )/pp.lpp;
%             Cx = state_sim_traj( i, 1 )/pp.lpp;
%             rot =  state_sim_traj(i,5);
%             %   Rotation of the coordinate
%             x2 = cos(rot) * X - sin(rot) * Y + Cx - Cx * cos(rot) + Cy * sin(rot);
%             y2 = sin(rot) * X + cos(rot) * Y + Cy - Cx * sin(rot) - Cy * cos(rot);           
%             plot( y2 , x2 ,'color',lightblue,'linestyle','-','linewidth',0.5)
%             
%             %plot EFD Coef. ship
%              Y = [state_sim_ori_traj( i, 3 ) - pp.breadth / 2.0d0, state_sim_ori_traj( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_sim_ori_traj( i, 3 ) + pp.breadth / 2.0d0, state_sim_ori_traj( i, 3 ) + pp.breadth / 2.0d0;...
%                  state_sim_ori_traj( i, 3 ) + pp.breadth / 2.0d0, state_sim_ori_traj( i, 3 );...
%                  state_sim_ori_traj( i, 3 ), state_sim_ori_traj( i, 3 ) - pp.breadth / 2.0d0;...
%                  state_sim_ori_traj( i, 3 ) - pp.breadth / 2.0d0, state_sim_ori_traj( i, 3 ) - pp.breadth / 2.0d0 ]/pp.lpp;
%             X = [state_sim_ori_traj( i, 1 ) - pp.lpp / 2.0d0, state_sim_ori_traj( i, 1 ) - pp.lpp / 2.0d0;...
%                  state_sim_ori_traj( i, 1 ) - pp.lpp / 2.0d0, state_sim_ori_traj( i, 1 ) + 1.0d0;...
%                  state_sim_ori_traj( i, 1 ) + 1.0d0, state_sim_ori_traj( i, 1 ) + pp.lpp / 2.0d0;...
%                  state_sim_ori_traj( i, 1 ) + pp.lpp / 2.0d0, state_sim_ori_traj( i, 1 ) + 1.0d0;...
%                  state_sim_ori_traj( i, 1 ) + 1.0d0, state_sim_ori_traj( i, 1 ) - pp.lpp / 2.0d0]/pp.lpp;
%             Cy = state_sim_ori_traj( i, 3 )/pp.lpp;
%             Cx = state_sim_ori_traj( i, 1 )/pp.lpp;
%             rot =  state_sim_ori_traj(i,5);
%             %   Rotation of the coordinate
%             x2 = cos(rot) * X - sin(rot) * Y + Cx - Cx * cos(rot) + Cy * sin(rot);
%             y2 = sin(rot) * X + cos(rot) * Y + Cy - Cx * sin(rot) - Cy * cos(rot);           
%             plot( y2 , x2 ,'color',orange,'linestyle','-','linewidth',0.5)
%             end
%         end
%     axis equal
%     xlabel ('$$\hat{y}/L_{pp}$$ ','interpreter','latex')
%     ylabel ('$$\hat{x}/L_{pp}$$ ','interpreter','latex')
% %     legend('Input','CMA-ES Coeff.', 'EFD Coef','location','best','fontsize',12)
%     legend('Input','CMA-ES Coeff.', 'EFD Coef.','location','best')
%     yl =  ylim;
%     xl =  xlim;
%     xlim([xl(1)-0.5, xl(2)+0.5])
%     ylim([yl(1)-0.5, yl(2)+0.5])
%     hold off
%         %%%
% %     ax = gca;
% %     outerpos = ax.OuterPosition;
% %     ti = ax.TightInset*1.1; 
% %     left = outerpos(1) + ti(1);
% %     bottom = outerpos(2) + ti(2);
% %     ax_width = outerpos(3) - ti(1) - ti(3);
% %     ax_height = outerpos(4) - ti(2) - ti(4);
% %     ax.Position = [left bottom ax_width ax_height];
%     fig = gcf;
% 
%     figurename=strcat('./',inputfilename,'_',num2str(batch),'_traj.',filetype);
%     saveas(gca,figurename,filetype)
%     
%     figurename_pdf=strcat('./',inputfilename,'_',num2str(batch),'_traj.pdf');
%     print(fig,figurename_pdf,'-dpdf')
%     
%  end
 
%% Function goes here
%







