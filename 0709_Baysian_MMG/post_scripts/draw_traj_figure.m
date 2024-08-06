function draw_traj_figure(save_directory_path, condition_name, inputfilename, start_time_step, usr_color, CmaSettings,...
                          state_input, state_sim, state_sim_ori, n_input, delta_input, ...
                          step_max, apparent_wind_velo, apparent_wind_dir, time_input,pp )
    %%% trajectory of whole test data
 
    figtrajall=figure(10);
    time_trajall = time_input(:)-time_input(1);
    S = get(0,'screensize');
    figtrajall.Color='white';
    figtrajall.Position=[ S(1) S(2)  S(3) S(4)];
    
     %%%
     minxaxis_t =min(min(state_input(:,1)),min(state_sim(:,1)))/pp.lpp-1;
     minyaxis_t =min(min(state_input(:,3)),min(state_sim(:,3)))/pp.lpp-1;
     maxxaxis_t =max(max(state_input(:,1)),max(state_sim(:,1)))/pp.lpp+1;
     maxyaxis_t =max(max(state_input(:,3)),max(state_sim(:,3)))/pp.lpp+1;
     %%%
     % u velocity plot
     
     minyaxis =min(min(state_input(:,2)),min(state_sim_ori(:,2)))*1.1;
     maxyaxis =max(max(state_input(:,2)),max(state_sim_ori(:,2)))*1.1;
    
     fig_u=subplot(4,2,5);
     plot(time_trajall,state_input(:,2),'k-')
     hold(fig_u,'on')
     plot(time_trajall,state_sim(:,2),'-','color',usr_color.lightblue)
     plot(time_trajall,state_sim_ori(:,2),'-','color',usr_color.orange)
     axis(fig_u,[time_trajall(1) time_trajall(end) minyaxis maxyaxis]);
     set(gca,'Xgrid','on','Ygrid','on');
     xtickformat('%.0f')
     ytickformat('%.1f')
     set(gca,'XMinorTick','on','YMinorTick','on');
     xlabel('Time \rm[s\rm]' );
     ylabel('\itu \rm[m/s]')
     hold(fig_u, 'off')
     
      % vm velocity plot
      
     minyaxis =min(min(state_input(:,4)),min(state_sim_ori(:,4)))*1.1;
     maxyaxis =max(max(state_input(:,4)),max(state_sim_ori(:,4)))*1.1;
    
     fig_vm=subplot(4,2,6);
     plot(time_trajall,state_input(:,4),'k-')
     hold(fig_vm,'on')
     plot(time_trajall,state_sim(:,4),'-','color',usr_color.lightblue)
     plot(time_trajall,state_sim_ori(:,4),'-','color',usr_color.orange)
     axis(fig_vm,[time_trajall(1) time_trajall(end) minyaxis maxyaxis]);
     set(gca,'Xgrid','on','Ygrid','on');
     xtickformat('%.0f')
     ytickformat('%.1f')
     set(gca,'XMinorTick','on','YMinorTick','on');
     xlabel('Time \rm[s\rm]' );
     ylabel('\itv_m \rm[m/s]')
     hold(fig_vm, 'off')
     
     % r plot
     
     minyaxis = rad2deg(min(min(state_input(:,6)),min(state_sim(:,6))))*1.1;
     maxyaxis = rad2deg(max(max(state_input(:,6)),max(state_sim(:,6))))*1.1;
     
     fig_r=subplot(4,2,7);
     plot(time_trajall,rad2deg(state_input(:,6)),'k-')
     hold(fig_r,'on')
     plot(time_trajall,rad2deg(state_sim(:,6)),'-','color',usr_color.lightblue)
     plot(time_trajall,rad2deg(state_sim_ori(:,6)),'-','color',usr_color.orange)
     axis(fig_r,[time_trajall(1) time_trajall(end) minyaxis maxyaxis]);
     set(gca,'Xgrid','on','Ygrid','on');
     xtickformat('%.0f')
     ytickformat('%.1f')
     set(gca,'XMinorTick','on','YMinorTick','on');
     xlabel('Time \rm[s\rm]' );
     ylabel('\itr \rm[degree/s]')
     hold(fig_r, 'off')
     
    
     % control plot
    
     fig_control=subplot(4,2,8);
     plot(time_trajall,rad2deg(delta_input(:)),'k-')
     hold(fig_control,'on')
     plot(time_trajall,n_input,'-','color',usr_color.lightblue)
     axis(fig_control,[time_trajall(1) time_trajall(end) -36 36]);
     set(gca,'Xgrid','on','Ygrid','on');
     xtickformat('%.0f')
     ytickformat('%.1f')
     set(gca,'XMinorTick','on','YMinorTick','on');
     xlabel('Time \rm[s\rm]' );
     ylabel('\delta \rm[degree], \itn \rm[rps]')
     legend('\delta', '\itn','location','northeastoutside')
     hold(fig_control, 'off')
    
    % trajectory
    fig_traj=subplot(4,2,[1, 2, 3, 4]);
    plot(state_input(:, 1)/pp.lpp, state_input(:, 3)/pp.lpp,'k:','linewidth',1.5)
    hold on
    plot(state_sim(:, 1)/pp.lpp, state_sim(:, 3)/pp.lpp,'color',usr_color.lightblue,'linestyle',':','linewidth',1.5)
    plot(state_sim_ori(:, 1)/pp.lpp, state_sim_ori(:, 3)/pp.lpp,'color',usr_color.orange,'linestyle',':','linewidth',1.5)
    
        for i=1:step_max-start_time_step+1
            if i==1 || mod(i, 200)==0 || i==step_max
            % plot input ship
            Y = [state_input( i, 3 ) - pp.breadth / 2.0d0, state_input( i, 3 ) + pp.breadth / 2.0d0;...
                 state_input( i, 3 ) + pp.breadth / 2.0d0, state_input( i, 3 ) + pp.breadth / 2.0d0;...
                 state_input( i, 3 ) + pp.breadth / 2.0d0, state_input( i, 3 );...
                 state_input( i, 3 ), state_input( i, 3 ) - pp.breadth / 2.0d0;...
                 state_input( i, 3 ) - pp.breadth / 2.0d0, state_input( i, 3 ) - pp.breadth / 2.0d0 ]/pp.lpp;
            X = [state_input( i, 1 ) - pp.lpp / 2.0d0, state_input( i, 1 ) - pp.lpp / 2.0d0;...
                 state_input( i, 1 ) - pp.lpp / 2.0d0, state_input( i, 1 ) + 1.0d0;...
                 state_input( i, 1 ) + 1.0d0, state_input( i, 1 ) + pp.lpp / 2.0d0;...
                 state_input( i, 1 ) + pp.lpp / 2.0d0, state_input( i, 1 ) + 1.0d0;...
                 state_input( i, 1 ) + 1.0d0, state_input( i, 1 ) - pp.lpp / 2.0d0]/pp.lpp;
            Cy = state_input( i, 3 )/pp.lpp;
            Cx = state_input( i, 1 )/pp.lpp;
            rot =  state_input(i,5);
            %   Rotation of the coordinate
            x2 = cos(rot) * X - sin(rot) * Y + Cx - Cx * cos(rot) + Cy * sin(rot);
            y2 = sin(rot) * X + cos(rot) * Y + Cy - Cx * sin(rot) - Cy * cos(rot);           
%             plot( y2 , x2 , 'k-', 'Linewidth', 0.5 );
            plot( x2 , y2 , 'k-', 'Linewidth', 0.5 );
            
            % apparent wind of input data
%             windend = [state_input(i, 3)/pp.lpp, state_input(i, 1)/pp.lpp];

            windend = [state_input(i, 1)/pp.lpp, state_input(i, 3)/pp.lpp];
            windstart = [1 1];
            arrowratio = 0.2;
            windstart(1) = windend(1) + arrowratio* apparent_wind_velo(i) * sin(apparent_wind_dir(i) +state_input(i, 5));
            windstart(2) = windend(2) + arrowratio* apparent_wind_velo(i) * cos(apparent_wind_dir(i) +state_input(i, 5));
            % plot wind
            quiver(windstart(1), windstart(2), windend(1)-windstart(1), windend(2)-windstart(2), 0, 'color','r' ,'linewidth',2.0,...
                'MaxHeadSize',1.5)

            %plot sim ship
            Y = [state_sim( i, 3 ) - pp.breadth / 2.0d0, state_sim( i, 3 ) + pp.breadth / 2.0d0;...
                 state_sim( i, 3 ) + pp.breadth / 2.0d0, state_sim( i, 3 ) + pp.breadth / 2.0d0;...
                 state_sim( i, 3 ) + pp.breadth / 2.0d0, state_sim( i, 3 );...
                 state_sim( i, 3 ), state_sim( i, 3 ) - pp.breadth / 2.0d0;...
                 state_sim( i, 3 ) - pp.breadth / 2.0d0, state_sim( i, 3 ) - pp.breadth / 2.0d0 ]/pp.lpp;
            X = [state_sim( i, 1 ) - pp.lpp / 2.0d0, state_sim( i, 1 ) - pp.lpp / 2.0d0;...
                 state_sim( i, 1 ) - pp.lpp / 2.0d0, state_sim( i, 1 ) + 1.0d0;...
                 state_sim( i, 1 ) + 1.0d0, state_sim( i, 1 ) + pp.lpp / 2.0d0;...
                 state_sim( i, 1 ) + pp.lpp / 2.0d0, state_sim( i, 1 ) + 1.0d0;...
                 state_sim( i, 1 ) + 1.0d0, state_sim( i, 1 ) - pp.lpp / 2.0d0]/pp.lpp;
            Cy = state_sim( i, 3 )/pp.lpp;
            Cx = state_sim( i, 1 )/pp.lpp;
            rot =  state_sim(i,5);
            %   Rotation of the coordinate
            x2 = cos(rot) * X - sin(rot) * Y + Cx - Cx * cos(rot) + Cy * sin(rot);
            y2 = sin(rot) * X + cos(rot) * Y + Cy - Cx * sin(rot) - Cy * cos(rot); 
            
%             plot( y2 , x2 ,'color',usr_color.lightblue,'linestyle','-','linewidth',0.5)
             plot( x2 , y2 ,'color',usr_color.lightblue,'linestyle','-','linewidth',0.5)
            %plot EFD Coef. ship
             Y = [state_sim_ori( i, 3 ) - pp.breadth / 2.0d0, state_sim_ori( i, 3 ) + pp.breadth / 2.0d0;...
                 state_sim_ori( i, 3 ) + pp.breadth / 2.0d0, state_sim_ori( i, 3 ) + pp.breadth / 2.0d0;...
                 state_sim_ori( i, 3 ) + pp.breadth / 2.0d0, state_sim_ori( i, 3 );...
                 state_sim_ori( i, 3 ), state_sim_ori( i, 3 ) - pp.breadth / 2.0d0;...
                 state_sim_ori( i, 3 ) - pp.breadth / 2.0d0, state_sim_ori( i, 3 ) - pp.breadth / 2.0d0 ]/pp.lpp;
            X = [state_sim_ori( i, 1 ) - pp.lpp / 2.0d0, state_sim_ori( i, 1 ) - pp.lpp / 2.0d0;...
                 state_sim_ori( i, 1 ) - pp.lpp / 2.0d0, state_sim_ori( i, 1 ) + 1.0d0;...
                 state_sim_ori( i, 1 ) + 1.0d0, state_sim_ori( i, 1 ) + pp.lpp / 2.0d0;...
                 state_sim_ori( i, 1 ) + pp.lpp / 2.0d0, state_sim_ori( i, 1 ) + 1.0d0;...
                 state_sim_ori( i, 1 ) + 1.0d0, state_sim_ori( i, 1 ) - pp.lpp / 2.0d0]/pp.lpp;
            Cy = state_sim_ori( i, 3 )/pp.lpp;
            Cx = state_sim_ori( i, 1 )/pp.lpp;
            rot =  state_sim_ori(i,5);
            %   Rotation of the coordinate
            x2 = cos(rot) * X - sin(rot) * Y + Cx - Cx * cos(rot) + Cy * sin(rot);
            y2 = sin(rot) * X + cos(rot) * Y + Cy - Cx * sin(rot) - Cy * cos(rot);           
%             plot( y2 , x2 ,'color',usr_color.orange,'linestyle','-','linewidth',0.5)
            plot( x2 , y2 ,'color',usr_color.orange,'linestyle','-','linewidth',0.5)
            
            end
        end
    axis equal
    axis ij
%     xlabel ('$$\hat{y}/L_{pp}$$ ','interpreter','latex')
%     ylabel ('$$\hat{x}/L_{pp}$$ ','interpreter','latex')
    minxaxis =min(min(state_input(:,1)),min(state_sim(:,1)))/pp.lpp-1;
    minyaxis =min(min(state_input(:,3)),min(state_sim(:,3)))/pp.lpp-1;
    maxxaxis =max(max(state_input(:,1)),max(state_sim(:,1)))/pp.lpp+1;
    maxyaxis =max(max(state_input(:,3)),max(state_sim(:,3)))/pp.lpp+1;
    axis([minxaxis, maxxaxis, minyaxis, maxyaxis])
    ylabel ('$$\hat{y}/L_{pp}$$ ','interpreter','latex')
    xlabel ('$$\hat{x}/L_{pp}$$ ','interpreter','latex')
    legend('Input','CMA-ES Coef.', 'EFD Coef.','location','best')
    hold off
    %%
    % save figure as image
    fig = gcf;

    figurename = strcat(save_directory_path,'/',condition_name,'_', inputfilename,'_DTeval=',num2str(CmaSettings.validation_period),'_all_traj');
    SaveFig(fig,figurename)
end