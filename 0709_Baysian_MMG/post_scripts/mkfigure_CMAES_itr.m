function mkfigure_CMAES_itr(itr,arffeas,sigma_sqrt_eigC,sigma_sqrt_daigC,dim,directory_path)
    %   Prepare the figure block

    S = get(0,'ScreenSize');
    figure('Color','white','Position',[2*S(4)/12 2*S(4)/12 5*S(3)/12 9*S(4)/12])
   
    %%%   Make drawing
    % Graph about arffeas and arf
    subplot(2,2,1);

%     semilogy(itr(:),arffeas(:),'k');
    plot(itr(:),arffeas(:),'k');
    grid on;
    plotsize1 = [0 max(itr) 0 2.0e4];
    axis([0 max(itr) 0 2.0e4]);
    set(gca,'Xgrid','on','Ygrid','on');
    xtickformat('%.1f');
    ytickformat('%.1f');
    set(gca,'XMinorTick','on','YMinorTick','on');
    xlabel( '\sl Iteration'); 
    ylabel('\it Obj. func.')
    pbaspect([1 1 1]);
    legend('\sl J','Location','best','Orientation','vertical');
    
    % Graph about arffeas-min and arf-min
    subplot(2,2,2);

    semilogy(itr(:),arffeas(:)-min(arffeas(:)),'k');
    set(gca,'Xgrid','on','Ygrid','on');
    xtickformat('%.1f');
    ytickformat('%.0f');
    xlim([0 max(itr)])
    set(gca,'XMinorTick','on','YMinorTick','on');
    xlabel( '\sl Iteration' ); 
      ylabel('\it Obj. func.')
    %ylabel('\slVelocity \rm[\slm / s\rm]', 'FontName','Times','FontSize',16 );
    grid on;
    pbaspect([1 1 1]);
    legend('\sl J\rm-\slmin\rm[\slJ\rm]','Location','best','Orientation','vertical');
 
    % Graph about sigma*sqrt(eig(C)) (Sqrt. of Eigenvalues of Cov.)
    subplot(2,2,3);
    for i=1:dim
       semilogy(itr(:),sigma_sqrt_eigC(:,i),'k');hold on; 
    end
    plotsize1 = [0 max(itr) 1.0e-10 1.0e+0];
    axis(plotsize1);
    set(gca,'Xgrid','on','Ygrid','on');
    xtickformat('%.1f');
    ytickformat('%.0f');
    set(gca,'XMinorTick','on','YMinorTick','on');
    xlabel( '\sl Iteration'); 
    ylabel('\sl Sqrt. of Eigenvalues of Cov.' );
    grid on;
    pbaspect([1 1 1]);

    
    % Graph about sigma*sqrt(diag(C)) (Sqrt. of Diagonal Elem. of Cov.)
    subplot(2,2,4);
    for i=1:dim
       semilogy(itr(:),sigma_sqrt_daigC(:,i),'k');hold on; 
    end
    plotsize1 = [0 max(itr) 1.0e-15 1.0e+0];
    axis(plotsize1);
    set(gca,'Xgrid','on','Ygrid','on');
    xtickformat('%.1f');
    ytickformat('%.0f');
    set(gca,'XMinorTick','on','YMinorTick','on');
    xlabel( '\sl Iteration', 'FontName','Times','FontSize',12 ); 

    ylabel('\slSqrt. of Diagonal Elem. of Cov.' );
    grid on;
    pbaspect([1 1 1]);

    savename = strcat(directory_path,'conv') ;
    SaveFig(gcf,savename)
    
    % J history only
    figJhist=figure(2);
    figJhist.Color='white';
    semilogy(itr(:),arffeas(:)-min(arffeas(:)),'k');
    set(gca,'Xgrid','on','Ygrid','on');
    xtickformat('%.1f');
    ytickformat('%.0f');
    xlim([0 max(itr)])
    set(gca,'XMinorTick','on','YMinorTick','on');
    xlabel( 'Iteration' ); 
    ylabel('Obj. func.')
    %ylabel('\slVelocity \rm[\slm / s\rm]', 'FontName','Times','FontSize',16 );
    grid on;
    pbaspect([1 1 1]);
    legend('\sl J\rm-\slmin\rm[\slJ\rm]','Location','SouthWest','Orientation','vertical');

    
    savename = strcat(directory_path,'J_minus_minJ') ;
    SaveFig(gcf,savename)
    
end