function [CmaSettings, itr, minj, minarffeas, arffeas,arf,df,dx,sigma,vecD,xmean,...
            sigma_sqrt_eigC,sigma_sqrt_daigC,dim,itrplot,Xopt] ...
            = read_outputcmaes(fid_CMAsetting, fid_CMAprob, fid_outputcmaes,CmaSettings)
    for i=1:5
        tline = fgetl(fid_CMAsetting);

        tline = extractBefore( tline, "!" );
        Data = str2double( tline );
        if i==2
            CmaSettings.switch_wind = Data;
        elseif i==3
            CmaSettings.type_si = Data;
        elseif i==4
            CmaSettings.number_files = Data;
        elseif i==5
            CmaSettings.integration_period = Data;
        end
    end
    
    initprob = textscan(fid_CMAprob,'%d %c %c %c');
    dim  = cell2mat(initprob(1)); % Number of variables
    % Data input and sort

    j=0;


    tline = fgetl(fid_outputcmaes);
    while ischar(tline)
%         disp(tline) % comment out if you dont need to show each itr
        tline = fgetl(fid_outputcmaes);
        if ischar(tline)==0
            break;
        end
        Data = str2num(tline);
        % Find Infinity and NaN output. 
        % In this case, utilize the previous values
        TF1 = contains(tline,"Infinity");
        TF2 = contains(tline,"NaN");
        %
        j=j+1;
        %
        if j==1
            [rows,cols] = size(Data);
            vecD=zeros(rows,dim);
            xmean=zeros(rows,dim);
            diagSqrtC=zeros(rows,dim);
        end
        %
        if TF1==1 || TF2==1
            itr    (j)=itr    (j-1);
            arffeas(j)=arffeas(j-1);
            arf    (j)=arf    (j-1);
            df     (j)=df     (j-1);
            dx     (j)=dx     (j-1);
            sigma  (j)=sigma  (j-1);
            %
            for i=1:dim
                vecD     (j,i)=vecD     (j-1,i);
                xmean    (j,i)=xmean    (j-1,i);
                diagSqrtC(j,i)=diagSqrtC(j-1,i);
            end
        else
            itr    (j)=Data(1,1);
            arffeas(j)=Data(1,2);
            arf    (j)=Data(1,3);
            df     (j)=Data(1,4);
            dx     (j)=Data(1,5);
            sigma  (j)=Data(1,6);
            %
            for i=1:dim
                vecD     (j,i)=Data(1,6      +i);
                xmean    (j,i)=Data(1,6+  dim+i);
                diagSqrtC(j,i)=Data(1,6+2*dim+i);
            end
        end
     end
    %
    % choose minimum J

    [ minarffeas, minj]=min(arffeas);
    [rows,cols] = size(itr);

    sigma_sqrt_eigC = zeros(cols, dim);
    sigma_sqrt_daigC= zeros(cols, dim);

    for i=1:cols
        sigma_sqrt_eigC (i,:)=sigma(i)*vecD(i,:);
        sigma_sqrt_daigC(i,:)=sigma(i)*diagSqrtC(i,:);
    end
    %
    %   Set itrplot (horizontal axis of the graph)
    Npower=fix(log10(itr(cols)));
    itrplot=round(itr(cols),-Npower);
    if itrplot<itr(cols)
        itrplot=itrplot+10^(Npower);
    end
    %
    fclose(fid_outputcmaes);
    Xopt=zeros(1,dim);
    for i=1:dim
        Xopt(i)=xmean(minj,i);
    end
end