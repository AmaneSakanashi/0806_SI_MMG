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

%% parameters
num_randomseed = 5;
num_trainobj   = 3;
num_traindata  = 3;
num_testobj    = 3;
num_testdata   = 5;
test_berth_num = 2;

num_totalcase =  num_randomseed*num_trainobj*num_traindata;
num_case      =  num_totalcase/num_randomseed;

csv_suffix = '*100_Objfunc.csv';
%% allocate struct 
EFDdata.testJ1 = zeros(1,num_testdata+2);
EFDdata.testJ2 = zeros(1,num_testdata+2);
EFDdata.testJ3 = zeros(1,num_testdata+2);


casedata(1).randomseed = 0;
casedata(1).trainobj   = 0;
casedata(1).traindata  = "";
casedata(1).testJ1     = zeros(1,num_testdata+2);
casedata(1).testJ2     = zeros(1,num_testdata+2);
casedata(1).testJ3     = zeros(1,num_testdata+2);
casedata(1).casename   = "";
casedata(1).use        = false;

stat(1).trainobj   = 0;
stat(1).traindata  = "R";
stat(1).testJ1     = zeros(num_randomseed,num_testdata+2);
stat(1).testJ2     = zeros(num_randomseed,num_testdata+2);
stat(1).testJ3     = zeros(num_randomseed,num_testdata+2);

%%
% search outfiles directiory
outputfiles = dir('../case*');

%%

% open output files
for cases = 1:length(outputfiles)
% for cases=2 % use this for sentense for debug
    % set attributes of each cases
    casesetting = extractAfter(outputfiles(cases).name,'100-');
    
    casedata(cases).randomseed = str2num(casesetting(1));
    casedata(cases).trainobj = str2num(casesetting(3));
    casedata(cases).traindata = casesetting(5:end);
    casedata(cases).casename = outputfiles(cases).name;
    
    directory_path = strcat(outputfiles(cases).folder,'/', outputfiles(cases).name, '/');
    save_directory_path = strcat(directory_path,'/','processed','/');
    temp_data = read_objfunc(save_directory_path,csv_suffix);
    for testno =1:num_testdata
            EFDdata.testJ1(testno) = temp_data(testno).obj(1,2).Variables;
            EFDdata.testJ2(testno) = temp_data(testno).obj(2,2).Variables;
            EFDdata.testJ3(testno) = temp_data(testno).obj(3,2).Variables;
            
            casedata(cases).testJ1(testno) = temp_data(testno).obj(1,1).Variables;
            casedata(cases).testJ2(testno) = temp_data(testno).obj(2,1).Variables;
            casedata(cases).testJ3(testno) = temp_data(testno).obj(3,1).Variables;
    end
    % total value of J
    EFDdata.testJ1(testno+1) = sum(EFDdata.testJ1(1:testno));
    EFDdata.testJ2(testno+1) = sum(EFDdata.testJ2(1:testno));
    EFDdata.testJ3(testno+1) = sum(EFDdata.testJ3(1:testno));
    casedata(cases).testJ1(testno+1) = sum(casedata(cases).testJ1(1:testno));
    casedata(cases).testJ2(testno+1) = sum(casedata(cases).testJ2(1:testno));
    casedata(cases).testJ3(testno+1) = sum(casedata(cases).testJ3(1:testno));
    
    % sub total of berthing maneuver
    EFDdata.testJ1(testno+2) = sum(EFDdata.testJ1(1:test_berth_num));
    EFDdata.testJ2(testno+2) = sum(EFDdata.testJ2(1:test_berth_num));
    EFDdata.testJ3(testno+2) = sum(EFDdata.testJ3(1:test_berth_num));
    casedata(cases).testJ1(testno+2) = sum(casedata(cases).testJ1(1:test_berth_num));
    casedata(cases).testJ2(testno+2) = sum(casedata(cases).testJ2(1:test_berth_num));
    casedata(cases).testJ3(testno+2) = sum(casedata(cases).testJ3(1:test_berth_num));
    
end


%% save case data to csv
vartype = {'string','double','double','double','double','double','double','double'};
varName = {'case_name','','','','','','total','subtotal_berthing'};
for testno=1:num_testdata
    varName{testno+1} = cell2mat(extractBetween(temp_data(testno).name,'_','.xls'));
end
varName = strrep(varName,'-','_');
savetable = table('size',[length(outputfiles)+1,num_testdata+3],'VariableTypes',vartype,...
    'VariableNames',varName);
savetable{1,1} = "EFD";
savetable{1,2:num_testdata+3} = EFDdata.testJ2;

for cases = 2:length(outputfiles)+1
    savetable{cases,1} =convertCharsToStrings(casedata(cases-1).casename);
    savetable{cases,2:num_testdata+3} = casedata(cases-1).testJ2;
end
writetable(savetable,"../statistics_eachcaseJ2.csv")


%% for J1 and J3
% J1
savetable{1,2:num_testdata+3} = EFDdata.testJ1;

for cases = 2:length(outputfiles)+1
%     savetable{cases,1} =convertCharsToStrings(casedata(cases-1).casename);
    savetable{cases,2:num_testdata+3} = casedata(cases-1).testJ1;
end

writetable(savetable,"../statistics_eachcaseJ1.csv")

% J3

savetable{1,2:num_testdata+3} = EFDdata.testJ3;

for cases = 2:length(outputfiles)+1
%     savetable{cases,1} =convertCharsToStrings(casedata(cases-1).casename);
    savetable{cases,2:num_testdata+3} = casedata(cases-1).testJ3;
end

writetable(savetable,"../statistics_eachcaseJ3.csv")

%% statistical analysis
for subcase = 1:num_case
    stat(subcase).trainobj   = ceil(subcase/num_trainobj);
    if mod(subcase, num_traindata)==1
        stat(subcase).traindata  = "R";
    elseif mod(subcase, num_traindata)==2
        stat(subcase).traindata  = "TR";
    else
        stat(subcase).traindata  = "TZR";
    end
    stat(subcase).testJ1     = zeros(num_randomseed,num_testdata+2);
    stat(subcase).testJ2     = zeros(num_randomseed,num_testdata+2);
    stat(subcase).testJ3     = zeros(num_randomseed,num_testdata+2);
    stat(subcase).testJ2_mean = zeros(1,num_testdata+2);
    stat(subcase).testJ2_std = zeros(1,num_testdata+2);
end
for cases = 1:length(outputfiles)
    for subcase = 1:num_case
        if (casedata(cases).trainobj == stat(subcase).trainobj &&...
               casedata(cases).traindata == stat(subcase).traindata)
            stat(subcase).testJ1(casedata(cases).randomseed+1,:) = casedata(cases).testJ1;
            stat(subcase).testJ2(casedata(cases).randomseed+1,:) = casedata(cases).testJ2;
            stat(subcase).testJ3(casedata(cases).randomseed+1,:) = casedata(cases).testJ3;
        end
    end
end

for subcase = 1:num_case
    non_zero_index = find(stat(subcase).testJ2(1:num_randomseed,1));
    stat(subcase).testJ2_mean = mean(stat(subcase).testJ2(non_zero_index,:));
    stat(subcase).testJ2_std  = std(stat(subcase).testJ2(non_zero_index,:));
    
    stat(subcase).testJ1_mean = mean(stat(subcase).testJ1(non_zero_index,:));
    stat(subcase).testJ1_std  = std(stat(subcase).testJ1(non_zero_index,:));
    
    stat(subcase).testJ3_mean = mean(stat(subcase).testJ3(non_zero_index,:));
    stat(subcase).testJ3_std  = std(stat(subcase).testJ3(non_zero_index,:));
end

%% save statistics to csv
statVarName = {'train_obj','traindata','mean_berth1','mean_berth2','mean_random', 'mean_turn','mean_z','mean_total','mean_subtotal_berthing',...
                'std_berth1','std_berth2','std_random', 'std_turn','std_z','std_total','std_subtotal_berthing'};
statvartype = {'double','string','double','double','double','double','double','double','double','double','double','double','double','double','double','double'};

statsavetable = table('size',[num_case,2+2*(num_testdata+2)],'VariableTypes',statvartype,...
    'VariableNames',statVarName);
for subcase = 1:num_case
    statsavetable{subcase,1} = stat(subcase).trainobj;
    statsavetable{subcase,2} = stat(subcase).traindata;
    statsavetable{subcase,3:num_testdata+4} =stat(subcase).testJ2_mean;
    statsavetable{subcase,num_testdata+5:end} =stat(subcase).testJ2_std;
end
writetable(statsavetable,"../statistics_meanstd.csv")
%%
tabVarName = {'train_obj','traindata', 'bearh1', 'bearth2', 'random', 'turn', 'z', 'total','subtotal_berth'};
tabVerType = {'double','string', 'string','string','string','string','string','string','string'};
tab_statsavetable = table('size',[num_case+1,2+num_testdata+2],'VariableType',tabVerType,...
                            'VariableName',tabVarName);
% EFD                        
tab_statsavetable{1,1} = 0;
tab_statsavetable{1,2} = "EFD";
tab_statsavetable{1,3:end}= strings(1,num_testdata+2);
tab_statsavetable{1,3:end}=transpose(string(num2str(EFDdata.testJ2(:),'%.1f')));

% tale of string style for journal paper
for subcase = 1:num_case
    tab_statsavetable{subcase+1,1} = stat(subcase).trainobj;
    tab_statsavetable{subcase+1,2} = stat(subcase).traindata;
    mean_str = strings(1,num_testdata+2);
    std_str  = strings(1,num_testdata+2);
    mean_str(:) = num2str(stat(subcase).testJ2_mean(:),'%.1f');
    std_str(:)  = num2str(stat(subcase).testJ2_std(:),'%.1f');
    tab_statsavetable{subcase+1,3:num_testdata+4} =transpose(strcat(mean_str(:),' (',std_str(:),')'))

end

writetable(tab_statsavetable,"../statistics_tableJ2.csv")
%% J1 

% EFD                        
tab_statsavetable{1,1} = 0;
tab_statsavetable{1,2} = "EFD";
tab_statsavetable{1,3:end}= strings(1,num_testdata+2);
tab_statsavetable{1,3:end}=transpose(string(num2str(EFDdata.testJ1(:),'%.1f')));

% tale of string style for journal paper
for subcase = 1:num_case
    tab_statsavetable{subcase+1,1} = stat(subcase).trainobj;
    tab_statsavetable{subcase+1,2} = stat(subcase).traindata;
    mean_str = strings(1,num_testdata+2);
    std_str  = strings(1,num_testdata+2);
    mean_str(:) = num2str(stat(subcase).testJ1_mean(:),'%.1f');
    std_str(:)  = num2str(stat(subcase).testJ1_std(:),'%.1f');
    tab_statsavetable{subcase+1,3:num_testdata+4} =transpose(strcat(mean_str(:),' (',std_str(:),')'))

end

writetable(tab_statsavetable,"../statistics_tableJ1.csv")

%% J3

% EFD                        
tab_statsavetable{1,1} = 0;
tab_statsavetable{1,2} = "EFD";
tab_statsavetable{1,3:end}= strings(1,num_testdata+2);
tab_statsavetable{1,3:end}=transpose(string(num2str(EFDdata.testJ3(:),'%.1f')));

% tale of string style for journal paper
for subcase = 1:num_case
    tab_statsavetable{subcase+1,1} = stat(subcase).trainobj;
    tab_statsavetable{subcase+1,2} = stat(subcase).traindata;
    mean_str = strings(1,num_testdata+2);
    std_str  = strings(1,num_testdata+2);
    mean_str(:) = num2str(stat(subcase).testJ3_mean(:),'%.1f');
    std_str(:)  = num2str(stat(subcase).testJ3_std(:),'%.1f');
    tab_statsavetable{subcase+1,3:num_testdata+4} =transpose(strcat(mean_str(:),' (',std_str(:),')'))

end

writetable(tab_statsavetable,"../statistics_tableJ3.csv")
%%
function[objfuncfiles]=read_objfunc(save_directory_path,csv_suffix)
    objfuncfiles=dir(strcat(save_directory_path,csv_suffix));
    
    for files=1:length(objfuncfiles)
        filepath = strcat(objfuncfiles(files).folder,'/',objfuncfiles(files).name);
        data = readtable(filepath);
        objfuncfiles(files).obj = data;
    end
end