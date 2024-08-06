clear
close all
clc

scripts = dir('**/postprocess/*Result*.m'); % use wild card (*) to squeeze the searching

for numberofscript=1:length(scripts)
    
    if numberofscript == 1
        % because of the "clear" command of each scpript, whole workspace
        % are cleared eachtime, so path of scprit need to read every time.
        save scripts.mat
    end    
    script_fullpath = strcat(scripts(numberofscript).folder, '/',scripts(numberofscript).name)
    run(script_fullpath)
    aaa =1
    load scripts.mat
end

