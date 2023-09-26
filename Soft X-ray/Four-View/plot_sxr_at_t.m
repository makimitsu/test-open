function plot_sxr_at_t(grid2D,data2D,date,shot,t,show_xpoint,show_localmax,show_flux_surface,start,interval,save,SXRfilename,filter,NL)

newProjectionNumber = 80;
newGridNumber = 100;

if filter & NL
    options = 'NLF_NLR';
elseif ~filter & NL
    options = 'LF_NLR';
elseif filter & ~NL
    options = 'NLF_LR';
else
    options = 'LF_LR';
end
matrixFolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/' ...
    ,options,'/',num2str(date),'/shot',num2str(shot));

if exist(matrixFolder,'dir') == 0
    doCalculation = true;
else
    doCalculation = false;
end

% 再構成計算に必要なパラメータを計算するなら読み込む、しない場合も範囲に関しては読み込む
parameterFile = '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/Four-View_Simulation/parameters.mat';

if doCalculation
    if isfile(parameterFile)
        load(parameterFile, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
            's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
        if newProjectionNumber ~= N_projection || newGridNumber ~= N_grid
            disp('Different parameters - Start calculation!');
            get_parameters(newProjectionNumber,newGridNumber,parameterFile);
            load(parameterFile, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
                's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
        end
    else
        disp('No parameter - Start calculation!');
        get_parameters(newProjectionNumber,newGridNumber,parameterFile);
        load(parameterFile, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
            's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');   
    end
else
    load(parameterFile,'range');
end
    
number = (t-start)/interval+1;
doPlot = false;

if doCalculation
    [Iwgn1,Iwgn2,Iwgn3,Iwgn4] = get_sxr_image(date,number,newProjectionNumber,SXRfilename,filter);
    
    EE1 = get_distribution(M,K,gm2d1,U1,s1,v1,Iwgn1,doPlot,NL);
    EE2 = get_distribution(M,K,gm2d2,U2,s2,v2,Iwgn2,doPlot,NL);
    EE3 = get_distribution(M,K,gm2d3,U3,s3,v3,Iwgn3,doPlot,NL);
    EE4 = get_distribution(M,K,gm2d4,U4,s4,v4,Iwgn4,doPlot,NL);
else
    matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
    load(matrixPath,'EE1','EE2','EE3','EE4');
end

f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.2,0.8,0.8];

EE = cat(3,EE1,EE2,EE3,EE4);
plot_save_sxr(grid2D,data2D,range,date,shot,t,EE,show_localmax,show_flux_surface,show_xpoint,save,filter,NL);

end