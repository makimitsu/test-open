function [] = plot_sxr_multi(grid2D,data2D,date,shot,show_xpoint,show_localmax,start,interval,doSave,SXRfilename,doFilter,NL)
% plot SXR emission on psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   integer: date, date of experiment
%   integer: shot, number of shot
%   boolean: show_xpoint, option for showing the x-point
%   boolean: show_localmax, option for showing the local maximum point
%   boolean: show_flux_surface, option for showing the flux_surface
%   integer: start, start time (us)
%   integer: interval, interval time of the framing camera (us)
%   boolean: save, option for saving the reconstruction result
%   string: SXRfilename, name of the SXR image file
%   boolean: filter, option for applying non-linear mean (NLM) filter
%   boolean: NL, option for using non-linear reconstruction

if doFilter & NL
    options = 'NLF_NLR';
elseif ~doFilter & NL
    options = 'LF_NLR';
elseif doFilter & ~NL
    options = 'NLF_LR';
else
    options = 'LF_LR';
end

dirPath = getenv('SXR_MATRIX_DIR');
matrixFolder = strcat(dirPath,'/',options,'/',num2str(date),'/shot',num2str(shot));
if exist(matrixFolder,'dir') == 0
    doCalculation = true;
    mkdir(matrixFolder);
else
    doCalculation = false; 
end

% 再構成計算に必要なパラメータを計算するなら読み込む
parameterFile = 'parameters.mat';
if doCalculation
    disp('No matrix data -- Start calculation');
    newProjectionNumber = 50;
    newGridNumber = 90;
    if isfile(parameterFile)
        disp(strcat('Loading parameter from ',which(parameterFile)))
        load(parameterFile,'gm2d1','gm2d2','gm2d3','gm2d4', ...
                'U1','U2','U3','U4','s1','s2','s3','s4', ...
                'v1','v2','v3','v4','M','K','range','N_projection','N_grid');
            
        if newProjectionNumber ~= N_projection || newGridNumber ~= N_grid
            disp('Different parameters - Start calculation!');
            get_parameters(newProjectionNumber,newGridNumber,parameterFile);
            load(parameterFile,'gm2d1','gm2d2','gm2d3','gm2d4', ...
                    'U1','U2','U3','U4','s1','s2','s3','s4', ...
                    'v1','v2','v3','v4','M','K','range');
        end
    else
        disp('No parameters - Start calculation!');
        get_parameters(newProjectionNumber,newGridNumber,parameterFile);
        load(parameterFile,'gm2d1','gm2d2','gm2d3','gm2d4', ...
                'U1','U2','U3','U4','s1','s2','s3','s4', ...
                'v1','v2','v3','v4','M','K','range');
    end

    % 生画像の取得
    rawImage = imread(SXRfilename);
    
    % 非線形フィルターをかける（必要があれば）
    if doFilter
        if date == 230920 && shot == 15
            load("230920\230920_shot15_denoisedTIF_w91_c5_d0625.mat","TIFImage");
            rawImage = TIFImage;
        else
            patch = rawImage(1:200,1:200);
            patchSq = patch.^2;
            edist = sqrt(sum(patchSq,3));
            patchSigma = sqrt(var(edist(:)));
            [rawImage,~] = imnlmfilt(rawImage,'SearchWindowSize',91,'ComparisonWindowSize',5,'DegreeOfSmoothing',patchSigma*0.625);
        end
    end
else
    disp(strcat('Loading matrix from :',matrixFolder))
    load(parameterFile,'range');
end


times = start:interval:(start+interval*7);
doPlot = false;

% f = figure;
% f.Units = 'normalized';
% f.Position = [0.1,0.2,0.8,0.4];


if doSave
    f = figure;
    f.Units = 'normalized';
    f.Position = [0.1,0.2,0.8,0.8];
end

for t = times
    number = (t-start)/interval+1;
    % disp(clc_flag);
    
    if doCalculation
%         ベクトル形式の画像データの読み込み
        [VectorImage1,VectorImage2, VectorImage3, VectorImage4] = get_sxr_image(date,number,newProjectionNumber,rawImage,doFilter);

%         再構成計算
        EE1 = get_distribution(M,K,gm2d1,U1,s1,v1,VectorImage1,doPlot,NL);
        EE2 = get_distribution(M,K,gm2d2,U2,s2,v2,VectorImage2,doPlot,NL);
        EE3 = get_distribution(M,K,gm2d3,U3,s3,v3,VectorImage3,doPlot,NL);
        EE4 = get_distribution(M,K,gm2d4,U4,s4,v4,VectorImage4,doPlot,NL);
        
%         再構成結果を保存するファイルを作成、保存
        
        matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
        save(matrixPath,'EE1','EE2','EE3','EE4');

        % savepath_one = strcat(savefolder,'/',num2str(number),'_one.txt');
        % savepath_two = strcat(savefolder,'/',num2str(number),'_two.txt');
        % savepath_three = strcat(savefolder,'/',num2str(number),'_three.txt');
        % savepath_four = strcat(savefolder,'/',num2str(number),'_four.txt');
        % writematrix(EE1,savepath_one);
        % writematrix(EE2,savepath_two);
        % writematrix(EE3,savepath_three);
        % writematrix(EE4,savepath_four);
        
    else
        matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
        % disp(strcat('Loading result matrix from ',which(matrixPath)));
        load(matrixPath,'EE1','EE2','EE3','EE4');

        % loadpath_one = strcat(savefolder,'/',num2str(number),'_one.txt');
        % loadpath_two = strcat(savefolder,'/',num2str(number),'_two.txt');
        % loadpath_three = strcat(savefolder,'/',num2str(number),'_three.txt');
        % loadpath_four = strcat(savefolder,'/',num2str(number),'_four.txt');
        % EE1 = readmatrix(loadpath_one);
        % EE2 = readmatrix(loadpath_two);
        % EE3 = readmatrix(loadpath_three);
        % EE4 = readmatrix(loadpath_four);

    end
    
    EE = cat(3,EE1,EE2,EE3,EE4);

    if ~doSave
        f = figure;
        f.Units = 'normalized';
        f.Position = [0.1,0.2,0.8,0.8];
    end

    plot_save_sxr(grid2D,data2D,range,date,shot,t,EE,show_localmax,show_xpoint,doSave,doFilter,NL);

end

if doSave
    close(f);
end

end