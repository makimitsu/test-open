function plot_sxr_at_t(Date,FiberImageNumber,SXRImageNumber)

TIFImagePath = strcat("G:\My Drive\X-ray\Data\TIF\" + num2str(Date));
% FiberImageNumber = 40;
% SXRImageNumber = 40;

% 目的の時間をタイミング番号に変換します
TimeSeriesIndex = (t-start)/interval+1;
doPlot = false;

if filter && NL
    options = 'NLF_NLR'
elseif ~filter && NL
    options = 'LF_NLR'
elseif filter && ~NL
    options = 'NLF_LR'
else
    options = 'LF_LR'
end

ResultFolder = strcat('G:\My Drive\X-ray\Data\SXROUT\Matrix\' ...
    ,options,'/',num2str(Date),'/shot',num2str(shot));

if exist(ResultFolder,'dir') == 0
    GenerateResult = true;
else
    GenerateResult = false;
end

% 再構成計算をする場合は，パラメータを読み込み、しない場合は，表示範囲に関してのみ読み込む
MatrixFilePath = '/Users/yuleo/Documents/GitHub/test-open/Soft X-ray/Four-View_Simulation/parameters.mat';

if GenerateResult
    if isfile(MatrixFilePath)
        load(MatrixFilePath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
            's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
        if FiberImageNumber ~= N_projection || SXRImageNumber ~= N_grid
            disp('Different parameters - Generating matrixes!');
            get_parameters(FiberImageNumber,SXRImageNumber,MatrixFilePath);
            load(MatrixFilePath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
                's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
        end
    else
        disp('Not found parameter - Generating matrixes!');
        get_parameters(FiberImageNumber,SXRImageNumber,MatrixFilePath);
        load(MatrixFilePath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
            's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');   
    end
else
    load(MatrixFilePath,'range');
end

if GenerateResult
    % ベクトルイメージをTIF画像から取得します
    [Iwgn1,Iwgn2,Iwgn3,Iwgn4] = get_sxr_image(Date,TimeSeriesIndex,FiberImageNumber,TIFImagePath,ApplyFilter);
    % 再構成像を計算します
    EE1 = get_distribution(M,K,gm2d1,U1,s1,v1,Iwgn1,doPlot,NL);
    EE2 = get_distribution(M,K,gm2d2,U2,s2,v2,Iwgn2,doPlot,NL);
    EE3 = get_distribution(M,K,gm2d3,U3,s3,v3,Iwgn3,doPlot,NL);
    EE4 = get_distribution(M,K,gm2d4,U4,s4,v4,Iwgn4,doPlot,NL);
else
    matrixPath = strcat(ResultFolder,'/',num2str(TimeSeriesIndex),'.mat');
    load(matrixPath,'EE1','EE2','EE3','EE4');
end

f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.2,0.8,0.4];

plot_save_sxr(range,EE1,EE2,EE3,EE4);

end