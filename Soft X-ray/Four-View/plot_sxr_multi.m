function [] = plot_sxr_multi(PCBdata,SXR)
% grid2D = PCBdata.grid2D;
% data2D = PCBdata.data2D;
date = SXR.date;
shot = SXR.shot;
% show_xpoint = SXR.show_xpoint;
% show_localmax = SXR.show_localmax;
start = SXR.start;
interval = SXR.interval;
doSave = SXR.doSave;
doFilter = SXR.doFilter;
ReconMethod = SXR.ReconMethod;
SXRfilename = SXR.SXRfilename;

addpath '/Users/shohgookazaki/Documents/GitHub/test-open/Soft X-ray/Machine_Learning/code'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス


if doFilter == 0
    if ReconMethod == 0
        options = 'NLF_TP';
    elseif ReconMethod == 1
        options = 'NLF_MFI';
    elseif ReconMethod == 2
        options = 'NLF_MEM';
    elseif ReconMethod == 3
        options  = 'NLF_cGAN';
    end
else
    if ReconMethod == 0
        options = 'LF_TP';
    elseif ReconMethod == 1
        options = 'LF_MFI';
    elseif ReconMethod == 2
        options = 'LF_MEM';
    elseif ReconMethod == 3
        options  = 'LF_cGAN';
    end
end

dirPath = getenv('SXR_MATRIX_DIR');
matrixFolder = strcat(dirPath,'/',options,'/',num2str(date),'/shot',num2str(shot));
if exist(matrixFolder,'dir') == 0
    doCalculation = true;
    mkdir(matrixFolder);
elseif length(dir(matrixFolder))-2 ~= 8 %フォルダが存在しても全結果がない場合は計算する
    doCalculation = true;
else
    doCalculation = false; 
end

newProjectionNumber = 30;
newGridNumber = 50;
% 再構成計算に必要なパラメータを計算するなら読み込む
parameterFile = 'parameters.mat';
if doCalculation
    disp('No matrix data -- Start calculation');
    
    
    if evalin('base', 'exist(''N_projection'', ''var'')')
        NP = evalin('base', 'N_projection');
        if NP ~= newProjectionNumber
            [gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
                      s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);
        end
    else
        [gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
                  s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);
    end

    % 生画像の取得
    rawImage = imread(SXRfilename);
    
    % 非線形フィルターをかける（必要があれば）
    if doFilter
        % figure;imagesc(rawImage);
        disp(size(rawImage));
        [rawImage,~] = imnlmfilt(rawImage,'SearchWindowSize',91,'ComparisonWindowSize',15);
        % figure;imagesc(rawImage);
    end
else
    disp(strcat('Loading matrix from :',matrixFolder))
    load(parameterFile,'range');
end


times = start:interval:(start+interval*7);
doPlot = false;

if doSave
    f = figure;
    f.Units = 'normalized';
    f.Position = [0.1,0.2,0.8,0.8];
end

for t = times
    number = (t-start)/interval+1;
    matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
    if ~exist(matrixPath,'file')%doCalculation
%         ベクトル形式の画像データの読み込み
        [VectorImage1,VectorImage2, VectorImage3, VectorImage4] = get_sxr_image(date,number,newProjectionNumber,rawImage);
        
        datadirPath = getenv('SXR_DATA_DIR');
        dataFolder = strcat(datadirPath,'/',num2str(date),'/shot',num2str(shot));
        if ~exist(dataFolder, 'dir')
            mkdir(dataFolder);
        end

        dataPath = strcat(dataFolder,'/',num2str(number),'.mat');
        n_p = N_projection;
        sxr1 = zeros(n_p);
        sxr2 = zeros(n_p);
        sxr3 = zeros(n_p);
        sxr4 = zeros(n_p);
        k=FindCircle(n_p/2);
        sxr1(k) = VectorImage1;
        sxr2(k) = VectorImage2;
        sxr3(k) = VectorImage3;
        sxr4(k) = VectorImage4;
        save(dataPath, 'sxr1', 'sxr2','sxr3','sxr4')
        

        


%         再構成計算

        EE1 = get_distribution(M,K,gm2d1,U1,s1,v1,VectorImage1,doPlot,ReconMethod);
        EE2 = get_distribution(M,K,gm2d2,U2,s2,v2,VectorImage2,doPlot,ReconMethod);
        EE3 = get_distribution(M,K,gm2d3,U3,s3,v3,VectorImage3,doPlot,ReconMethod);
        EE4 = get_distribution(M,K,gm2d4,U4,s4,v4,VectorImage4,doPlot,ReconMethod);
        
%         再構成結果を保存するファイルを作成、保存
        
        % matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
        save(matrixPath,'EE1','EE2','EE3','EE4');
        
    else
        % matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
        % disp(strcat('Loading result matrix from ',which(matrixPath)));
        load(matrixPath,'EE1','EE2','EE3','EE4');
    end
    
    EE = cat(3,EE1,EE2,EE3,EE4);

    if ~doSave
        f = figure;
        f.Units = 'normalized';
        f.Position = [0.1,0.2,0.8,0.8];
    end

    SXRdata.t = t;
    SXRdata.range = range;

    % plot_save_sxr(grid2D,data2D,range,date,shot,t,EE,show_localmax,show_xpoint,doSave,doFilter,ReconMethod);
    
    if ReconMethod ~= 3
        SXRdata.EE = EE;
        plot_save_sxr(PCBdata,SXR,SXRdata);
    else    
        cGANPath = strcat(dirPath,'/cGAN_large/',num2str(date),'/shot',num2str(shot),'/',num2str(number),'.mat');
        load(cGANPath,'EE1','EE2','EE3','EE4');
        EE = cat(3,EE1,EE2,EE3,EE4);
        SXRdata.EE = EE;
        plot_save_sxr(PCBdata,SXR,SXRdata);
    end


end

if doSave
    close(f);
end

threed = true;
newProjectionNumber3d = 30;
newGridNumber3d = 20;
for t = times
    number = (t-start)/interval+1;
    if threed
        if evalin('base', 'exist(''N_projection'', ''var'')')
            NP = evalin('base', 'N_projection');
            if NP ~= newProjectionNumber3d
                [N_projection, N_grid3d, gm3d,U3d,s3d, v3d, M3d, K3d] = parametercheck_3d(newProjectionNumber3d, newGridNumber3d);
            end
        else
            [N_projection, N_grid3d, gm3d,U3d,s3d, v3d, M3d, K3d] = parametercheck_3d(newProjectionNumber3d, newGridNumber3d);
        end
    
        rawImage = imread(SXRfilename);
        [VectorImage1,VectorImage2, VectorImage3, VectorImage4] = get_sxr_image(date,number,newProjectionNumber3d,rawImage);
    
        l = cat(2,VectorImage1,VectorImage2,VectorImage4); %20240721ではVectorImage3はなし
        
        EE = get_distribution_3d(M3d, K3d, U3d, s3d, v3d, l, N_grid3d);
        EE = permute(EE, [2,1,3]);


        rmin3d = 0; 
        rmax3d = 375; 
        dmin3d = -375; 
        dmax3d = 375; 
        zmin3d = -300; 
        zmax3d = 300;
        r = linspace(rmin3d, rmax3d, newGridNumber3d+1);
        z = linspace(zmin3d, zmax3d, newGridNumber3d+1);
        d = linspace(dmin3d, dmax3d, newGridNumber3d+1);
    
        % Create a grid for the data
        [Z, R, D] = meshgrid(z, r, d);
        
        figure;
        slice(Z, R, D, EE, [-100,100], [100,200], [-100,100]); % Example slice positions along Z

        xlabel('z軸');
        ylabel('r軸');
        zlabel('depth');
        title(strcat('shot',num2str(shot),' at ',num2str(t),'s'));
    
        % Adjust color and viewing angle
        colormap(jet); % Use jet colormap
        colorbar;      % Add a color bar
        view(3);       % Set to 3D view
    end

end
end


function k = FindCircle(L)
    R = zeros(2*L);
    for i = 1:2*L
        for j = 1:2*L
            R(i,j) = sqrt((L-i+0.5)^2+(j-L-0.5)^2);
        end
    end
    % figure;imagesc(R)
    k = find(R<L);
end