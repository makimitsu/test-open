function IdData = get_emission_downstream(PCB,SXR,pathname)

trange = PCB.trange;
[grid2D,data2D] = process_PCBdata_200ch(PCB,pathname);
[magAxisList,xPointList] = get_axis_x_multi(grid2D,data2D);
mergingRatio = xPointList.psi./mean(magAxisList.psi);
trange_sxr = SXR.start:SXR.interval:SXR.start+SXR.interval*7;

mergingRatioSXR = mergingRatio(ismember(trange,trange_sxr));

Imax = zeros(4,8);
Imean = zeros(4,8);
Istd = zeros(4,8);
IdData = struct('max_l',Imax,'mean_l',Imean,'std_l',Istd,'max_r',Imax,'mean_r',Imean,'std_r',Istd,'MR',mergingRatioSXR,'t',trange_sxr);

date = SXR.date;
shot = SXR.shot;
% start = SXR.start;
% interval = SXR.interval;
doFilter = SXR.doFilter;
doNLR = SXR.doNLR;
% SXRfilename = SXR.SXRfilename;

if doFilter & doNLR
    options = 'NLF_NLR';
elseif ~doFilter & doNLR
    options = 'LF_NLR';
elseif doFilter & ~doNLR
    options = 'NLF_LR';
else
    options = 'LF_LR';
end

dirPath = getenv('SXR_MATRIX_DIR');
matrixFolder = strcat(dirPath,filesep,options,filesep,num2str(date),'/shot',num2str(shot));
if exist(matrixFolder,'dir') == 0
    disp('No emission distribution calculation results');
    return
end

parameterFile = 'parameters.mat';
disp(strcat('Loading matrix from :',matrixFolder))
load(parameterFile,'range');

% xPointListに含まれるr,z座標からインデックスを抽出（rq,zqから）
% それぞれの時間の画像でX点座標（インデックス）を中心に5*5で発光強度を抽出
% 最大値、平均値、標準偏差を保存？

range = range./1000;
zmin1 = range(1);
zmax1 = range(2);
zmin2 = range(3);
zmax2 = range(4);
rmin = range(5);
rmax = range(6);
% r_space_SXR = linspace(rmin,rmax,size(EE,1));
% z_space_SXR1 = linspace(zmin1,zmax1,size(EE,2));
% z_space_SXR2 = linspace(zmin2,zmax2,size(EE,2));

% Imax = zeros(4,8);
% Imean = zeros(4,8);
% Istd = zeros(4,8);
% IdData = struct('max_l',Imax,'mean_l',Imean,'std_l',Istd,'max_r',Imax,'mean_r',Imean,'std_r',Istd,'MR',mergingRatioSXR,'t',trange_sxr);

for i = 1:8
    matrixPath = strcat(matrixFolder,filesep,num2str(i),'.mat');
    load(matrixPath,'EE1','EE2','EE3','EE4');
    r_space_SXR = linspace(rmin,rmax,size(EE1,1));
    z_space_SXR1 = linspace(zmin1,zmax1,size(EE1,2));
    z_space_SXR2 = linspace(zmin2,zmax2,size(EE1,2));
    EE = cat(3,EE1,EE2,EE3,EE4);
    EE(EE<0) = 0;
    t = SXR.start+SXR.interval*(i-1);
    t_idx = find(trange==t);
    x_r = xPointList.r(t_idx);
    x_z = xPointList.z(t_idx);
    if ~isnan(x_r)
        for j = 1:4
            % x点の左下・右下にあたる範囲をインデックスで取得したい
            r_idx = knnsearch(r_space_SXR.',x_r);
            if j <= 2
                z_idx = knnsearch(z_space_SXR2.',x_z);
            else
                z_idx = knnsearch(z_space_SXR1.',x_z);
            end
            % 左下の切り出し
            EE_l = EE(1:r_idx,1:z_idx,j);
            % 右下
            EE_r = EE(r_idx:end,z_idx:end,j);
            IdData.max_l(j,i) = max(EE_l,[],'all');
            IdData.mean_l(j,i) = mean(EE_l,'all');
            IdData.std_l(j,i) = std(EE_l,0,'all');
            IdData.max_r(j,i) = max(EE_r,[],'all');
            IdData.mean_r(j,i) = mean(EE_r,'all');
            IdData.std_r(j,i) = std(EE_r,0,'all');
        end
    end

end

end