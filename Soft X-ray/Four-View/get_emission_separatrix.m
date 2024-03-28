function IsepData = get_emission_separatrix(PCB,SXR,pathname)

trange = PCB.trange;
[grid2D,data2D] = process_PCBdata_200ch(PCB,pathname);
[magAxisList,xPointList] = get_axis_x_multi(grid2D,data2D);
mergingRatio = xPointList.psi./mean(magAxisList.psi);
trange_sxr = SXR.start:SXR.interval:SXR.start+SXR.interval*7;

mergingRatioSXR = mergingRatio(ismember(trange,trange_sxr));

Imax_l = zeros(4,8);
Imean_l = zeros(4,8);
Istd_l = zeros(4,8);
Imax_r = zeros(4,8);
Imean_r = zeros(4,8);
Istd_r = zeros(4,8);
IsepData = struct('max_l',Imax_l,'mean_l',Imean_l,'std_l',Istd_l, ...
    'max_r',Imax_r,'mean_r',Imean_r,'std_r',Istd_r, ...
    'MR',mergingRatioSXR,'t',trange_sxr);

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

% range = range./1000;
% zmin1 = range(1);
% zmax1 = range(2);
% zmin2 = range(3);
% zmax2 = range(4);
% rmin = range(5);
% rmax = range(6);
% r_space_SXR = linspace(rmin,rmax,size(EE,1));
% z_space_SXR1 = linspace(zmin1,zmax1,size(EE,2));
% z_space_SXR2 = linspace(zmin2,zmax2,size(EE,2));

% Imax_l = zeros(4,8);
% Imean_l = zeros(4,8);
% Istd_l = zeros(4,8);
% Imax_r = zeros(4,8);
% Imean_r = zeros(4,8);
% Istd_r = zeros(4,8);
% IsepData = struct('max_l',Imax_l,'mean_l',Imean_l,'std_l',Istd_l, ...
%     'max_r',Imax_r,'mean_r',Imean_r,'std_r',Istd_r, ...
%     'MR',mergingRatioSXR,'t',trange_sxr);

for i = 1:8
    matrixPath = strcat(matrixFolder,filesep,num2str(i),'.mat');
    load(matrixPath,'EE1','EE2','EE3','EE4');
    % r_space_SXR = linspace(rmin,rmax,size(EE1,1));
    % z_space_SXR1 = linspace(zmin1,zmax1,size(EE1,2));
    % z_space_SXR2 = linspace(zmin2,zmax2,size(EE1,2));
    EE = cat(3,EE1,EE2,EE3,EE4);
    EE(EE<0) = 0;
    t = SXR.start+SXR.interval*(i-1);
    SXRdata.EE = EE;
    SXRdata.t = t;
    SXRdata.range = range;
    t_idx = find(trange==t);
    x_r = xPointList.r(t_idx);
    x_z = xPointList.z(t_idx);
    x_psi = xPointList.psi(t_idx);
    [data2D_q,zq,rq] = interp_data(grid2D,data2D,SXRdata);
    if ~isnan(x_r)
        % 要電流シートとの切り分け
        % Bp = sqrt((data2D_q.Br).^2+(data2D_q.Bz).^2);
        sep_idx_l = find(((data2D_q.psi<=x_psi+5e-4&data2D_q.psi>=x_psi-5e-5)|(data2D_q.psi<=x_psi*1.1&data2D_q.psi>=x_psi*0.9))&zq<=x_z&rq<=x_r);
        sep_idx_r = find(((data2D_q.psi<=x_psi+5e-4&data2D_q.psi>=x_psi-5e-5)|(data2D_q.psi<=x_psi*1.1&data2D_q.psi>=x_psi*0.9))&zq>=x_z&rq<=x_r);
        % sep_idx_l = find(data2D_q.psi<=x_psi*1.1&data2D_q.psi>=x_psi*0.9&zq<=x_z&rq<=x_r&Bp>=0.02);
        % sep_idx_r = find(data2D_q.psi<=x_psi*1.1&data2D_q.psi>=x_psi*0.9&zq>=x_z&rq<=x_r&Bp>=0.02);
        if isempty(sep_idx_r) || isempty(sep_idx_l)
            continue
        end
        % A = zeros(size(zq));
        % A(sep_idx_l) = 1;
        % B = zeros(size(zq));
        % B(sep_idx_r) = 1;        
        % figure;
        % subplot(1,2,1);contourf(zq,rq,A);hold on;contour(zq,rq,data2D_q.psi,'Fill','off');
        % subplot(1,2,2);contourf(zq,rq,B);hold on;contour(zq,rq,data2D_q.psi,'Fill','off');
        for j = 1:4
            EE_j = data2D_q.EE(:,:,j);
            IsepData.max_l(j,i) = max(EE_j(sep_idx_l),[],'all');
            IsepData.mean_l(j,i) = mean(EE_j(sep_idx_l),'all');
            IsepData.std_l(j,i) = std(EE_j(sep_idx_l),0,'all');
            IsepData.max_r(j,i) = max(EE_j(sep_idx_r),[],'all');
            IsepData.mean_r(j,i) = mean(EE_j(sep_idx_r),'all');
            IsepData.std_r(j,i) = std(EE_j(sep_idx_r),0,'all');
        end
    end
end
end