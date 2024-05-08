% function plot_sxr_at_t(grid2D,data2D,date,shot,t,show_xpoint,show_localmax,start,interval,save,SXRfilename,filter,NL)
function plot_sxr_at_t(PCBdata,SXR)

% grid2D = PCBdata.grid2D;
% data2D = PCBdata.data2D;
date = SXR.date;
shot = SXR.shot;
% show_xpoint = SXR.show_xpoint;
% show_localmax = SXR.show_localmax;
start = SXR.start;
interval = SXR.interval;
doFilter = SXR.doFilter;
doNLR = SXR.doNLR;
SXRfilename = SXR.SXRfilename;
t = SXR.t;
doCheck = SXR.doCheck;


% save = false;

newProjectionNumber = 50;
newGridNumber = 90;

if doFilter & doNLR
    options = 'NLF_NLR';
elseif ~doFilter & doNLR
    options = 'LF_NLR';
elseif doFilter & ~doNLR
    options = 'NLF_LR';
else
    options = 'LF_LR';
end
matrixFolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/' ...
    ,options,'/',num2str(date),'/shot',num2str(shot));

if exist(matrixFolder,'dir') == 0 || doCheck
    if ~doCheck
        mkdir(matrixFolder)
    end
    doCalculation = true;
else
    doCalculation = false;
end

% 再構成計算に必要なパラメータを計算するなら読み込む、しない場合も範囲に関しては読み込む
parameterFile = 'parameters.mat';

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

    % if ~doCheck
        % 生画像の取得
        rawImage = imread(SXRfilename);
        
        % 非線形フィルターをかける（必要があれば）
        if doFilter
            % figure;imagesc(rawImage);
            [rawImage,~] = imnlmfilt(rawImage,'SearchWindowSize',91,'ComparisonWindowSize',15);
            % figure;imagesc(rawImage);
        end
    % end
else
    load(parameterFile,'range');
end
    
number = (t-start)/interval+1;
doPlot = false;

if doCheck
    [Iwgn_old1,Iwgn_old2,Iwgn_old3,Iwgn_old4] = get_sxr_image(date,number,newProjectionNumber,rawImage);
    get_gamma_min(M,U1,U2,U3,U4,s1,s2,s3,s4,Iwgn_old1,Iwgn_old2,Iwgn_old3,Iwgn_old4);

    doPlot = false;
    matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
    load(matrixPath,'EE1','EE2','EE3','EE4');
    E1 = reshape(flipud(EE1),[],1);
    E2 = reshape(flipud(EE2),[],1);
    E3 = reshape(flipud(EE3),[],1);
    E4 = reshape(flipud(EE4),[],1);

    % f = figure;
    % f.Units = 'normalized';
    % f.Position = [0.1,0.2,0.8,0.8];
    % EE = cat(3,EE1,EE2,EE3,EE4);
    % plot_save_sxr(grid2D,data2D,range,date,shot,t,EE,show_localmax,show_xpoint,save,doFilter,doNLR);  

    Iwgn1 = (gm2d1*E1).';
    Iwgn2 = (gm2d2*E2).';
    Iwgn3 = (gm2d3*E3).';
    Iwgn4 = (gm2d4*E4).';

    % Iwgn1 = awgn((gm2d1*E1).',10,'measured');
    % Iwgn2 = awgn((gm2d2*E2).',10,'measured');
    % Iwgn3 = awgn((gm2d3*E3).',10,'measured');
    % Iwgn4 = awgn((gm2d4*E4).',10,'measured');

    get_gamma_min(M,U1,U2,U3,U4,s1,s2,s3,s4,Iwgn1,Iwgn2,Iwgn3,Iwgn4);

    figure;
    subplot(2,2,1);plot(Iwgn_old1);hold on;plot(Iwgn1);
    subplot(2,2,2);plot(Iwgn_old2);hold on;plot(Iwgn2);
    subplot(2,2,3);plot(Iwgn_old3);hold on;plot(Iwgn3);
    subplot(2,2,4);plot(Iwgn_old4);hold on;plot(Iwgn4);

    EE_new1 = get_distribution(M,K,gm2d1,U1,s1,v1,Iwgn1,doPlot,doNLR);
    EE_new2 = get_distribution(M,K,gm2d2,U2,s2,v2,Iwgn2,doPlot,doNLR);
    EE_new3 = get_distribution(M,K,gm2d3,U3,s3,v3,Iwgn3,doPlot,doNLR);
    EE_new4 = get_distribution(M,K,gm2d4,U4,s4,v4,Iwgn4,doPlot,doNLR);

    f = figure;
    f.Units = 'normalized';
    f.Position = [0.1,0.2,0.8,0.8];
    EE_new = cat(3,EE_new1,EE_new2,EE_new3,EE_new4);
    SXRdata.EE = EE_new;
    SXRdata.t = t;
    SXRdata.range = range;
    % plot_save_sxr(grid2D,data2D,range,date,shot,t,EE_new,show_localmax,show_xpoint,save,doFilter,doNLR);
    plot_save_sxr(PCBdata,SXR,SXRdata);

    error1 = sum((EE_new1-EE1).^2/max(EE1,[],'all'),'all')/numel(EE1);
    error2 = sum((EE_new2-EE2).^2/max(EE2,[],'all'),'all')/numel(EE2);
    error3 = sum((EE_new3-EE3).^2/max(EE3,[],'all'),'all')/numel(EE3);
    error4 = sum((EE_new4-EE4).^2/max(EE4,[],'all'),'all')/numel(EE4);

    disp(error1);
    disp(error2);
    disp(error3);
    disp(error4);

    return
end

matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
if doCalculation

    [Iwgn1,Iwgn2,Iwgn3,Iwgn4] = get_sxr_image(date,number,newProjectionNumber,rawImage);
    
    EE1 = get_distribution(M,K,gm2d1,U1,s1,v1,Iwgn1,doPlot,doNLR);
    EE2 = get_distribution(M,K,gm2d2,U2,s2,v2,Iwgn2,doPlot,doNLR);
    EE3 = get_distribution(M,K,gm2d3,U3,s3,v3,Iwgn3,doPlot,doNLR);
    EE4 = get_distribution(M,K,gm2d4,U4,s4,v4,Iwgn4,doPlot,doNLR);
    save(matrixPath,'EE1','EE2','EE3','EE4');
else
    % matrixPath = strcat(matrixFolder,'/',num2str(number),'.mat');
    load(matrixPath,'EE1','EE2','EE3','EE4');
end

f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.2,0.8,0.8];

EE = cat(3,EE1,EE2,EE3,EE4);

SXRdata.EE = EE;
SXRdata.t = t;
SXRdata.range = range;

% plot_save_sxr(grid2D,data2D,range,date,shot,t,EE,show_localmax,show_xpoint,save,doFilter,doNLR);
plot_save_sxr(PCBdata,SXR,SXRdata);

end

function get_gamma_min(M,U1,U2,U3,U4,s1,s2,s3,s4,Iwgn1,Iwgn2,Iwgn3,Iwgn4)
plot_flag = true;
gmin=-15;gmax=15;
% gmin=-30;gmax=0;
lg_gamma=gmin:1:gmax;
l_g = numel(lg_gamma);
gamma=10.^(lg_gamma);

Z1=U1'*Iwgn1.';
Vgamma_1 = zeros(1,l_g);
Z2=U2'*Iwgn2.';
Vgamma_2 = zeros(1,l_g);
Z3=U3'*Iwgn3.';
Vgamma_3 = zeros(1,l_g);
Z4=U4'*Iwgn4.';
Vgamma_4 = zeros(1,l_g);

for n=1:l_g
    rho1 = M*gamma(n)./(s1.^2+M*gamma(n));
    Vgamma_1(n)=(M*sum((rho1.*(Z1.')).^2))/(sum(rho1))^2;
    rho2 = M*gamma(n)./(s2.^2+M*gamma(n));
    Vgamma_2(n)=(M*sum((rho2.*(Z2.')).^2))/(sum(rho2))^2;
    rho3 = M*gamma(n)./(s3.^2+M*gamma(n));
    Vgamma_3(n)=(M*sum((rho3.*(Z3.')).^2))/(sum(rho3))^2;
    rho4 = M*gamma(n)./(s4.^2+M*gamma(n));
    Vgamma_4(n)=(M*sum((rho4.*(Z4.')).^2))/(sum(rho4))^2;
end

if plot_flag
    figure;
    subplot(2,2,1);semilogy(lg_gamma,Vgamma_1,'*');
    xlabel('logγ');ylabel('GCV');
    subplot(2,2,2);semilogy(lg_gamma,Vgamma_2,'*');
    xlabel('logγ');ylabel('GCV');
    subplot(2,2,3);semilogy(lg_gamma,Vgamma_3,'*');
    xlabel('logγ');ylabel('GCV');
    subplot(2,2,4);semilogy(lg_gamma,Vgamma_4,'*');
    xlabel('logγ');ylabel('GCV');
end

[~,gamma_index1]=min(Vgamma_1);
disp(gamma_index1);
[~,gamma_index2]=min(Vgamma_2);
disp(gamma_index2);
[~,gamma_index3]=min(Vgamma_3);
disp(gamma_index3);
[~,gamma_index4]=min(Vgamma_4);
disp(gamma_index4);
end