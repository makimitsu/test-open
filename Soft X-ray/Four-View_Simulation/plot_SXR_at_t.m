function [EE_high,EE_low] = plot_SXR_at_t(grid2D,data2D,date,shot,t,show_xpoint,show_localmax,start,interval,save,SXRfilename,filter,NL)
% plot SXR emission on psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   integer: date, date of experiment
%   integer: shot, number of shot
%   integer: t, time of interest (us)
%   boolean: show_xpoint, option for showing the x-point
%   boolean: show_localmax, option for showing the local maximum point
%   integer: start, start time (us)
%   integer: interval, interval time of the framing camera (us)
%   boolean: save, option for saving the reconstruction result
%   string: SXRfilename, name of the SXR image file
%   boolean: filter, option for applying non-linear mean (NLM) filter
%   boolean: NL, option for using non-linear reconstruction

% if filter & NL
%     options = 'NLF_NLR';
% elseif ~filter & NL
%     options = 'LF_NLR';
% elseif filter & ~NL
%     options = 'NLF_LR';
% else
%     options = 'LF_LR';
% end
% savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/',options,'/',num2str(date),'/shot',num2str(shot));
% 非線形再構成では計算精度を落としたい→グリッド数を落とす
% 投影数は変えない？非線形フィルタなら落とす？
% グリッド数以下の投影数だと正則化する意味もなさそう
% もっと投影数落としていいのかな→その分ノイズ落とす？
% 光ファイバーの本数にしてるけどそこまでする必要もなさそう
% savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/',num2str(date),'/shot',num2str(shot));
% savefolder = strcat('/Users/shinjirotakeda/Documents/Github/SXR_diagnostics/result_matrix/',num2str(date),'/shot',num2str(shot));
% savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/',options,'/',num2str(date),'/shot',num2str(shot));
% filepath = strcat('/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/parameters/',options,'.mat');
% 
% if exist(savefolder,'dir') == 0
%     clc_flag = true;
% %     mkdir(savefolder);
% else
%     clc_flag = false;
% end

% 再構成計算に必要なパラメータを計算するなら読み込む、しない場合も範囲に関しては読み込む
filepath = '/Users/yuleo/Documents/GitHub/test-open/Soft X-ray/Four-View_Simulation/parameters.mat';
% if clc_flag
%     N_projection_new = 80;
%     N_grid_new = 100;
%     if isfile(filepath)
%         load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range','N_projection', 'N_grid');
% %         load(filepath);
%         if N_projection_new ~= N_projection || N_grid_new ~= N_grid
%             disp('Different parameters - Start calculation!');
%             clc_parameters(N_projection_new,N_grid_new,filepath);
%             load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range');
%         end
%     else
%         disp('No parameter - Start calculation!');
%         clc_parameters(N_projection_new,N_grid_new,filepath);
%         load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range');
%     end
% else
%     load(filepath,'range');
% end
% load(filepath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
    % 's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
if isfile(filepath)
    load(filepath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
        's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
    if N_projection_new ~= N_projection || N_grid_new ~= N_grid
        disp('Different parameters - Start calculation!');
        clc_parameters(N_projection_new,N_grid_new,filepath);
        load(filepath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
            's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
    end
else
    disp('No parameter - Start calculation!');
    clc_parameters(N_projection_new,N_grid_new,filepath);
    load(filepath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
        's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');   
end
    
% number = (t-start)/interval+1;
plot_flag = false;
[~,Iwgn] = Assumption(N_projection,gm2d,true);
EE = clc_distribution(M,K,gm2d1,U1,s1,v1,Iwgn,plot_flag,NL);
% [VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename,filter);
% EE_high = clc_distribution(M,K,gm2d1,U1,s1,v1,VectorImage1,plot_flag,NL);
% plot_save_SXR(grid2D,data2D,range,date,shot,t,EE_high,EE_low,show_localmax,show_xpoint,save,filter,NL);

% if clc_flag
% %         ベクトル形式の画像データの読み込み
% %     if date <= 210924
% %         [VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename,filter);
% % %         return
% %     else
% %         [VectorImage2,VectorImage1] = get_SXRImage(date,number,SXRfilename,filter);
% %     end
% [VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename,filter);
% 
% %         再構成計算
%     EE_high = clc_distribution(M,K,gm2d1,U1,s1,v1,VectorImage1,plot_flag,NL);
%     EE_low = clc_distribution(M,K,gm2d2,U2,s2,v2,VectorImage2,plot_flag,NL);
% else
%     loadpath_high = strcat(savefolder,'/',num2str(number),'_high.txt');
%     loadpath_low = strcat(savefolder,'/',num2str(number),'_low.txt');
%     EE_high = readmatrix(loadpath_high);
%     EE_low = readmatrix(loadpath_low);
% end

% subplot(1,2,1);
% imagesc(EE_high);
% subplot(1,2,2);
% imagesc(EE_low);

f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.2,0.8,0.4];

plot_save_SXR(grid2D,data2D,range,date,shot,t,EE_high,EE_low,show_localmax,show_xpoint,save,filter,NL);

end