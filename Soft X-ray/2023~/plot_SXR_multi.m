function [] = plot_SXR_multi(B_z,r_probe,z_probe,date,shot,show_xpoint,show_localmax,start,interval,save,SXRfilename,filter,NL)
% plot SXR emission on psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   integer: date, date of experiment
%   integer: shot, number of shot
%   boolean: show_xpoint, option for showing the x-point
%   boolean: show_localmax, option for showing the local maximum point
%   integer: start, start time (us)
%   integer: interval, interval time of the framing camera (us)
%   boolean: save, option for saving the reconstruction result
%   string: SXRfilename, name of the SXR image file
%   boolean: filter, option for applying non-linear mean (NLM) filter
%   boolean: NL, option for using non-linear reconstruction

% 実行結果（行列）を保存するフォルダの確認
% なければ作成＆計算、あれば読み込みsave
% if filter & NL
%     savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/NLF_NLR/',num2str(date),'/shot',num2str(shot));
% elseif ~filter & NL
%     savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/LF_NLR/',num2str(date),'/shot',num2str(shot));
% elseif filter & ~NL
%     savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/NLF_LR/',num2str(date),'/shot',num2str(shot));
% else
%     savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/LF_LR/',num2str(date),'/shot',num2str(shot));
% end
% savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/',num2str(date),'/shot',num2str(shot));

if filter & NL
    options = 'NLF_NLR';
elseif ~filter & NL
    options = 'LF_NLR';
elseif filter & ~NL
    options = 'NLF_LR';
else
    options = 'LF_LR';
end
savefolder = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/result_matrix/',options,'/',num2str(date),'/shot',num2str(shot));
if exist(savefolder,'dir') == 0
    clc_flag = true;
    mkdir(savefolder);
else
    clc_flag = false;
end

% 再構成計算に必要なパラメータを計算するなら読み込む
filepath = '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/2023~/parameters.mat';
if clc_flag
    N_projection_new = 80;
    N_grid_new = 100;
    if isfile(filepath)
        load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range','N_projection', 'N_grid');
        if N_projection_new ~= N_projection || N_grid_new ~= N_grid
            disp('Different parameters - Start calculation!');
            clc_parameters(N_projection_new,N_grid_new,filepath);
            load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range');
        end
    else
        disp('No parameter - Start calculation!');
        clc_parameters(N_projection_new,N_grid_new,filepath);
        load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range');
    end
else
    load(filepath,'range');
end

times = start:interval:(start+interval*7);
plot_flag = false;

for t = times
    number = (t-start)/interval+1;
    
    if clc_flag
%         ベクトル形式の画像データの読み込み
        if date <= 210924
            [VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename,filter);
        else
            [VectorImage2,VectorImage1] = get_SXRImage(date,number,SXRfilename,filter);
        end
        
%         再構成計算
        EE_high = clc_distribution(M,K,gm2d1,U1,s1,v1,VectorImage1,plot_flag,NL);
        EE_low = clc_distribution(M,K,gm2d2,U2,s2,v2,VectorImage2,plot_flag,NL);
        
%         再構成結果を保存するファイルを作成、保存
        savepath_high = strcat(savefolder,'/',num2str(number),'_high.txt');
        writematrix(EE_high,savepath_high);
        savepath_low = strcat(savefolder,'/',num2str(number),'_low.txt');
        writematrix(EE_low,savepath_low);
    else
        loadpath_high = strcat(savefolder,'/',num2str(number),'_high.txt');
        loadpath_low = strcat(savefolder,'/',num2str(number),'_low.txt');
        EE_high = readmatrix(loadpath_high);
        EE_low = readmatrix(loadpath_low);
    end
    
    plot_save_SXR(B_z,r_probe,z_probe,range,date,shot,t,EE_high,EE_low,show_localmax,show_xpoint,save,filter,NL);

end

end