clearvars -except date IDXlist times

clearvars -except date IDXlist doSave doFilter doNLR
addpath '/Users/shohgookazaki/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
addpath '/Users/shohgookazaki/Documents/GitHub/test-open'/'Soft X-ray'/Four-View; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
% ~/Documents/MATLAB にてstartup.mを作って、その中でsetenv('パス名','アドレス')していくと自動になる。
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.pre_processed_directory_path=getenv('pre_processed_directory_path');


%%%%実験オペレーションの取得
prompt = {'Date:','Shot number:','times:'};
definput = {'','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
if exist('IDXlist','var')
    definput{2} = num2str(IDXlist);
end
if exist('times','var')
    definput{3} = num2str(times);
end
dlgtitle = 'Input';
dims = [1 35];
% if exist('date','var') && exist('IDX','var') && exist('times','var')
%     definput = {num2str(date),num2str(IDX),num2str(times)};
% else
%     definput = {'','',''};
% end
% definput = {'','',''};
% definput = {num2str(date),num2str(IDXlist),num2str(doCheck)};
answer = inputdlg(prompt,dlgtitle,dims,definput);
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2))); 
times = str2num(cell2mat(answer(3))); 


%-----------スプレッドシートからデータ抜き取り--------------------%
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);

T=searchlog(T,'date',date);
if isnan(T.shot(1))
    T(1, :) = [];
end
n_data=numel(IDXlist);%計測データ
shotlist = [T.a039(IDXlist), T.a040(IDXlist)];
tfshotlist = [T.a039_TF(IDXlist), T.a040_TF(IDXlist)];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1); % 39が計測データ数だけ縦に並ぶ。
startlist = T.SXRStart(IDXlist);
intervallist = T.SXRInterval(IDXlist);


PCB.trange=400:800;%【input】計算時間範囲
PCB.n=50; %【input】rz方向のメッシュ数
PCB.restart = 0;

figure;hold on
xlabel('time [us]');ylabel('Merging ratio [%]');
legendList = cell(1,n_data);
ax=gca;ax.FontSize=18;

% % エクセルファイルの保存先とファイル名を指定
% outputFile = 'merging_rate.xlsx';
% 
% % 書き込み対象のデータを初期化
% all_merging_ratios = [];
% % 最長のmerging_ratioの長さを記録する変数
% max_length = 0;


for i = 1:n_data
    shot = IDXlist(i);
    start = startlist(i);
    interval = intervallist(i);

    PCB.idx = IDXlist(i);
    PCB.shot=shotlist(i,:);
    PCB.tfshot=tfshotlist(i,:);
    if PCB.shot == PCB.tfshot
        PCB.tfshot = [0,0];
    end
    PCB.i_EF=EFlist(i);
    PCB.date = date;
    TF=TFlist(i);
    [grid2D,data2D] = process_PCBdata_280ch(PCB,pathname); %process_PCBdata_200ch.mに行く
    if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
        return
    end
    
    merging_ratio = get_merging_ratio(data2D,grid2D,times);
    legendList(i) = cellstr(strcat('shot',num2str(IDXlist(i))));
    
    % % 最長のmerging_ratioの長さを更新
    % max_length = max(max_length, length(merging_ratio));
    % % merging_ratioのデータを保存用の配列に追加
    % all_merging_ratios = [all_merging_ratios; merging_ratio];

    % figure;
    plot(times,merging_ratio,'LineWidth',2);
    % xlabel('time [us]');ylabel('Merging ratio [%]');
    % title('Merging ratio');
    % ax=gca;ax.FontSize=18;

end
legend(legendList,'Location','northwest');

% % Initialize padded array with NaN and set proper dimensions
% padded_merging_ratios = NaN(n_data, max_length);
% 
% for i = 1:n_data
%     % Pad each row of all_merging_ratios to match max_length
%     padded_merging_ratios(i, 1:length(all_merging_ratios(i, :))) = all_merging_ratios(i, :);
% end
% 
% % Verify sizes of IDXlist and padded_merging_ratios for concatenation
% if length(IDXlist) ~= size(padded_merging_ratios, 1)
%     error('Length of IDXlist (%d) does not match the number of rows in padded_merging_ratios (%d).', length(IDXlist), size(padded_merging_ratios, 1));
% end
% 
% % Concatenate IDXlist with padded_merging_ratios
% outputData = [IDXlist, padded_merging_ratios];
% 
% % Create header
% header = [{'Shot Number'}, arrayfun(@(t) sprintf('Time_%dus', t), times(1:max_length), 'UniformOutput', false)];
% 
% % Write to output file
% outputData = [header; num2cell(outputData)];
% writecell(outputData, outputFile);








function merging_ratio = get_merging_ratio(data2D,grid2D,times)

merging_ratio = zeros(numel(times),1);
i = 1;

for time = times
    psi_timeseries = data2D.psi;
    psi = psi_timeseries(:,:,data2D.trange==time);
    merging_ratio(i) = clc_merging_ratio(psi,grid2D);
    i = i+1;
end
idx_nan=find(isnan(merging_ratio));
for j = 1:numel(idx_nan)
    i = idx_nan(j);
    if times(i)<470
        merging_ratio(i) = 0;
    elseif times(i)>480
        merging_ratio(i) = 1;
    end
end
merging_ratio(merging_ratio<0)=0;
merging_ratio=merging_ratio*100;
% figure('Position',[100 100 600 200]);
% if plot_flag
%     figure;
%     plot(times,merging_ratio,'k-','LineWidth',2);
%     xlabel('time [us]');ylabel('Merging ratio [%]');
%     title('Merging ratio');
%     ax=gca;ax.FontSize=18;
% end

% merging_speed = diff(merging_ratio);
% figure;
% plot(times(2:end),merging_speed,'k-','LineWidth',2);
% xlabel('time [us]');ylabel('Merging speed [a.u.]');
% title('Merging speed');
% ax=gca;ax.FontSize=18;

end

function merging_ratio = clc_merging_ratio(psi,grid2D)

[pos_oz,pos_or,pos_xz,~,psi_o,psi_x] = search_xo(psi,grid2D);

% if sum(pos_oz<pos_xz)==1 && sum(isnan(pos_oz))==0 % pos_xzがpos_ozの最大値と最小値の間にあることを判別＆O点が2つあることを判別
% 必要な条件：O点が少なくとも一つある＋X点が（あれば）その二つの間
if max(sum(pos_oz<pos_xz),sum(pos_oz>pos_xz))==1 % pos_xzがpos_ozの最大値と最小値の間にあることを判別（NaNが入ってもいいように）
    common_flux = psi_x;
    private_flux = mean(psi_o,'omitnan');
elseif sum(isnan(pos_oz))==0 && sum(pos_oz>0)==1 && isnan(pos_xz)% O点二つ（値が違う）あるけどX点なし
    pos_xz = mean(pos_oz);
    pos_xr = mean(pos_or);
    r_space = grid2D.rq(:,1);
    z_space = grid2D.zq(1,:);
    z_idx = knnsearch(z_space',pos_xz);
    r_idx = knnsearch(r_space,pos_xr);
    common_flux = psi(r_idx,z_idx);
    private_flux = min(psi_o);
else
    common_flux = NaN;
    private_flux = NaN;
end

merging_ratio = common_flux/private_flux;

end

function [pos_oz,pos_or,pos_xz,pos_xr,psi_o,psi_x] = search_xo(psi,grid2D)

r_space = grid2D.rq(:,1);
z_space = grid2D.zq(1,:);
mesh_r = numel(r_space);

% psi_or=zeros(2,length(time));
pos_or=zeros(2,1);
pos_oz=pos_or;
psi_o = zeros(2,1);
% psi_xr=zeros(1,length(time));
% psi_xz=psi_xr;
% [max_psi,max_psi_r]=max(psi_store,[],1); %各時間、各列(z)ごとのpsiの最大値
[max_psi,max_psi_r]=max(psi,[],1); %各列(z)ごとのpsiの最大値
max_psi=squeeze(max_psi);

% for i = time
%     ind=i-time(1)+1;
%     max_psi_ind=find(islocalmax(smooth(max_psi(:,i)),'MaxNumExtrema', 2));
    max_psi_ind=find(islocalmax(smooth(max_psi),'MaxNumExtrema', 2));
%     r_ind=max_psi_r(:,:,i);
    r_ind=max_psi_r;
%     if numel(max_psi(min(max_psi_ind),i))==0
    if numel(max_psi(min(max_psi_ind)))==0
%         psi_or(1,ind)=NaN;
%         psi_oz(1,ind)=NaN;
        pos_or(1)=NaN;
        pos_oz(1)=NaN;
        psi_o(1) = NaN;
    else
        o1r=r_ind(min(max_psi_ind));
%         psi_or(1,ind)=r_space(o1r);
%         psi_oz(1,ind)=z_space(min(max_psi_ind));
        pos_or(1)=r_space(o1r);
        pos_oz(1)=z_space(min(max_psi_ind));
%         psi_o(1) = psi_store(min(max_psi_ind),o1r);
        psi_o(1) = psi(o1r,min(max_psi_ind));
    end
%     if numel(max_psi(max(max_psi_ind),i))==0
    if numel(max_psi(max(max_psi_ind)))==0
%         psi_or(2,ind)=NaN;
%         psi_oz(2,ind)=NaN;
        pos_or(2)=NaN;
        pos_oz(2)=NaN;
        psi_o(2) = NaN;
    else
        o2r=r_ind(max(max_psi_ind));
%         psi_or(2,ind)=r_space(o2r);
%         psi_oz(2,ind)=z_space(max(max_psi_ind));
        pos_or(2)=r_space(o2r);
        pos_oz(2)=z_space(max(max_psi_ind));
%         psi_o(2) = psi_store(min(max_psi_ind),o2r);
        psi_o(2) = psi(o2r,max(max_psi_ind));
    end
%     if numel(find(islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1)))==0
    if numel(find(islocalmin(smooth(max_psi),'MaxNumExtrema', 1)))==0
%         psi_xr(1,ind)=NaN;
%         psi_xz(1,ind)=NaN;
        pos_xr=NaN;
        pos_xz=NaN;
        psi_x = NaN;
    else
%         min_psi_ind=islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1);
        min_psi_ind=islocalmin(smooth(max_psi),'MaxNumExtrema', 1);
        xr=r_ind(min_psi_ind);
        if xr==1 || xr==mesh_r %r両端の場合は検知しない
%             psi_xr(1,ind)=NaN;
%             psi_xz(1,ind)=NaN;
            pos_xr=NaN;
            pos_xz=NaN;
            psi_x = NaN;
        else
%             psi_xr(1,ind)=r_space(xr);
%             psi_xz(1,ind)=z_space(min_psi_ind);
            pos_xr=r_space(xr);
            pos_xz=z_space(min_psi_ind);
%             psi_x = psi_store(min_psi_ind,xr);
            psi_x = psi(xr,min_psi_ind);
        end
    end
%     if max_psi(1,i)==max(max_psi(1:end/2,i)) %z両端の場合はpsi中心が画面外として検知しない
    mid_idx = round(numel(max_psi)/2);
    if max_psi(1)==max(max_psi(1:mid_idx)) %z両端の場合はpsi中心が画面外として検知しない
%         psi_or(1,ind)=NaN; %r_space(r_ind(1));
%         psi_oz(1,ind)=NaN; %z_space(1);
        pos_or(1)=NaN; %r_space(r_ind(1));
        pos_oz(1)=NaN; %z_space(1);
        psi_o(1) = NaN;
    end
%     if max_psi(end,i)==max(max_psi(end/2:end,i))
    if max_psi(end)==max(max_psi(mid_idx:end))
%         psi_or(2,ind)=NaN; %r_space(r_ind(1));
%         psi_oz(2,ind)=NaN; %z_space(end);
        pos_or(2)=NaN; %r_space(r_ind(1));
        pos_oz(2)=NaN; %z_space(end);
        psi_o(2) = NaN;
    end
end