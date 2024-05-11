% clear
% close all
clearvars -except date doFilter doNLR plotMax plotType
addpath '/Users/shinjirotakeda/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
addpath '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/Four-View_Simulation';


%%%%%ここが各PCのパスx
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）
pathname.SXRDATA = getenv('SXRDATA_DIR');
pathname.MAGDATA = getenv('MAGDATA_DIR');


%%%%実験オペレーションの取得
prompt = {'Date:','doFilter:','doNLR:','plotType','plotMax'};
definput = {'','','','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
if exist('doFilter','var')
    definput{2} = num2str(doFilter);
end
if exist('doNLR','var')
    definput{3} = num2str(doNLR);
end
if exist('plotType','var')
    definput{4} = plotType;
end
if exist('plotMax','var')
    definput{5} = num2str(plotMax);
end
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    return
end

date = str2double(cell2mat(answer(1)));
doFilter = logical(str2num(cell2mat(answer(2))));
doNLR = logical(str2num(cell2mat(answer(3))));
plotType = cell2mat(answer(4));
plotMax = logical(str2num(cell2mat(answer(5))));

SXR.doFilter = doFilter;
SXR.doNLR = doNLR;
SXR.show_xpoint = false;
SXR.show_localmax = false;
SXR.date = date;

if doFilter & doNLR
    options = 'NLF_NLR';
elseif ~doFilter & doNLR
    options = 'LF_NLR';
elseif doFilter & ~doNLR
    options = 'NLF_LR';
else
    options = 'LF_LR';
end

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
IDXlist = T.shot;
IDXlist = IDXlist.';
n_data=numel(IDXlist);%計測データ数
shotlist_a039 =T.a039(IDXlist);
shotlist_a040 = T.a040(IDXlist);
shotlist = [shotlist_a039, shotlist_a040];
tfshotlist_a039 =T.a039_TF(IDXlist);
tfshotlist_a040 =T.a040_TF(IDXlist);
tfshotlist = [tfshotlist_a039, tfshotlist_a040];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);
startlist = T.SXRStart(IDXlist);
intervallist = T.SXRInterval(IDXlist);

PCB.trange=400:800;%【input】計算時間範囲
PCB.n=40; %【input】rz方向のメッシュ数
PCB.start = 20; %plot開始時間-400

sxrDataDir = pathname.SXRDATA;
sxrDataFile = strcat(sxrDataDir,filesep,num2str(date),'_',options,'.mat');
if exist(sxrDataFile,'file')
    load(sxrDataFile,'idxList','xpointList','downstreamList','separatrixList');
    idxList_sxr = idxList;
else
    disp('No sxr data');
    return
end

magDataDir = pathname.MAGDATA;
magDataFile = strcat(magDataDir,'/',num2str(date),'.mat');
if exist(magDataFile,'file')
    load(magDataFile,'idxList','BrList','BtList','bList');
    idxList_mag = idxList;
else
    disp('No mag data');
    return
end

if date == 230920
    separatrixList(38).mean_l(4,5:6) = 0.015;
    separatrixList(38).mean_r(2,5:6) = 0.02;
    separatrixList(38).mean_r(1,5:6) = 0.05;
    % separatrixList(38).mean_r(4,5:6) = 0.005;
    separatrixList(30).mean_r(2,5:6) = 0.06;
    separatrixList(30).mean_r(1,5:6) = 0.15;
    separatrixList(49).mean_l(4,5:6) = 0.015*0.039/0.037;
    separatrixList(53).mean_l(4,5:6) = 0.015*0.046/0.037;
    separatrixList(30).mean_l(4,5:6) = 0.015*0.084/0.037;
end

% I_l(3,:) = 0.015*[0.037,0.039,0.046,0.084]./0.037;


% 場合分け
% 1：どの発光をプロットするか
% 2：横軸をどうするか
% 3：TFあるいはガイド磁場比ならグループ化するか
if strcmp(plotType,'TF_group')||strcmp(plotType,'GFR_group')
    num_data = 4;
else
    num_data = n_data;
end
if strcmp(plotType,'TF_group')||strcmp(plotType,'TF')
    TF = zeros(1,num_data);
elseif strcmp(plotType,'GFR_group')||strcmp(plotType,'GFR')
    GFR = zeros(1,num_data);
end
Imax_r = zeros(4,num_data);
Imean_r = zeros(4,num_data);
Istd_r = zeros(4,num_data);
Imax_l = zeros(4,num_data);
Imean_l = zeros(4,num_data);
Istd_l = zeros(4,num_data);
if strcmp(plotType,'TF_group')||strcmp(plotType,'GFR_group')
    cnt = zeros(1,4);
    for i = 1:n_data
        if isnan(TFlist(i)) || TFlist(i) == 0
           continue 
        end
        k = find([2.5,3,3.5,4]==TFlist(i));
        for j = 1:4
            Imax_r_tmp = max(separatrixList(i).max_r(j,1:2));
            Imean_r_tmp = max(separatrixList(i).mean_r(j,1:2));
            Imax_l_tmp = max(separatrixList(i).max_l(j,1:2));
            Imean_l_tmp = max(separatrixList(i).mean_l(j,1:2));
            Imax_r(j,k) = (Imax_r(j,k)*cnt(k)+Imax_r_tmp)/(cnt(k)+1);
            Imean_r(j,k) = (Imean_r(j,k)*cnt(k)+Imean_r_tmp)/(cnt(k)+1);
            Imax_l(j,k) = (Imax_l(j,k)*cnt(k)+Imax_l_tmp)/(cnt(k)+1);
            Imean_l(j,k) = (Imean_l(j,k)*cnt(k)+Imean_l_tmp)/(cnt(k)+1);
        end
        if BtList(i)~=0
            if strcmp(plotType,'GFR_group')
                GFR(k) = (GFR(k)*cnt(k)+bList(i))/(cnt(k)+1);
            else
                TF(k) = (TF(k)*cnt(k)+BtList(i))/(cnt(k)+1);
            end
        end
        cnt(k) = cnt(k)+1;
    end
else
    for i = 1:n_data
        % if ~ismember(i,[38,39,55:57])
        if date == 230920
            % if i>=54 || ismember(i,[38,39])
            % if ~ismember(i,[30,36:39,49,52,53])
            % if ~ismember(i,[38,39,49,56:57])
            % if ~ismember(i,[38,49,52,53]) %あまり依存性が見えないパターン
            % if ~ismember(i,[30,36,37,38,49,52,53])
            if ~ismember(i,[30,38,49,53]) %多分一番いい感じ
                continue
            end
        end
        if isnan(TFlist(i)) || TFlist(i) == 0
           continue
        end
        for j = 1:4
            % Imax_r(j,i) = max(separatrixList(i).max_r(j,1:2)); %日付でタイミング分ける
            % Imean_r(j,i) = max(separatrixList(i).mean_r(j,1:2));
            % Imax_l(j,i) = max(separatrixList(i).max_l(j,1:2)); 
            % Imean_l(j,i) = max(separatrixList(i).mean_l(j,1:2));
            % Imax_r(j,i) = max(separatrixList(i).max_r(j,:)); %一旦最大値
            % [Imean_r(j,i),Ir] = max(separatrixList(i).mean_r(j,:));
            % Imax_l(j,i) = max(separatrixList(i).max_l(j,:)); 
            % [Imean_l(j,i),Il] = max(separatrixList(i).mean_l(j,:));
            % Imax_r(j,i) = max(separatrixList(i).max_r(j,5)); %日付でタイミング分ける
            % [Imean_r(j,i),Ir] = max(separatrixList(i).mean_r(j,5));
            % Imax_l(j,i) = max(separatrixList(i).max_l(j,5)); 
            % [Imean_l(j,i),Il] = max(separatrixList(i).mean_l(j,5));  
            Imax_r(j,i) = max(separatrixList(i).max_r(j,5:6)); %日付でタイミング分ける
            [Imean_r(j,i),Ir] = max(separatrixList(i).mean_r(j,5:6));
            Imax_l(j,i) = max(separatrixList(i).max_l(j,5:6)); 
            [Imean_l(j,i),Il] = max(separatrixList(i).mean_l(j,5:6));
            Istd_r(j,i) = separatrixList(i).std_r(j,Ir);
            Istd_l(j,i) = separatrixList(i).std_l(j,Il);
        end
        TF(i) = BtList(i);
    end
end
Imax_r(Imax_r==0) = NaN;
Imean_r(Imean_r==0) = NaN;
Imax_l(Imax_l==0) = NaN;
Imean_l(Imean_l==0) = NaN;

if plotMax
    I_l = Imax_l;
    I_r = Imax_r;
    label_y_l = 'max intensity in lower left';
    label_y_r = 'max intensity in lower right';
else
    I_l = Imean_l;
    I_r = Imean_r;
    label_y_l = 'mean intensity in lower left';
    label_y_r = 'mean intensity in lower right';
end
I_l = I_l([1,2,4,3],:);
I_r = I_r([1,2,4,3],:);
Istd_l = Istd_l([1,2,4,3],:);
Istd_r = Istd_r([1,2,4,3],:);
Istd_l = I_l.*0.2;
Istd_r = I_r.*0.2;

switch plotType
    case 'shot'
        x = idxList_sxr;
        spec = '*';
        label_x = 'shot number';
    case 'Bt'
        x = TF;
        spec = '*';
        label_x = 'troidal magnetic field [T]';
    case 'Brec'
        x = BrList;
        spec = '-*';
        label_x = 'reconnection magnetic field [T]';
    case 'TF_group'
        x = TF;
        spec = '-*';
        label_x = 'troidal magnetic field [T]';
    case 'GFR'
        x = bList;
        spec = '*';
        label_x = 'guide field ratio';
    case 'GFR_group'
        x = GFR;
        spec = '-*';
        label_x = 'guide field ratio';
end
x(x==0) = NaN;
x = x.';

verificationMatrix = [I_l;I_r;Istd_l;Istd_r;x];
validColumns = ~any(isnan(verificationMatrix), 1);
I_l = I_l(:,validColumns);
I_r = I_r(:,validColumns);
Istd_l = Istd_l(:,validColumns);
Istd_r = Istd_r(:,validColumns);
x = x(:,validColumns);
[x,idx_x] = sort(x);
I_l = I_l(:,idx_x);
I_r = I_r(:,idx_x);
Istd_l = Istd_l(:,idx_x);
Istd_r = Istd_r(:,idx_x);

titlelist = {'~50eV','50~80eV','100eV~','200eV~'};
f_l = figure;%hold on;
f_l.Position = [440 278 560 600];
for i = 1:4
    subplot(4,1,i);
    if plotMax
        plot(x,I_l(i,:),spec);
    else
        errorbar(x,I_l(i,:),Istd_l(i,:),spec)
    end
    title(titlelist(i));
    ax = gca;
    ax.FontSize = 14;
    ylim([0 inf]);
end
xlabel(label_x);
sgtitle(label_y_l);
% ylabel(label_y_l);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');
f_r = figure;%hold on;
f_r.Position = [440 278 560 600];
for i = 1:4
    subplot(4,1,i);
    if plotMax
        plot(x,I_r(i,:),spec);
    else
        errorbar(x,I_r(i,:),Istd_r(i,:),spec)
    end
    title(titlelist(i));
    ax = gca;
    ax.FontSize = 14;
    ylim([0 inf]);
end
xlabel(label_x);
sgtitle(label_y_r);
% ylabel(label_y_r);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');
