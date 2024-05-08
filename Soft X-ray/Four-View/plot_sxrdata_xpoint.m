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
Imax = zeros(4,num_data);
Imean = zeros(4,num_data);
Istd = zeros(4,num_data);

if strcmp(plotType,'TF_group')||strcmp(plotType,'GFR_group')
    cnt = zeros(1,4);
    for i = 1:n_data
        if isnan(TFlist(i)) || TFlist(i) == 0
           continue 
        end
        k = find([2.5,3,3.5,4]==TFlist(i));
        for j = 1:4
            Imax_tmp = max(xpointList(i).max(j,1:2));
            Imean_tmp = max(xpointList(i).mean(j,1:2));
            Imax(j,k) = (Imax(j,k)*cnt(k)+Imax_tmp)/(cnt(k)+1);
            Imean(j,k) = (Imean(j,k)*cnt(k)+Imean_tmp)/(cnt(k)+1);
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
            % if ~ismember(i,[38,39,56:57])
            if i>=54
                continue
            end
        elseif date == 240111
            if ~ismember(i,[7,11,15,19])
                continue
            end
        end
        if isnan(TFlist(i)) || TFlist(i) == 0
           continue
        end
        for j = 1:4
            Imax(j,i) = max(xpointList(i).max(j,:)); %一旦最大値
            [Imean(j,i),I] = max(xpointList(i).mean(j,:));
            Istd(j,i) = xpointList(i).std(j,I);
        end
        TF(i) = BtList(i);
    end
end
Imax(Imax==0) = NaN;
Imean(Imean==0) = NaN;
if plotMax
    I_plot = Imax;
    label_y = 'max intensity at x-point';
else
    I_plot = Imean;
    label_y = 'mean intensity at x-point';
end
I_plot = I_plot([1,2,4,3],:);
% Istd = Istd([1,2,4,3],:);
Istd = I_plot.*0.2;
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
        spec = '*';
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
titlelist = {'~50eV','50~80eV','100eV~','200eV~'};
figure;%hold on;
for i = 1:4
    subplot(4,1,i);
    if plotMax
        plot(x,I_plot(i,:),spec);
    else
        errorbar(x,I_plot(i,:),Istd(i,:),spec)
    end
    title(titlelist(i));
    ylim([0 inf]);
end
xlabel(label_x);
sgtitle(label_y);
% ylabel(label_y);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');
