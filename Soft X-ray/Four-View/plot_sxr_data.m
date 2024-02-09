% clear
% close all
clearvars -except date doFilter doNLR
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
prompt = {'Date:','doFilter:','doNLR:'};
definput = {'','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
if exist('doFilter','var')
    definput{2} = num2str(doFilter);
end
if exist('doNLR','var')
    definput{3} = num2str(doNLR);
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
    load(sxrDataFile,'idxList','xpointList');
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

Imax = zeros(1,n_data);
Imean = zeros(1,n_data);
TF = zeros(1,n_data);
for i = 1:n_data
    Imax(i) = xpointList(i).max(2,2);
    Imean(i) = xpointList(i).mean(2,2);
    if xpointList(i).MR(2) > 0.5
        Imax(i) = xpointList(i).max(2,1);
        Imean(i) = xpointList(i).mean(2,1);
    end
    if ismember(idxList(i),[6,11])
        Imax(i) = NaN;
        Imean(i) = NaN;
    end
    TF(i) = BtList(i);
end
TF(TF==0) = NaN;
Imax(Imax==0) = NaN;
Imean(Imean==0) = NaN;
% figure;
% subplot(1,2,1);plot(TF,Imax,'*');
% subplot(1,2,2);plot(TF,Imean,'*');
% figure;
% subplot(1,2,1);plot(idxList,Imax,'*');
% subplot(1,2,2);plot(idxList,Imean,'*');

figure;plot(TF,Imax,'*');
xlabel('troidal field');ylabel('intensity');

% MR = zeros(1,n_data);
% TF = zeros(1,n_data);
% for i = 1:n_data
%     MR(i) = xpointList(i).MR(2);
%     if MR(i) > 0.5
%         MR(i) = xpointList(i).MR(1);
%     end
%     TF(i) = BtList(i);
% end
% TF(TF==0) = NaN;
% MR(MR==0) = NaN;
% figure;plot(TF,MR,'*');
% figure;plot(idxList,MR,'*');


% MR40 = zeros(8,numel(find(TFlist==4)));
% t40 = zeros(8,numel(find(TFlist==4)));
% MR35 = zeros(8,numel(find(TFlist==3.5)));
% t35 = zeros(8,numel(find(TFlist==3.5)));
% MR30 = zeros(8,numel(find(TFlist==3)));
% t30 = zeros(8,numel(find(TFlist==3)));
% MR25 = zeros(8,numel(find(TFlist==2.5)));
% t25 = zeros(8,numel(find(TFlist==2.5)));
% i40 = 1;
% i35 = 1;
% i30 = 1;
% i25 = 1;
% 
% for i = 1:numel(idxList_sxr)
%     if TFlist(i) == 4
%         MR40(:,i40) = (xpointList(i).MR).';
%         t40(:,i40) = (xpointList(i).t).';
%         i40 = i40+1;
%     end
%     if TFlist(i) == 3.5
%         MR35(:,i35) = (xpointList(i).MR).';
%         t35(:,i35) = (xpointList(i).t).';
%         i35 = i35+1;
%     end
%     if TFlist(i) == 3
%         MR30(:,i30) = (xpointList(i).MR).';
%         t30(:,i30) = (xpointList(i).t).';
%         i30 = i30+1;
%     end
%     if TFlist(i) == 2.5
%         MR25(:,i25) = (xpointList(i).MR).';
%         t25(:,i25) = (xpointList(i).t).';
%         i25 = i25+1;
%     end
% end

% figure;
% subplot(2,2,1);
% plot(t40,MR40,'*');
% xlim([460 480]);
% ylim([0 1]);
% subplot(2,2,2);
% plot(t35,MR35,'*');
% xlim([460 480]);
% ylim([0 1]);
% subplot(2,2,3);
% plot(t30,MR30,'*');
% xlim([460 480]);
% ylim([0 1]);
% subplot(2,2,4);
% plot(t25,MR25,'*');
% xlim([460 480]);
% ylim([0 1]);

% Imax = NaN(numel(idxList_sxr),4);
% Imean = NaN(numel(idxList_sxr),4);
% for i = 1:n_data
%     for j = 1:4
%         Imax(i,j) = max(xpointList(i).max(j,:));
%         Imean(i,j) = max(xpointList(i).mean(j,:));
%     end
% end
% 
% figure;
% for i = 1:4
%     subplot(2,2,i);
%     plot(BtList,Imax(:,i),'*');
%     xlabel('guide field');
% end
% 
% figure;
% for i = 1:4
%     subplot(2,2,i);
%     plot(BtList,Imean(:,i),'*');
%     xlabel('guide field');
% end

% fmax = figure;
% 
% fmean = figure;
% 
% fMR = figure;
% 
% for i=1:n_data
%     figure(fmax);
%     for j = 1:4
%         subplot(2,2,j);
%         hold on
%         plot(xpointList(i).t,xpointList(i).max(j,:));
%     end
%     figure(fmean);
%     for j = 1:4
%         subplot(2,2,j);
%         hold on
%         plot(xpointList(i).t,xpointList(i).mean(j,:));
%     end
%     figure(fMR);
%     hold on
%     plot(xpointList(i).t,xpointList(i).MR);
% end
% 
% figure(fmax);
% for i = 1:4
%     subplot(2,2,i)
%     xlabel('time [us]');
%     ylabel('max intensity [a.u.]');
%     xlim([450 480]);
% end
% 
% figure(fmean);
% for i = 1:4
%     subplot(2,2,i)
%     xlabel('time [us]');
%     ylabel('mean intensity [a.u.]');
%     xlim([450 480]);
% end
% 
% figure(fMR);
% xlabel('time [us]');
% ylabel('merging ratio');
% xlim([450 480]);