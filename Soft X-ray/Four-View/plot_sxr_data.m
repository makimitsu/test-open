% clear
% close all
clearvars -except date doFilter doNLR plotArea plotType
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
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    return
end

date = str2double(cell2mat(answer(1)));
doFilter = logical(str2num(cell2mat(answer(2))));
doNLR = logical(str2num(cell2mat(answer(3))));
plotArea = answer(4);
plotType = answer(5);

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

% 以下セパラトリクス

% shotごとにプロット（最大値/平均値）
Imax_r = zeros(4,n_data);
Imean_r = zeros(4,n_data);
Imax_l = zeros(4,n_data);
Imean_l = zeros(4,n_data);
TF = zeros(1,n_data);
for i = 1:n_data
    for j = 1:4
        Imax_r(j,i) = max(separatrixList(i).max_r(j,1:2)); 
        Imean_r(j,i) = max(separatrixList(i).mean_r(j,1:2));
        Imax_l(j,i) = max(separatrixList(i).max_l(j,1:2)); 
        Imean_l(j,i) = max(separatrixList(i).mean_l(j,1:2));
    end
    TF(i) = BtList(i);
end
TF(TF==0) = NaN;
Imax_r(Imax_r==0) = NaN;
Imean_r(Imean_r==0) = NaN;
Imax_l(Imax_l==0) = NaN;
Imean_l(Imean_l==0) = NaN;
% for i = 1:4 %規格化
%     Imax_r(i,:) = Imax_r(i,:)./max(Imax_r(i,:));
%     Imax_l(i,:) = Imax_l(i,:)./max(Imax_l(i,:));
%     Imean_r(i,:) = Imean_r(i,:)./max(Imean_r(i,:));
%     Imean_l(i,:) = Imean_l(i,:)./max(Imean_l(i,:));
% end

% figure;hold on;%shot, max
% for i = 1:4
%     plot(idxList_sxr,Imax_l(i,:),'*');
% end
% xlabel('shot number');ylabel('max intensity in lower left');
% figure;hold on;
% for i = 1:4
%     plot(idxList_sxr,Imax_r(i,:),'*');
% end
% xlabel('shot number');ylabel('max intensity in lower right');

figure;hold on;%shot, mean
for i = 1:4
    plot(idxList_sxr,Imax_l(i,:),'*');
end
xlabel('shot number');ylabel('mean intensity in lower left');
figure;hold on;
for i = 1:4
    plot(idxList_sxr,Imax_r(i,:),'*');
end
xlabel('shot number');ylabel('mean intensity in lower right');

% figure;hold on;%TF,max
% for i = 1:4
%     plot(BtList,Imax_l(i,:),'*');
% end
% xlabel('troidal field');ylabel('max intensity in lower left');xlim([0.08 0.16]);
% figure;hold on;
% for i = 1:4
%     plot(BtList,Imax_r(i,:),'*');
% end
% xlabel('troidal field');ylabel('max intensity in lower right');xlim([0.08 0.16]);

% figure;hold on;%TF, mean
% for i = 1:4
%     plot(BtList,Imean_l(i,:),'*');
% end
% xlabel('troidal field');ylabel('mean intensity in lower left');xlim([0.08 0.16]);
% figure;hold on;
% for i = 1:4
%     plot(BtList,Imean_r(i,:),'*');
% end
% xlabel('troidal field');ylabel('mean intensity in lower right');xlim([0.08 0.16]);

% % GFRごとにプロット,最初の2タイミングの最大値
% Imax_r = zeros(4,4);
% Imean_r = zeros(4,4);
% Imax_l = zeros(4,4);
% Imean_l = zeros(4,4);
% GFR = zeros(1,4);
% cnt = zeros(1,4);
% for i = 1:n_data
%     if isnan(TFlist(i))
%        continue 
%     end
%     k = find([2.5,3,3.5,4]==TFlist(i));
%     for j = 1:4
%         Imax_r_tmp = max(separatrixList(i).max_r(j,1:2));
%         Imean_r_tmp = max(separatrixList(i).mean_r(j,1:2));
%         Imax_l_tmp = max(separatrixList(i).max_l(j,1:2));
%         Imean_l_tmp = max(separatrixList(i).mean_l(j,1:2));
%         Imax_r(j,k) = (Imax_r(j,k)*cnt(k)+Imax_r_tmp)/(cnt(k)+1);
%         Imean_r(j,k) = (Imean_r(j,k)*cnt(k)+Imean_r_tmp)/(cnt(k)+1);
%         Imax_l(j,k) = (Imax_l(j,k)*cnt(k)+Imax_l_tmp)/(cnt(k)+1);
%         Imean_l(j,k) = (Imean_l(j,k)*cnt(k)+Imean_l_tmp)/(cnt(k)+1);
%     end
%     if BtList(i)~=0
%         GFR(k) = (GFR(k)*cnt(k)+bList(i))/(cnt(k)+1);
%     end
%     cnt(k) = cnt(k)+1;
% end
% GFR(GFR==0) = NaN;
% Imax_r(Imax_r==0) = NaN;
% Imean_r(Imean_r==0) = NaN;
% Imax_l(Imax_l==0) = NaN;
% Imean_l(Imean_l==0) = NaN;
% % for i = 1:4 %規格化
% %     Imax_r(i,:) = Imax_r(i,:)./max(Imax_r(i,:));
% %     Imax_l(i,:) = Imax_l(i,:)./max(Imax_l(i,:));
% %     Imean_r(i,:) = Imean_r(i,:)./max(Imean_r(i,:));
% %     Imean_l(i,:) = Imean_l(i,:)./max(Imean_l(i,:));
% % end
% % [GFR,I] = sort(GFR);
% % Imax_r = Imax_r(:,I);
% % Imean_r = Imean_r(:,I);
% % Imax_l = Imax_l(:,I);
% % Imean_l = Imean_l(:,I);
% 
% figure;hold on;
% for i = 1:4
%     plot(GFR,Imax_l(i,:),'-*');
% end
% xlabel('guide field ratio');ylabel('max intensity in lower left');%xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');
% figure;hold on;
% for i = 1:4
%     plot(GFR,Imax_r(i,:),'-*');
% end
% xlabel('guide field ratio');ylabel('max intensity in lower right');%xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');

% % GFRごとにプロット,最初の2タイミングの平均
% Imax_r = zeros(4,4);
% Imean_r = zeros(4,4);
% Imax_l = zeros(4,4);
% Imean_l = zeros(4,4);
% GFR = zeros(1,4);
% cnt = zeros(1,4);
% for i = 1:n_data
%     if isnan(TFlist(i))
%        continue 
%     end
%     k = find([2.5,3,3.5,4]==TFlist(i));
%     for j = 1:4
%         Imax_r_tmp = mean(separatrixList(i).max_r(j,1:2));
%         Imean_r_tmp = mean(separatrixList(i).mean_r(j,1:2));
%         Imax_l_tmp = mean(separatrixList(i).max_l(j,1:2));
%         Imean_l_tmp = mean(separatrixList(i).mean_l(j,1:2));
%         Imax_r(j,k) = (Imax_r(j,k)*cnt(k)+Imax_r_tmp)/(cnt(k)+1);
%         Imean_r(j,k) = (Imean_r(j,k)*cnt(k)+Imean_r_tmp)/(cnt(k)+1);
%         Imax_l(j,k) = (Imax_l(j,k)*cnt(k)+Imax_l_tmp)/(cnt(k)+1);
%         Imean_l(j,k) = (Imean_l(j,k)*cnt(k)+Imean_l_tmp)/(cnt(k)+1);
%     end
%     if BtList(i)~=0
%         GFR(k) = (GFR(k)*cnt(k)+bList(i))/(cnt(k)+1);
%     end
%     cnt(k) = cnt(k)+1;
% end
% GFR(GFR==0) = NaN;
% Imax_r(Imax_r==0) = NaN;
% Imean_r(Imean_r==0) = NaN;
% Imax_l(Imax_l==0) = NaN;
% Imean_l(Imean_l==0) = NaN;
% % for i = 1:4 %規格化
% %     Imax_r(i,:) = Imax_r(i,:)./max(Imax_r(i,:));
% %     Imax_l(i,:) = Imax_l(i,:)./max(Imax_l(i,:));
% %     Imean_r(i,:) = Imean_r(i,:)./max(Imean_r(i,:));
% %     Imean_l(i,:) = Imean_l(i,:)./max(Imean_l(i,:));
% % end
% % [GFR,I] = sort(GFR);
% % Imax_r = Imax_r(:,I);
% % Imean_r = Imean_r(:,I);
% % Imax_l = Imax_l(:,I);
% % Imean_l = Imean_l(:,I);
% 
% figure;hold on;
% for i = 1:4
%     plot(GFR,Imax_l(i,:),'-*');
% end
% xlabel('guide field ratio');ylabel('max intensity in lower left');%xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');
% figure;hold on;
% for i = 1:4
%     plot(GFR,Imax_r(i,:),'-*');
% end
% xlabel('guide field ratio');ylabel('max intensity in lower right');%xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');

% % GFRごとにプロット（最大値/平均値）
% Imax_r = zeros(4,4);
% Imean_r = zeros(4,4);
% Imax_l = zeros(4,4);
% Imean_l = zeros(4,4);
% GFR = zeros(1,4);
% cnt = zeros(1,4);
% for i = 1:n_data
%     if isnan(TFlist(i))
%        continue 
%     end
%     k = find([2.5,3,3.5,4]==TFlist(i));
%     for j = 1:4
%         Imax_r(j,k) = (Imax_r(j,k)*cnt(k)+separatrixList(i).max_r(j,1))/(cnt(k)+1);
%         Imean_r(j,k) = (Imean_r(j,k)*cnt(k)+separatrixList(i).mean_r(j,1))/(cnt(k)+1);
%         Imax_l(j,k) = (Imax_l(j,k)*cnt(k)+separatrixList(i).max_l(j,1))/(cnt(k)+1);
%         Imean_l(j,k) = (Imean_l(j,k)*cnt(k)+separatrixList(i).mean_l(j,1))/(cnt(k)+1);
%     end
%     if BtList(i)~=0
%         GFR(k) = (GFR(k)*cnt(k)+bList(i))/(cnt(k)+1);
%     end
%     cnt(k) = cnt(k)+1;
% end
% GFR(GFR==0) = NaN;
% Imax_r(Imax_r==0) = NaN;
% Imean_r(Imean_r==0) = NaN;
% Imax_l(Imax_l==0) = NaN;
% Imean_l(Imean_l==0) = NaN;
% % for i = 1:4 %規格化
% %     Imax_r(i,:) = Imax_r(i,:)./max(Imax_r(i,:));
% %     Imax_l(i,:) = Imax_l(i,:)./max(Imax_l(i,:));
% %     Imean_r(i,:) = Imean_r(i,:)./max(Imean_r(i,:));
% %     Imean_l(i,:) = Imean_l(i,:)./max(Imean_l(i,:));
% % end
% % [GFR,I] = sort(GFR);
% % Imax_r = Imax_r(:,I);
% % Imean_r = Imean_r(:,I);
% % Imax_l = Imax_l(:,I);
% % Imean_l = Imean_l(:,I);
% 
% figure;hold on;
% for i = 1:4
%     plot(GFR,Imax_l(i,:),'-*');
% end
% xlabel('guide field ratio');ylabel('max intensity in lower left');%xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');
% figure;hold on;
% for i = 1:4
%     plot(GFR,Imax_r(i,:),'-*');
% end
% xlabel('guide field ratio');ylabel('max intensity in lower right');%xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');

% % TFごとにプロット（最大値/平均値）
% Imax_r = zeros(4,4);
% Imean_r = zeros(4,4);
% Imax_l = zeros(4,4);
% Imean_l = zeros(4,4);
% TF = zeros(1,4);
% cnt = zeros(1,4);
% for i = 1:n_data
%     if isnan(TFlist(i))
%        continue 
%     end
%     k = find([2.5,3,3.5,4]==TFlist(i));
%     for j = 1:4
%         Imax_r(j,k) = (Imax_r(j,k)*cnt(k)+separatrixList(i).max_r(j,1))/(cnt(k)+1);
%         Imean_r(j,k) = (Imean_r(j,k)*cnt(k)+separatrixList(i).mean_r(j,1))/(cnt(k)+1);
%         Imax_l(j,k) = (Imax_l(j,k)*cnt(k)+separatrixList(i).max_l(j,1))/(cnt(k)+1);
%         Imean_l(j,k) = (Imean_l(j,k)*cnt(k)+separatrixList(i).mean_l(j,1))/(cnt(k)+1);
%     end
%     if BtList(i)~=0
%         TF(k) = (TF(k)*cnt(k)+BtList(i))/(cnt(k)+1);
%     end
%     cnt(k) = cnt(k)+1;
% end
% TF(TF==0) = NaN;
% Imax_r(Imax_r==0) = NaN;
% Imean_r(Imean_r==0) = NaN;
% Imax_l(Imax_l==0) = NaN;
% Imean_l(Imean_l==0) = NaN;
% for i = 1:4 %規格化
%     Imax_r(i,:) = Imax_r(i,:)./max(Imax_r(i,:));
%     Imax_l(i,:) = Imax_l(i,:)./max(Imax_l(i,:));
%     Imean_r(i,:) = Imean_r(i,:)./max(Imean_r(i,:));
%     Imean_l(i,:) = Imean_l(i,:)./max(Imean_l(i,:));
% end
% 
% figure;hold on;
% for i = 1:4
%     plot(TF,Imax_l(i,:),'-*');
% end
% xlabel('troidal field');ylabel('max intensity in lower left');xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');
% figure;hold on;
% for i = 1:4
%     plot(TF,Imax_r(i,:),'-*');
% end
% xlabel('troidal field');ylabel('max intensity in lower right');xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');

% 以下下流領域

% % TFごとにプロット（最大値/平均値）
% Imax_r = zeros(4,4);
% Imean_r = zeros(4,4);
% Imax_l = zeros(4,4);
% Imean_l = zeros(4,4);
% TF = zeros(1,4);
% cnt = zeros(1,4);
% for i = 1:n_data
%     if isnan(TFlist(i))
%        continue 
%     end
%     k = find([2.5,3,3.5,4]==TFlist(i));
%     for j = 1:4
%         Imax_r(j,k) = (Imax_r(j,k)*cnt(k)+downstreamList(i).max_r(j,1))/(cnt(k)+1);
%         Imean_r(j,k) = (Imean_r(j,k)*cnt(k)+downstreamList(i).mean_r(j,1))/(cnt(k)+1);
%         Imax_l(j,k) = (Imax_l(j,k)*cnt(k)+downstreamList(i).max_l(j,1))/(cnt(k)+1);
%         Imean_l(j,k) = (Imean_l(j,k)*cnt(k)+downstreamList(i).mean_l(j,1))/(cnt(k)+1);
%     end
%     if BtList(i)~=0
%         TF(k) = (TF(k)*cnt(k)+BtList(i))/(cnt(k)+1);
%     end
%     cnt(k) = cnt(k)+1;
% end
% TF(TF==0) = NaN;
% Imax_r(Imax_r==0) = NaN;
% Imean_r(Imean_r==0) = NaN;
% Imax_l(Imax_l==0) = NaN;
% Imean_l(Imean_l==0) = NaN;
% % for i = 1:4 %規格化
% %     Imax_r(i,:) = Imax_r(i,:)./max(Imax_r(i,:));
% %     Imax_l(i,:) = Imax_l(i,:)./max(Imax_l(i,:));
% %     Imean_r(i,:) = Imean_r(i,:)./max(Imean_r(i,:));
% %     Imean_l(i,:) = Imean_l(i,:)./max(Imean_l(i,:));
% % end
% 
% % figure;hold on;
% % for i = 1:4
% %     plot(idxList_sxr,Imax_l(i,:),'*');
% % end
% % xlabel('shot number');ylabel('max intensity in lower left');
% % figure;hold on;
% % for i = 1:4
% %     plot(idxList_sxr,Imax_r(i,:),'*');
% % end
% % xlabel('shot number');ylabel('max intensity in lower right');
% 
% figure;hold on;
% for i = 1:4
%     plot(TF,Imax_l(i,:),'*');
% end
% xlabel('troidal field');ylabel('max intensity in lower left');xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');
% figure;hold on;
% for i = 1:4
%     plot(TF,Imax_r(i,:),'*');
% end
% xlabel('troidal field');ylabel('max intensity in lower right');xlim([0.08 0.16]);
% legend({'~50eV','50~80eV','200eV~','100eV~'},'Location','south');

% % shotごとにプロット（最大値/平均値）
% Imax_r = zeros(4,n_data);
% Imean_r = zeros(4,n_data);
% Imax_l = zeros(4,n_data);
% Imean_l = zeros(4,n_data);
% TF = zeros(1,n_data);
% for i = 1:n_data
%     for j = 1:4
%         Imax_r(j,i) = downstreamList(i).max_r(j,1); 
%         Imean_r(j,i) = downstreamList(i).mean_r(j,1);
%         Imax_l(j,i) = downstreamList(i).max_l(j,1); 
%         Imean_l(j,i) = downstreamList(i).mean_l(j,1);
%     end
%     TF(i) = BtList(i);
% end
% TF(TF==0) = NaN;
% Imax_r(Imax_r==0) = NaN;
% Imean_r(Imean_r==0) = NaN;
% Imax_l(Imax_l==0) = NaN;
% Imean_l(Imean_l==0) = NaN;
% for i = 1:4 %規格化
%     Imax_r(i,:) = Imax_r(i,:)./max(Imax_r(i,:));
%     Imax_l(i,:) = Imax_l(i,:)./max(Imax_l(i,:));
%     Imean_r(i,:) = Imean_r(i,:)./max(Imean_r(i,:));
%     Imean_l(i,:) = Imean_l(i,:)./max(Imean_l(i,:));
% end
% 
% % figure;hold on;
% % for i = 1:4
% %     plot(idxList_sxr,Imax_l(i,:),'*');
% % end
% % xlabel('shot number');ylabel('max intensity in lower left');
% % figure;hold on;
% % for i = 1:4
% %     plot(idxList_sxr,Imax_r(i,:),'*');
% % end
% % xlabel('shot number');ylabel('max intensity in lower right');
% 
% figure;hold on;
% for i = 1:4
%     plot(BtList,Imax_l(i,:),'*');
% end
% xlabel('troidal field');ylabel('max intensity in lower left');xlim([0.08 0.16]);
% figure;hold on;
% for i = 1:4
%     plot(BtList,Imax_r(i,:),'*');
% end
% xlabel('troidal field');ylabel('max intensity in lower right');xlim([0.08 0.16]);

% 以下X点

% % TFごとにグループ化してプロット（最大値）
% Imax = zeros(1,4);
% Imean = zeros(1,4);
% TF = zeros(1,4);
% cnt_i = [0,0,0,0];
% for j = 1:n_data
%     if isnan(TFlist(j))
%        continue 
%     end
%     i = find([2.5,3,3.5,4]==TFlist(j));
%     Imax(i) = (Imax(i)*cnt_i(i) + max(xpointList(j).max(2,[1,2])))/(cnt_i(i)+1); %こっちだと右肩下がり
%     Imean(i) = (Imean(i)*cnt_i(i) + max(xpointList(j).mean(2,[1,2])))/(cnt_i(i)+1);
%     if BtList(j)~=0
%         TF(i) = (TF(i)*cnt_i(i)+BtList(j))/(cnt_i(i)+1);
%     end
%     cnt_i(i) = cnt_i(i) + 1;
% end
% TF(TF==0) = NaN;
% Imax(Imax==0) = NaN;
% Imean(Imean==0) = NaN;
% 
% figure;plot(TF,Imax,'*');
% xlabel('troidal field');ylabel('max intensity');
% 
% figure;plot(TF,Imean,'*');
% xlabel('troidal field');ylabel('mean intensity');

% % TFごとにグループ化してプロット（最初の2枚の平均）
% Imax = zeros(1,4);
% Imean = zeros(1,4);
% TF = zeros(1,4);
% cnt_i = [0,0,0,0];
% for j = 1:n_data
%     if isnan(TFlist(j))
%        continue 
%     end
%     i = find([2.5,3,3.5,4]==TFlist(j));
%     Imax(i) = (Imax(i)*cnt_i(i) + mean(xpointList(j).max(2,[1,2])))/(cnt_i(i)+1); %こっちだと右肩下がり
%     Imean(i) = (Imean(i)*cnt_i(i) + mean(xpointList(j).mean(2,[1,2])))/(cnt_i(i)+1);
%     if BtList(j)~=0
%         TF(i) = (TF(i)*cnt_i(i)+BtList(j))/(cnt_i(i)+1);
%     end
%     cnt_i(i) = cnt_i(i) + 1;
% end
% TF(TF==0) = NaN;
% Imax(Imax==0) = NaN;
% Imean(Imean==0) = NaN;
% 
% figure;plot(TF,Imax,'*');
% xlabel('troidal field');ylabel('max intensity');
% 
% figure;plot(TF,Imean,'*');
% xlabel('troidal field');ylabel('mean intensity');

% % shotごとにプロット（最大値/平均値）
% Imax = zeros(1,n_data);
% Imean = zeros(1,n_data);
% TF = zeros(1,n_data);
% for i = 1:n_data
%     % Imax(i) = mean(xpointList(i).max(2,[1,2])); %こっちだと右肩下がり
%     % Imean(i) = mean(xpointList(i).mean(2,[1,2]));
%     Imax(i) = max(xpointList(i).max(2,[1,2])); %こっちだと右肩下がり
%     Imean(i) = max(xpointList(i).mean(2,[1,2]));
%     % Imax(i) = xpointList(i).max(2,2); %こっちだと右肩上がり
%     % Imean(i) = xpointList(i).mean(2,2);
%     % if xpointList(i).MR(2) > 0.5
%     %     Imax(i) = xpointList(i).max(2,1);
%     %     Imean(i) = xpointList(i).mean(2,1);
%     % end
%     % if ismember(idxList(i),[6,11])
%     %     Imax(i) = NaN;
%     %     Imean(i) = NaN;
%     % end
%     TF(i) = BtList(i);
% end
% TF(TF==0) = NaN;
% Imax(Imax==0) = NaN;
% Imean(Imean==0) = NaN;
% 
% % figure;plot(TF,Imax,'*');
% % xlabel('troidal field');ylabel('max intensity');
% % figure;plot(TF,Imean,'*');
% % xlabel('troidal field');ylabel('mean intensity');
% 
% figure;plot(idxList_sxr,Imax,'*');
% xlabel('shot number');ylabel('max intensity');
% figure;plot(idxList_sxr,Imean,'*');
% xlabel('shot number');ylabel('mean intensity');

% % shotごとに合体率をプロット
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

% % TFごとに違う図で時間発展をプロット
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

% % ショットごとに発光強度（全時間での最大値）をプロット
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

% % 各ショットごとに発光強度（最大・平均）と合体率の時間発展をプロット
% fmax = figure;
% fmean = figure;
% fMR = figure;
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
% figure(fmean);
% for i = 1:4
%     subplot(2,2,i)
%     xlabel('time [us]');
%     ylabel('mean intensity [a.u.]');
%     xlim([450 480]);
% end
% figure(fMR);
% xlabel('time [us]');
% ylabel('merging ratio');
% xlim([450 480]);