clearvars -except date doFilter doNLR
addpath '/Users/shinjirotakeda/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
addpath '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/Four-View_Simulation';

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

sxrDataDir = getenv('SXRDATA_DIR');
sxrDataFile = strcat(sxrDataDir,filesep,num2str(date),'_',options,'.mat');
if exist(sxrDataFile,'file')
    load(sxrDataFile,'idxList','xpointList','downstreamList');
    idxList_sxr = idxList;
else
    disp('No sxr data');
    return
end

magDataDir = getenv('MAGDATA_DIR');
magDataFile = strcat(magDataDir,'/',num2str(date),'.mat');
if exist(magDataFile,'file')
    load(magDataFile,'idxList','BrList','BtList','bList');
    idxList_mag = idxList;
else
    disp('No mag data');
    return
end

tablePath = '/Users/shinjirotakeda/Library/CloudStorage/OneDrive-TheUniversityofTokyo/新規論文_2023/230920_separatrix.xlsx';
T2 = readtable(tablePath);

I_r = zeros(numel(T2.shot),4);
I_l = I_r;
for i = 1:numel(T2.shot)
    I_r(i,1) = T2.Al_10_r(i);
    I_l(i,1) = T2.Al_10_l(i);
    I_r(i,2) = T2.Al_25_r(i);
    I_l(i,2) = T2.Al_25_l(i);
    I_r(i,3) = T2.My_10_r(i);
    I_l(i,3) = T2.My_10_l(i);
    I_r(i,4) = T2.My_20_r(i);
    I_l(i,4) = T2.My_20_l(i);
end

I_l(I_l==0) = NaN;
I_r(I_r==0) = NaN;

f_r = figure;
% hold on;plot(BrList(T2.shot),I_r,'o');
for i = 1:4
    subplot(4,1,i);plot(BrList(T2.shot),I_r(:,i),'o');
end
f_l = figure;
% hold on;plot(BrList(T2.shot),I_l,'o');
for i = 1:4
    subplot(4,1,i);plot(BrList(T2.shot),I_l(:,i),'o');
end