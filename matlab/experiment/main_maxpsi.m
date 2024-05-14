%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 磁気プローブによるpsi最大値の時間変化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
clearvars -except PCB
addpath '/Users/rsomeya/Documents/lab/matlab/common'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
run define_path.m

prompt = {'Date:','Shot Num.:','Start [us]:','End [us]:'};
dlgtitle = 'Input';
dims = [1 35];
if exist('PCB','var') && exist('FIG','var')
    definput = {num2str(PCB.date),num2str(PCB.IDX),num2str(FIG.start),num2str(FIG.end)};
else
    definput = {'230830','22','450','490'};
end
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
   % User clicked cancel. Bail out.
   return;
end
PCB.date = str2double(cell2mat(answer(1)));
PCB.IDX = str2double(cell2mat(answer(2)));
FIG.start = str2double(cell2mat(answer(3)));
FIG.end = str2double(cell2mat(answer(4)));

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,PCB.date);
n_data=numel(PCB.IDX);%計測データ数
shotlist_a039 =T.a039(PCB.IDX);
shotlist_a040 = T.a040(PCB.IDX);
shotlist = [shotlist_a039, shotlist_a040];
tfshotlist_a039 =T.a039_TF(PCB.IDX);
tfshotlist_a040 =T.a040_TF(PCB.IDX);
tfshotlist = [tfshotlist_a039, tfshotlist_a040];
EFlist=T.EF_A_(PCB.IDX);
TFlist=T.TF_kV_(PCB.IDX);

PCB.trange=400:800;%【input】計算時間範囲
PCB.mesh=40; %【input】rz方向のメッシュ数

for i=1:n_data
    PCB.shot=shotlist(i,:);
    PCB.tfshot=tfshotlist(i,:);
    if PCB.tfshot == PCB.shot
        PCB.tfshot = [0,0];
    end
    PCB.i_EF=EFlist(i);
    [PCBgrid2D,PCBdata2D] = cal_psi(PCB,pathname);
    plot_maxpsi(PCBgrid2D,PCBdata2D,FIG);
end
