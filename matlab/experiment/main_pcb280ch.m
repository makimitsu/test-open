%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ダイアログボックス入力で
% 磁気プローブによる磁気面をプロット
% カラープロット選択肢
% psi,Br,Bz,Bt,Bt_ext,Bt_plasma,Et,Jt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
clearvars -except PCB colorplot FIG doCheck
addpath '/Users/rsomeya/Documents/lab/matlab/common'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
run define_path.m

prompt = {'Date:','Shot Num.:','Colorplot Type:','Num. Plot (Vertical):','Num. Plot (Horizontal):','Start [us]:','dt [us]:','Check Signal (1):'};
dlgtitle = 'Input';
dims = [1 35];
if exist('PCB','var') && exist('colorplot','var') && exist('FIG','var') && exist('doCheck','var')
    definput = {num2str(PCB.date),num2str(PCB.IDX),colorplot,num2str(FIG.tate),num2str(FIG.yoko),num2str(FIG.start),num2str(FIG.dt),num2str(doCheck)};
else
    definput = {'','','psi','2','2','470','4','0'};
end
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
   % User clicked cancel. Bail out.
   return;
end
PCB.date = str2double(cell2mat(answer(1)));
PCB.IDX = str2double(cell2mat(answer(2)));
colorplot = cell2mat(answer(3));
FIG.tate = str2double(cell2mat(answer(4)));
FIG.yoko = str2double(cell2mat(answer(5)));
FIG.start = str2double(cell2mat(answer(6)));
FIG.dt = str2double(cell2mat(answer(7)));
doCheck = str2double(cell2mat(answer(8)));

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
IDSPminZ=T.IDSPZ_cm_(PCB.IDX);
IDSPminR=T.IDSPMinR_cm_(PCB.IDX);

PCB.trange=400:800;%【input】計算時間範囲
PCB.mesh=40; %【input】rz方向のメッシュ数

IDSP.z = IDSPminZ*ones(7,1)*1E-2;
IDSP.r = (IDSPminR:2.5:IDSPminR+6*2.5)*1E-2;

for i=1:n_data
    PCB.shot=shotlist(i,:);
    PCB.tfshot=tfshotlist(i,:);
    if PCB.tfshot == PCB.shot
        PCB.tfshot = [0,0];
    end
    PCB.i_EF=EFlist(i);
    if doCheck == 1
        check_signal(PCB,pathname);
    else
        [PCBgrid2D,PCBdata2D] = cal_psi(PCB,pathname);
        plot_psi(PCBgrid2D,PCBdata2D,IDSP,FIG,colorplot);
    end
end
