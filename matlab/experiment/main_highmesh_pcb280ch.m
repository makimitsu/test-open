%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ダイアログボックス入力で
% 磁気プローブによる磁気面をプロット
% カラープロット選択肢
% psi,Br,Bz,Bt,Bt_ext,Bt_plasma,Et,Jt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
clearvars -except PCB
addpath '/Users/rsomeya/Documents/lab/matlab/common'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
run define_path.m

prompt = {'Date:','Shot Num.:','Num. Mesh:','Start cal. Time:','End cal. Time:'};
dlgtitle = 'Input';
dims = [1 35];
if exist('PCB','var')
    definput = {num2str(PCB.date),num2str(PCB.IDX),num2str(PCB.mesh),num2str(PCB.trange(1)),num2str(PCB.trange(end))};
else
    definput = {'','','500','460','490'};
end
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
   % User clicked cancel. Bail out.
   return;
end
PCB.date = str2double(cell2mat(answer(1)));
PCB.IDX = str2double(cell2mat(answer(2)));
PCB.mesh = str2double(cell2mat(answer(3)));
PCB.trange=str2double(cell2mat(answer(4))):str2double(cell2mat(answer(5)));%計算時間範囲

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

for i=1:n_data
    PCB.shot=shotlist(i,:);
    PCB.tfshot=tfshotlist(i,:);
    if PCB.tfshot == PCB.shot
        PCB.tfshot = [0,0];
    end
    PCB.i_EF=EFlist(i);
    [highmesh_PCBgrid2D,highmesh_PCBdata2D] = cal_highmesh_psi(PCB,pathname);
end
