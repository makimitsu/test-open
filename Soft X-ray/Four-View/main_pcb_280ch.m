clear
addpath '/Users/yuleo/Documents/GitHub/test-open/pcb_experiment'; %getMDSdata.mとcoeff200ch.xlsxのあるフォルダへのパス
addpath 'C:\Users\yuleo\Documents\GitHub\test-open'

%%%%%%%%%%%%%%%%%%%%%%%%
%280ch用新規pcbプローブのみでの磁気面（Bz）
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.pre_processed_directory = getenv('pre_processed_directory_path');%計算結果の保存先（どこでもいい）

%%%%実験オペレーションの取得
prompt = {'Date:','Shot number:','start:','dt'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'230920','15','470','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
date = str2double(cell2mat(answer(1)));
IDXlist = str2num(cell2mat(answer(2)));
start = str2num(cell2mat(answer(3)));
dt = str2num(cell2mat(answer(4)));

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
n_data=numel(IDXlist);% 計測データ数
shotlist_a039 =T.a039(IDXlist);
shotlist_a040 = T.a040(IDXlist);
shotlist = [shotlist_a039, shotlist_a040];
tfshotlist_a039 =T.a039_TF(IDXlist);
tfshotlist_a040 =T.a040_TF(IDXlist);
tfshotlist = [tfshotlist_a039, tfshotlist_a040];
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);

trange=start:600;%【input】計算時間範囲
n=40; %【input】rz方向のメッシュ数

doCheck = false;
% doCheck = true;

figure('Position', [0 0 1500 1500],'visible','on');
for i=1:n_data
    dtacq_num=dtacqlist;
    shot=shotlist(i,:);
    tfshot=tfshotlist(i,:);
    if shot == tfshot
        tfshot = [0,0];
    end
    i_EF=EFlist(i);
    TF=TFlist(i);
    if doCheck
        check_signal(date, dtacq_num, shot, tfshot, pathname);
    else
        plot_psi200ch(date, shot, tfshot, pathname,n,i_EF,trange,IDXlist,dt);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function plot_psi200ch(date, shot, tfshot, pathname, n,i_EF,trange, IDXlist, dt)

filename = strcat(pathname.pre_processed_directory,'/a039_',num2str(shot(1)),'.mat');
if exist(filename,'file') == 0
    doCalculation = true;
else
    doCalculation = false;
end

if doCalculation
%較正係数のバージョンを日付で判別
sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));
C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
r_shift = 0.00;
ok = logical(C(:,14));
dtacq_num_list = C(:,1);
dtaq_ch = C(:,2);
polarity=C(:,13);
coeff=C(:,12);
zpos=C(:,9);
rpos=C(:,10)+r_shift;
ch=C(:,7);

if ismember(39,dtacq_num_list)
    filename1 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
    if exist(filename1,"file")==0
        disp('No rawdata file of a039 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
    end
    a039_raw = importdata(filename1);
end
if ismember(40,dtacq_num_list)
    filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
    if exist(filename2,"file")==0
        disp('No rawdata file of a040 -- Start generating!')
        rawdataPath = pathname.rawdata;
        save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
    end
    a040_raw = importdata(filename2);
end

raw = zeros(1000,length(dtaq_ch));
for i = 1:length(dtaq_ch)
    if dtacq_num_list(i) == 39
        raw(:,i) = a039_raw(:,dtaq_ch(i));
    elseif dtacq_num_list(i) == 40
        raw(:,i) = a040_raw(:,dtaq_ch(i));
    end
end

b=raw.*coeff';% 較正係数RC/NS
b=b.*polarity';% 極性揃え

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,100);
bt=bz;
ok_bz=false(100,1);
ok_bt=ok_bz;
zpos_bz=zeros(100,1);
rpos_bz=zpos_bz;
zpos_bt=zpos_bz;
rpos_bt=zpos_bz;

%digital filter
windowSize = 8;
bb = (1/windowSize)*ones(1,windowSize);
aa = 1;

for i=1:length(ch)
    b(:,i) = filter(bb,aa,b(:,i));
    b(:,i) = b(:,i) - mean(b(1:40,i));
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);
        zpos_bz(ceil(ch(i)/2))=zpos(i);
        rpos_bz(ceil(ch(i)/2))=rpos(i);
    elseif rem(ch(i),2)==0
        bt(:,ch(i)/2)=b(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);
        zpos_bt(ceil(ch(i)/2))=zpos(i);
        rpos_bt(ceil(ch(i)/2))=rpos(i);
    end
end

% zprobepcb    = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
zprobepcb    = [-0.2975,-0.255,-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17,0.255,0.2975];
rprobepcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift;
rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift;
[zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
% [zq,rq]      = meshgrid(zprobepcb,rprobepcb);
[zq_probepcb,rq_probepcb]=meshgrid(zprobepcb,rprobepcb);
ok_bt_matrix = false(length(rprobepcb),length(zprobepcb));
ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
for i = 1:length(ok_bt)
    if rpos_bt(i) > (r_shift)
        index_r = (abs(rpos_bt(i)-rprobepcb_t)<0.001);index_z = (zpos_bt(i)==zprobepcb);
        ok_bt_matrix = ok_bt_matrix + rot90(index_r,-1)*index_z*ok_bt(i);
    end
    index_r = (abs(rpos_bz(i)-rprobepcb)<0.001);index_z = (zpos_bz(i)==zprobepcb);
    ok_bz_matrix = ok_bz_matrix + rot90(index_r,-1)*index_z*ok_bz(i);
end

grid2D=struct(...
    'zq',zq,...
    'rq',rq,...
    'zprobepcb',zprobepcb,...
    'rprobepcb',rprobepcb,...
    'rprobepcb_t',rprobepcb_t,...1
    'ok_bz_matrix',ok_bz_matrix,...
    'ok_bt_matrix',ok_bt_matrix);
grid2D_probe = struct('zq',zq_probepcb,'rq',rq_probepcb,'rq_t',rprobepcb_t);

clear zq rq zprobepcb rprobepcb zq_probepcb rq_probepcb rprobepcb_t ok_bz_matrix ok_bt_matrix

% probecheck_script;

%data2Dcalc.m
r_EF   = 0.5 ;
n_EF   = 234. ;

if date<221119
    z1_EF   = 0.68;
    z2_EF   = -0.68;
else
    z1_EF   = 0.78;
    z2_EF   = -0.78;
end
[Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
clear EF r_EF n_EF i_EF z_EF

data2D=struct(...
    'psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Bt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jr',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Lambda',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'trange',trange);

% ******************* no angle correction ********************
for i=1:size(trange,2)
    t=trange(i);

    %Bzの二次元補間(線形fit)
    vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
    B_z = -Bz_EF+vq;
    B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);

    % PSI計算
    data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
    % data2D.psi(:,:,i) = flip(get_psi(flip(B_z,1),flip(grid2D.rq(:,1)),1),1);
    % このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bt(:,:,i)=B_t;
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
end

else
    load(filename,'data2D','grid2D');
end
% ***********************************************

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

% プロット部分
% figure('Position', [0 0 1500 1500],'visible','on');
start=0;
%  t_start=470+start;
 for m=1:16 %図示する時間
     i=start+m.*dt; %end
     t=trange(i);
     subplot(4,4,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none')
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none');% clim([0,max(data2D.Jt,[],'all')]);climError出た
    hold on;contour(grid2D.zq(1,:),grid2D.rq(:,1),-data2D.Jt(:,:,i),[2.3e5:1e5:7.5e5],'white','LineWidth',1);

%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')
    colormap("turbo")
    axis image
    axis tight manual
%     caxis([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
%     caxis([-0.01,0.01])%Bz
     % clim([-0.1,0.1])%Bt
    % clim([-5e-3,5e-3])%psi
%     caxis([-500,400])%Et
%     colorbar('Location','eastoutside')
    %カラーバーのラベル付け
%     c = colorbar;
%     c.Label.String = 'Jt [A/m^{2}]';
    hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.8e-3:40e-3],'black','LineWidth',0.5);xlim([-0.15,0.15]);
    % contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),'black','LineWidth',0.5);xlim([-0.15,0.15]);    
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
     % plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
    hold off
    title(string(t)+' us@shot' + string(IDXlist))
    drawnow
%     xlabel('z [m]')
%     ylabel('r [m]')
 end
 set(gcf,'Name','shot' + string(IDXlist),'NumberTitle','off');
 sgtitle('shot' + string(IDXlist))

if doCalculation
    clearvars -except data2D grid2D shot pathname;
    filename = strcat(pathname.pre_processed_directory,'/a039_',num2str(shot(1)),'.mat');
    save(filename);
end

end

function save_dtacq_data(dtacq_num,shot,tfshot,rawdataPath)

[rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
if tfshot==0
    [rawdata]=getMDSdata(dtacq_num,shot,0);
    save(strcat(rawdataPath,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata');
end

end

function check_signal(date, dtacq_num, shot, tfshot, pathname)
filename=strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
% filename=strcat(pathname.rawdata,'rawdata_noTF_dtacq',num2str(d_tacq),'.mat');
if exist(filename,"file")==0
    disp('No rawdata file -- Start generating!')
    rawdataPath = pathname.rawdata;
    save_dtacq_data(dtacq_num, shot, tfshot,rawdataPath)
    % return
end
load(filename,'rawdata');%1000×192

%正しくデータ取得できていない場合はreturn
if numel(rawdata)< 500
    return
end

%較正係数のバージョンを日付で判別
sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));

C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
ok = logical(C(:,14));
P=C(:,13);
coeff=C(:,12);
zpos=C(:,9);
rpos=C(:,10);
probe_num=C(:,5);
probe_ch=C(:,6);
ch=C(:,7);
d2p=C(:,15);
d2bz=C(:,16);
d2bt=C(:,17);

b=rawdata.*coeff';%較正係数RC/NS
b=b.*P';%極性揃え
b=smoothdata(b,1);

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,100);
bt=bz;
ok_bz=true(1,100);
ok_bt=ok_bz;

for i=1:192
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);
    elseif rem(ch(i),2)==0
        bt(:,ceil(ch(i)/2))=b(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);
    end
end
ok_bz_plot=ok_bz;
% ok_bt([4 5 6 7 8 9 10 15 21 27 30 42 43 49 53 69 84 87 92 94 95 96 97 98 99 100]) = false;

% 221219ver
% Pcheck=[1	-1	1	1	-1	-1	-1	1	1	1	1	1	1	-1	-1	0	-1	-1	1	-1	1	0	1	-1	-1	1	-1	-1	1	-1	-1	1	-1	1	-1	1	-1	-1	-1	1	-1	-1	1	1	-1	-1	-1	-1	-1	-1	1	-1	-1	-1	-1	-1	-1	-1	-1	-1	1	1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	1	-1	1	-1	-1	-1	-1	1	-1	-1	1	-1	-1	1	-1	-1	1	1	1	-1	1	-1	-1	1	1	1	-1];
% bz=bz.*Pcheck;
% bz(:,[37 47 57])=-bz(:,[37 47 57]);
% bz(:,70)=-bz(:,70);
% for i=0:9
%     bz(:,[7 8 9 10]+10.*i)=-bz(:,[7 8 9 10]+10.*i);
% end
% ok_bz([5 6 11 16 20 22 31 33 39 49 63 66 71 72 79 80 95 100])=false;
% ok_bz(49)=true;
% bz(:,49)=-bz(:,49);

%生信号描画用パラメータ
r = 5;%プローブ本数＝グラフ出力時の縦に並べる個数
col = 10;%グラフ出力時の横に並べる個数
y_upper_lim = 0.4;%3e-3;%0.1;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.4;%3e-3;%-0.1;%縦軸プロット領域（b_z下限）
t_start=1;%430;%455;%横軸プロット領域（開始時間）
t_end=1000;%550;%横軸プロット領域（終了時間）
% r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数

f1=figure;
f1.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz_plot(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe1-5')

f2=figure;
f2.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz_plot(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe6-10')

f3=figure;
f3.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)));
        xticks([t_start t_end]);
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bt signal probe1-5')

f4=figure;
f4.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)));
        xticks([t_start t_end]);
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bt signal probe6-10')

% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_02','.png'))
% close

%横軸z, 縦軸Bzのプロット
% f5=figure;
% f5.WindowState = 'maximized';
% t=465;
% for i=1:10
%     zline=(1:10:91)+(i-1);
%     bz_zline=bz(t,zline);
%     bz_zline(ok_bz(zline)==false)=NaN;
%     plot([-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17],bz_zline,'-*')
%     clear bz_zline
%     hold on
% end
% hold off
% xlabel('z [m]')
% ylabel('Bz')
% yline(0,'k--')
% title(strcat('t=',num2str(t),' us'))
% legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')

% bz1=[bz(t,1) bz(t,11) bz(t,21) bz(t,31) bz(t,41) bz(t,51) bz(t,61) bz(t,71) bz(t,81) bz(t,91)];
% bz2=[bz(t,2) bz(t,12) bz(t,22) bz(t,32) bz(t,42) bz(t,52) bz(t,62) bz(t,72) bz(t,82) bz(t,92)];
% bz3=[bz(t,3) bz(t,13) bz(t,23) bz(t,33) bz(t,43) bz(t,53) bz(t,63) bz(t,73) bz(t,83) bz(t,93)];
% bz4=[bz(t,4) bz(t,14) bz(t,24) bz(t,34) bz(t,44) bz(t,54) bz(t,64) bz(t,74) bz(t,84) bz(t,94)];
% bz5=[bz(t,5) bz(t,15) bz(t,25) bz(t,35) bz(t,45) bz(t,55) bz(t,65) bz(t,75) bz(t,85) bz(t,95)];
% bz6=[bz(t,6) bz(t,16) bz(t,26) bz(t,36) bz(t,46) bz(t,56) bz(t,66) bz(t,76) bz(t,86) bz(t,96)];
% bz7=[bz(t,7) bz(t,17) bz(t,27) bz(t,37) bz(t,47) bz(t,57) bz(t,67) bz(t,77) bz(t,87) bz(t,97)];
% bz8=[bz(t,8) bz(t,18) bz(t,28) bz(t,38) bz(t,48) bz(t,58) bz(t,68) bz(t,78) bz(t,88) bz(t,98)];
% bz9=[bz(t,9) bz(t,19) bz(t,29) bz(t,39) bz(t,49) bz(t,59) bz(t,69) bz(t,79) bz(t,89) bz(t,99)];
% bz10=[bz(t,10) bz(t,20) bz(t,30) bz(t,40) bz(t,50) bz(t,60) bz(t,70) bz(t,80) bz(t,90) bz(t,100)];

end

