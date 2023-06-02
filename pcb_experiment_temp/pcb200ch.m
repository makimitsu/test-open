%%%%%%%%%%%%%%%%%%%%%%%%
%200ch用新規pcbプローブのみでの磁気面（Bz）
%dtacqのshot番号を直接指定する場合
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



% %直接入力の場合
 dtacqlist=39;
 shotlist=645;%【input】dtacqの保存番号
 tfshotlist=584;%【input】dtacqのTFのみ番号;
 date_list = 230127;%【input】計測日
 n_data=numel(shotlist);%計測データ数
 EFlist = 160;%【input】EF電流

 trange=400:600;%【input】計算時間範囲
 n=40; %【input】rz方向のメッシュ数

%{
% %自動入力の場合
log_cell = readcell('/Users/yunhancai/Google Drive/Data/log2022-2023.xlsx','Sheet','TOTAL');
dtacqlist=39;
shotlist=896:1398;%【input】dtacqの保存番号
shotlist_index = cell2mat(log_cell(2:end,4));
[sharedvals,idx] = intersect(shotlist_index,shotlist,'stable');
tfshotlist=cell2mat(log_cell(idx+1,6));% TF offset
date_list = cell2mat(log_cell(idx+1,2));%【input】計測日
n_data=numel(sharedvals);%計測データ数
EFlist = cell2mat(log_cell(idx+1,19));%【input】EF電流
%}



for i=1:n_data
    dtacq_num=dtacqlist;
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    i_EF=EFlist(i);
    date = date_list(i);
    plot_psi200ch(date, dtacq_num, shot, tfshot, pathname,n,i_EF,trange);
    disp(['pcb:',num2str(i),'/',num2str(n_data)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function plot_psi200ch(date, dtacq_num, shot, tfshot, pathname, n,i_EF,trange)
filename=strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
if exist(filename,"file")==0
    disp('File does not exit');
    return
end
load(filename,'rawdata');%1000×192

%正しくデータ取得できていない場合はreturn
if numel(rawdata)< 500
    disp('Unable to extract data from file');
    return
end

%較正係数のバージョンを日付で判別
sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));

C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
r_shift = 0.00;
ok = logical(C(:,14));
P=C(:,13);
coeff=C(:,12);
zpos=C(:,9);
rpos=C(:,10)+r_shift;
probe_num=C(:,5);
probe_ch=C(:,6);
ch=C(:,7);
d2p=C(:,15);
d2bz=C(:,16);
d2bt=C(:,17);

b=rawdata.*coeff';%較正係数RC/NS
b=b.*P';%極性揃え
b=smoothdata(b,1);%Smoothing

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
windowSize = 3;
bb = (1/windowSize)*ones(1,windowSize);
aa = 1;

for i=1:192
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

ok_bt([4 5 6 7 8 9 10 94 95 96 97 98 99 100]) = false;
zprobepcb    = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
rprobepcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift;
rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift;
[zq,rq]      = meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
% [zq,rq]      = meshgrid(zprobepcb,rprobepcb);
[zq_probepcb,rq_probepcb]=meshgrid(zprobepcb,rprobepcb);
ok_bt_matrix = false(length(zprobepcb),length(rprobepcb));
ok_bz_matrix = false(length(zprobepcb),length(rprobepcb));
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
    'rprobepcb_t',rprobepcb_t,...
    'ok_bz_matrix',ok_bz_matrix,...
    'ok_bt_matrix',ok_bt_matrix);
grid2D_probe = struct('zq',zq_probepcb,'rq',rq_probepcb,'rq_t',rprobepcb_t);

clear zq rq zprobepcb rprobepcb zq_probepcb rq_probepcb rprobepcb_t ok_bz_matrix ok_bt_matrix

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
    'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Lambda',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'trange',trange);

% **************** angle correction **************** 
sheets_angle = sheetnames('angle.xlsx');sheets_angle = str2double(sheets_angle);
sheet_angle_date=max(sheets_angle(sheets_angle<=date));
angle_file = readmatrix('angle.xlsx','Sheet',num2str(sheet_angle_date));angle = angle_file(2,2:end);
sin_matrix = repmat(sin(angle),10,1);%sin_matrix(1:10,:) = zeros(10,10);
cos_matrix = repmat(cos(angle),10,1);%cos_matrix(1:10,:) = ones(10,10);
B_z_calibrated = zeros(length(grid2D.rprobepcb),length(grid2D.zprobepcb),length(trange));
B_t_calibrated = zeros(length(grid2D.rprobepcb_t),length(grid2D.zprobepcb),length(trange));
for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間と角度調整(線形fit)
    B_z_noncalib = bz_rbfinterp(rpos_bz, zpos_bz, grid2D_probe, bz, ok_bz, t);
    B_t_noncalib = bz_rbfinterp(rpos_bt, zpos_bt, grid2D_probe, bt, ok_bt, t);
    B_t_calibrated(:,:,i) = B_t_noncalib.*cos_matrix+B_z_noncalib.*sin_matrix;
    B_z_calibrated(:,:,i) = -B_t_noncalib.*sin_matrix+B_z_noncalib.*cos_matrix;
end
B_z_calibrated_restored = bz;B_t_calibrated_restored = bt;
for i=1:192
    if rem(ch(i),2)==1
        z_index = zpos(i)==grid2D.zprobepcb;
        r_index = rpos(i)==grid2D.rprobepcb;
        B_z_calibrated_restored(trange,ceil(ch(i)/2)) = reshape(B_z_calibrated(r_index,z_index,:),[],1);
    elseif rem(ch(i),2)==0
        z_index = zpos(i)==grid2D.zprobepcb;
        r_index = (0.001>abs(rpos(i)-(grid2D.rprobepcb_t)));
        B_t_calibrated_restored(trange,ch(i)/2) = reshape(B_t_calibrated(r_index,z_index,:),[],1);
    end
end
for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間(線形fit)
    vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, B_z_calibrated_restored, ok_bz, t);
    B_z = -Bz_EF+vq;
    B_t = bz_rbfinterp(rpos_bt-0.01, zpos_bt, grid2D, B_t_calibrated_restored, ok_bt, t);

    %%PSI計算
%     data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
    data2D.psi(:,:,i) = flip(get_psi(flip(B_z,1),flip(grid2D.rq(:,1)),1),1);
    
    %このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bt(:,:,i)=B_t;
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
    data2D.Lambda(:,:,i) = (2*pi*grid2D.rq.*data2D.Bt(:,:,i))./(data2D.psi(:,:,i));
end
% ************************************************


% ******************* old way ********************
% for i=1:size(trange,2)
%     t=trange(i);
% 
%     %Bzの二次元補間(線形fit)
%     vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
%     B_z = -Bz_EF+vq;
%     B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);
% 
%     %PSI計算
%     data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
%     このままだと1/2πrが計算されてないので
%     [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
%     data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
%     data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
%     data2D.Bt(:,:,i)=B_t;
%     data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
%     data2D.Lambda(:,:,i) = (2*pi*grid2D.rq.*data2D.Bt(:,:,i))./(data2D.psi(:,:,i));
% end
% ***********************************************

%時間方向(次元3)に差分をとる。dt = 1e-6で割る。
data2D.Et=diff(data2D.psi,1,3).*1e+6;
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D.Et=-1.*data2D.Et./(2.*pi.*grid2D.rq);


if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

clearvars -except data2D grid2D shot data2D.Et;
filename = strcat('C:\Users\uswk0\OneDrive - g.ecc.u-tokyo.ac.jp\data\before_picture\a039_',num2str(shot),'.mat');
save(filename)
end