%%%2022/11/16~

pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fouier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path'); %保存先

%【input】date, i_EF(EF電流), d_tacq
date=221119;%T.date(IDX);
i_EF=150;%150;
d_tacq=10194;%T.d_tacq(IDX);

%【input】とりあえず変更しなくて良い
shot='';%T.shot(IDX);
TF_shot='';%T.TFoffset(IDX);
offset_TF='';%isfinite(TF_shot);

d_tacqTF='';%T.TFdtacq(IDX);
trange=460:520;
n=50; %rz方向のメッシュ数

%pcbdata.m
% pathname.rawdata='C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\rawdata_a038_noTF\'; %rawdataの保管場所
% filename=strcat(pathname.rawdata,'rawdata_noTF_dtacq',num2str(d_tacq),'.mat');
pathname.rawdata='C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\rawdata_a038\'; %rawdataの保管場所
filename=strcat(pathname.rawdata,'rawdata_dtacq',num2str(d_tacq),'.mat');
load(filename,'rawdata');
if numel(rawdata)< 500
    grid2D=NaN;
    data2D=NaN;
    return
end
load('rc_coeff2020.mat','coeff')
% if exist(coeff)==0
%     load('rc_coeff2020.mat')
% end

%getpcbbz.m
% データの並べ替え（積分器のチャンネルとデジタイザの対応）
d=reshape([1:32],2,16);%デジタイザの順番
c=[d;d+32];%積分器の順番
c=c(:);
d_ch=[c; c+64]';
[ch, c2d, d2c] = unique(d_ch);%d2cはデジタイザ順の配列を積分器順に直し,c2dはその逆
clear c d
% 無視するチャンネル(積分器＆デジタイザ)
ng1=var(rawdata(1:50,:),1);%0~50usの中で分散がゼロはショートしている？
%400~600usの中で明らかに巨大な信号はノイズ？
ng2 =sum(abs(rawdata(420:490,:)),1)>=10;
ng4 = sum(abs(rawdata(420:490,:)),1) <=0.3;
%構成係数の外れ値を検出して積分器のおかしいチャンネルを抜く
ng3= isoutlier(coeff,'median');
ok= ng1~=0 & ng2==0 & ng3 == 0 & ng4==0;
clear ng1 ng2 ng3 ng4
% 校正係数をかけてbzに直す。今回の構成係数はデジタイザごとになっているので変換の必要はない
bz=rawdata.*coeff;
clear rawdata coeff

% ch位置の読み込みとチャンネルと位置の対応
r_ch=25;
rind=0.055+0.005*[1,5,9,13,17,[19:41],43,47];
%load('r_pcb2019.mat')
[rpos,zpos] = meshgrid(rind(1:r_ch),[-0.0525 -0.021 0 0.021 0.0525]);
%p_ch(1,:)=[1:125];
for i=1:5
    p_ch(i,:)=[1:r_ch]+r_ch*(i-1);
end
p_ch(p_ch>=63)=p_ch(p_ch>=63)+1;
[ch_p, p2c, c2p] = unique(p_ch(:));
rpos=rpos(p2c);
zpos=zpos(p2c);
clear rind

% 無視するチャンネル(プローブ)
% ok=ok(d2c);
% ok=ok(ch(ch_p));
Probe=readmatrix('probe_version.xlsx','Sheet','22111802');
Probe(63,:)=[];
ok=Probe(1:125,4);
ok=logical(ok);

bz=bz(:,d2c);
bz=bz(:,ch(ch_p));
clear d2c ch_p okp c2p d2c d_ch p2c

% Bzの極性をそろえる（TF only）:Bz(time1000×ch125)
P=Probe(1:125,5);%125*1の±1の極性
for i=1:125
    bz(:,i)=P(i).*bz(:,i);
end

% チャンネルごとの生信号のプロット
bz=smoothdata(bz,1);%各ch内のデータを時間に関してsmoothing

ok([9 10 51 84 87 93 109 110 111])=false;%積分器故障？
% bz(:,9)=(2.*bz(:,8)+bz(:,11))./3;
% bz(:,10)=(bz(:,8)+2.*bz(:,11))./3;
% bz(:,84)=(bz(:,83)+bz(:,85))./2;
% bz(:,87)=(bz(:,86)+bz(:,88))./2;
% bz(:,93)=(bz(:,92)+bz(:,94))./2;
% bz(:,109)=(bz(:,108)+3.*bz(:,112))./4;
% bz(:,110)=2.*(bz(:,108)+2.*bz(:,112))./4;
% bz(:,111)=(bz(:,108)+3.*bz(:,112))./4;
% 
% bz(:,4)=(bz(:,3)+bz(:,5))./2;
% bz(:,39)=(bz(:,38)+bz(:,40))./2;
% bz(:,53)=(bz(:,52)+bz(:,54))./2;
% bz(:,64)=(bz(:,63)+bz(:,65))./2;
% bz(:,71)=(bz(:,70)+bz(:,72))./2;
% bz(:,104)=(bz(:,103)+bz(:,105))./2;
% bz(:,107)=(bz(:,106)+bz(:,108))./2;
% 
% bz(:,35)=(2.*bz(:,34)+bz(:,37))./3;
% bz(:,36)=(bz(:,34)+2.*bz(:,37))./3;
% bz(:,58)=(2.*bz(:,57)+bz(:,60))./3;
% bz(:,59)=(bz(:,57)+2.*bz(:,60))./3;
% 
% bz(:,78)=(bz(:,76)+bz(:,80))./2;
% bz(:,79)=(bz(:,76)+3.*bz(:,80))./4;
% 
% bz(:,50)=(bz(:,25)+bz(:,75))./2;
% bz(:,49)=(bz(:,24)+bz(:,74))./2;

%移動平均
% bz=movmean(bz,5,1);

% r = 5;%プローブ本数＝グラフ出力時の縦に並べる個数
% col1 = 12;%1枚目のグラフ出力時の横に並べる個数
% col2 = 13;%2枚目のグラフ出力時の横に並べる個数
% y_upper_lim = 0.03;%0.1;%縦軸プロット領域（b_z上限）
% y_lower_lim = -0.03;%-0.1;%縦軸プロット領域（b_z下限）
% t_start=330;%455;%横軸プロット領域（開始時間）
% t_end=600;%横軸プロット領域（終了時間）
% r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数
% plotbzsignal(y_upper_lim, col2, col1, t_end, p_ch, y_lower_lim, t_start, bz, ok, r_ch, r);
% clear r col1 col2 y_upper_lim y_lower_lim t t_start t_end

[zq,rq]=meshgrid(linspace(min(zpos),max(zpos),n),linspace(0.1,max(rpos),n));%計測からr=0.06と0.08mを除く
% [zq,rq]=meshgrid(linspace(min(zpos),max(zpos),n),linspace(min(rpos),max(rpos),n));
grid2D=struct('zq',zq,'rq',rq);
clear zq rq

%data2Dcalc.m
r_EF   = 0.5 ;
n_EF   = 234. ;
%i_EF    =i_EF;

if date<221119
    z1_EF   = 0.875;%0.68;
    z2_EF   = -0.830;%-0.68;
else
    z1_EF   = 0.78;
    z2_EF   = -0.78;
end

[Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
clear EF r_EF n_EF i_EF z_EF

data2D=struct('psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'trange',trange);
minpsi=zeros(size(grid2D.rq,1),size(trange,2));
minind=zeros(size(grid2D.rq,1),size(trange,2));
for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間（rbfinterp）
    %計測からr=0.06と0.08mを除く
    ok([1 2 26 27 51 52 76 77 101 102])=false;
    %ok([1 2 3 4 26 27 28 29 51 52 53 54 76 77 78 79 101 102 103 104])=false;
    vq = bz_rbfinterp(rpos, zpos, grid2D, bz, ok, t);

    B_z = -Bz_EF+vq;
    %%PSI計算
    data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),B_z,1);
    %このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
end
data2D.Et=diff(data2D.psi,1,3).*1e+6; 
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D.Et=data2D.Et./(2.*pi.*grid2D.rq);
%data2D=struct('psi',psi,'Bz',Bz,'Br',Br,'Jt',Jt,'trange',trange);

ok_z = zpos(ok); %z方向の生きているチャンネル
ok_r = rpos(ok); %r方向の生きているチャンネル

clear vq bz Bz_EF B_z ok ng rpos zpos i t

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

figure('Position', [0 0 1500 1500],'visible','on');
start=10; %460+?
%  t_start=470+start;
 for m=1:10 %図示する時間
     i=start+m.*2; %end
     t=trange(i);
     subplot(2,5,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),10,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    caxis([-2.5*1e+6,2.5*1e+6]) %カラーバーの軸の範囲
%     caxis([-1e-1,1e-1])
    %caxis([-maxrange,maxrange])
    colorbar('Location','eastoutside')
    %カラーバーのラベル付け
%     c = colorbar;
%     c.Label.String = 'Jt [A/m^{2}]';
    hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),30,'black')
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
%     plot(ok_z,ok_r,"k.",'MarkerSize', 7)%測定位置
    hold off
    title(string(t)+' us')
    xlabel('z [m]')
    %ylabel('r [m]')
%     ylim([0.1 grid2D.rq(end,1)])
%     xlim([-0.04 0.04])
%     ax = gca; %y軸を消す
%     ax.YTickLabel = cell(size(ax.YTickLabel)); 
 end

