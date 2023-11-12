pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先

pathname.rawdata038=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所

pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所

%%%%実験オペレーションの取得
%直接入力の場合
dtacqlist=38;
shotlist=12344;%【input】dtacqの保存番号
tfshotlist=12335;
date = 230929;%【input】計測日
n_data=numel(shotlist);%計測データ数
i_EF = 210;%【input】EF電流
trange=450:500;%【input】計算時間範囲
n=50; %【input】rz方向のメッシュ数

for i=1:n_data
    dtacq_num=dtacqlist(i);
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    plot_psi125ch(date, dtacq_num, shot, tfshot, pathname,n,i_EF,trange); 
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function plot_psi125ch(date, dtacq_num, shot, tfshot, pathname, n,i_EF,trange)

% filename1=strcat(pathname.woTFdata,'rawdata_dtacq',num2str(shot(1)),'.mat');
filename1=strcat(pathname.rawdata,'rawdata_dtacq',num2str(dtacq_num(1)),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
load(filename1,'rawdata');%a038
rawdata1=rawdata;
clear rawdata

%正しくデータ取得できていない場合はreturn
if numel(rawdata1)< 500
    return
end

%較正係数のバージョンを日付で判別
sheets1 = sheetnames('coeff125ch.xlsx');
sheets1 = str2double(sheets1);
sheet_date1=max(sheets1(sheets1<=date));
C1 = readmatrix('coeff125ch.xlsx','Sheet',num2str(sheet_date1));

%a038
ok1 = logical(C1(:,14));
P1=C1(:,13);
coeff1=C1(:,12);
zpos1=C1(:,9);
rpos1=C1(:,10);
probe_num1=C1(:,5);
probe_ch1=C1(:,6);
ch1=C1(:,7);
p_ch= readmatrix('coeff125ch.xlsx','Sheet','p_ch');

b1=rawdata1.*coeff1';%較正係数RC/NS
b1=b1.*P1';%極性揃え
b1=double(b1);
% b1=smoothdata(b1,1);%b1=smoothdata(b1,1,'lowess',3);
%デジタイザchからプローブ通し番号順への変換
bz1=zeros(1000,126);
ok_bz1=true(100,1);
rpos_bz1=zeros(126,1);
zpos_bz1=rpos_bz1;

for i=1:128
    if ch1(i)>0
        bz1(:,ch1(i))=b1(:,i);
        ok_bz1(ch1(i))=ok1(i);
        rpos_bz1(ch1(i))=rpos1(i);
        zpos_bz1(ch1(i))=zpos1(i);
    end
end
bz1(:,63)=[];
ok_bz1(63)=[];
rpos_bz1(63)=[];
zpos_bz1(63)=[];

% ok_bz1([1 26 51 76 101])=false;

% for i=1:125
%     bz1(:,i)=lowpass(bz1(:,i),0.5e5,1e6);
% end

bz1=smoothdata(bz1,1);
% ok_bz1([5 54 31 8 61 37 62 13 38 16 40 20 45 47 23 48 25 28 29 30])=false;

[zq,rq]=meshgrid(linspace(min(zpos_bz1),max(zpos_bz1),n),linspace(min(rpos_bz1),max(rpos_bz1),n));
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

for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間(線形)
    vq = bz_rbfinterp(rpos_bz1, zpos_bz1, grid2D, bz1, ok_bz1, t);
    B_z = -Bz_EF+vq;
    %%PSI計算
    data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2.*pi.*grid2D.rq(:,1).*B_z,1);
    %このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
end
data2D.Et=diff(data2D.psi,1,3).*1e+6; 
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D.Et=-1.*data2D.Et./(2.*pi.*grid2D.rq);

ok_z = zpos_bz1(ok_bz1); %z方向の生きているチャンネル
ok_r = rpos_bz1(ok_bz1); %r方向の生きているチャンネル

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

figure('Position', [0 0 1500 1500],'visible','on');
start=15;
%  t_start=470+start;
 for m=1:10 %図示する時間
     i=start+m.*2; %end
     t=trange(i);
     subplot(2,5,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    caxis([-2*1e+6,1.0*1e+6]) %カラーバーの軸の範囲
%     caxis([-0.1,0.1])
    colorbar('Location','eastoutside')
    %カラーバーのラベル付け
%     c = colorbar;
%     c.Label.String = 'Jt [A/m^{2}]';
    hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
    % contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),50,'black')
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
    plot(ok_z,ok_r,"k.",'MarkerSize', 7)%測定位置
    hold off
    title(string(t)+' us')
    xlabel('z [m]')
    ylabel('r [m]')
 end

%  saveas(gcf,strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230119\jt_a038\a038_',num2str(shot),'_jt_t466us.png'))
%  close

end
