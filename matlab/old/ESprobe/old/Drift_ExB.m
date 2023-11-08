
%%% cal & plot ExB drift velocity %%%

clear all

pathname.rawdata='C:\Users\Desktop\mat\pcb_raw';

% %直接入力の場合
dtacqlist=39;
shotlist=1918;              % dtacqの保存番号
tfshotlist=0;               
tfshotlist2=1903;           
date = 230712;              % 計測日
n_data=numel(shotlist);     % 計測データ数
EFlist = 150;               % EF電流
TFlist = 0;
trange=430:590;%【input】計算時間範囲(ほぼ固定)
n=50; %【input】rz方向のメッシュ数
cmap = false;%【input】カラーマップの有無

for i=1:n_data
    dtacq_num=dtacqlist;
    shot=shotlist(i);
    tfshot=tfshotlist;
    tfshot2=tfshotlist2;
    i_EF=EFlist;
    TF=TFlist;
    plot_psi200ch(date, dtacq_num, shot, tfshot,tfshot2,pathname,n,i_EF,trange,TF,cmap);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%


function plot_psi200ch(date, dtacq_num, shot, tfshot,tfshot2, pathname, n,i_EF,trange,TF,cmap)
Nplot =1;%【input】プロット枚数(16以下の4の倍数or3以下)
start = 474;%【input】プロット開始時間%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1;%【input】プロット時間間隔

filename=strcat(pathname.rawdata,'/',num2str(date),'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
if exist(filename,"file")==0
    disp(strcat(filename,' does not exist.'));
    return
end
load(filename,'rawdata0');%真空磁場含む

filename2=strcat(pathname.rawdata,'/',num2str(date),'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot2),'.mat');
if exist(filename2,"file")==0
    disp(strcat(filename2,' does not exist.'));
    return
end
load(filename2,'rawdata');%真空磁場差し引き


%===============================================================================
% phi, Ez, Er
%===============================================================================
z = linspace(-0.15,0.15,50);
r = linspace(0.11,0.285,50);
z_probe=linspace(-0.15,0.15,21);
r_probe=linspace(0.11,0.285,8);
filepath = strcat...
    ('C:\Users\msi\Desktop\electroprobe\230712\shot\toroidal_mode20230712_No_pb','.xlsx');

phi_new = zeros(8,20); %% （スキャン回数、生きたチャンネル数）
for i = 1:8
    phi_new(i,:)= readmatrix(filepath,'Sheet',num2str(i),...
        'Range','B4740:U4740');
end
phi = fliplr(phi_new).*50;
z_probe = z_probe([1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]);%%% 生きたチャンネル
[phi_mesh_z,phi_mesh_r] = meshgrid(z,r);
phi_grid = griddata(z_probe,r_probe,phi,phi_mesh_z,phi_mesh_r);

Ez = zeros(8,20); 
Er = zeros(8,20); 

for i = 1:8
    Ez(i,:) = -gradient(phi(i,:),z_probe);
end
Ez_grid = griddata(z_probe,r_probe,Ez,phi_mesh_z,phi_mesh_r);

for j = 1:20 
    Er(:,j) = -gradient(phi(:,j),r_probe);
end
Er_grid = griddata(z_probe,r_probe,Er,phi_mesh_z,phi_mesh_r);
%===============================================================================


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

b=rawdata0.*coeff';%較正係数RC/NS TF
b=b.*P';%極性揃え
b=smoothdata(b,1);

b2=rawdata.*coeff';%較正係数RC/NS PSI
b2=b2.*P';%極性揃え
b2=smoothdata(b2,1);

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,100);
bt=bz;
bz2=zeros(1000,100);
bt2=bz2;

ok_bz=false(100,1);
ok_bt=ok_bz;

zpos_bz=zeros(100,1);
rpos_bz=zpos_bz;
zpos_bt=zpos_bz;
rpos_bt=zpos_bz;

for i=1:192
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        bz2(:,ceil(ch(i)/2))=b2(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);

        zpos_bz(ceil(ch(i)/2))=zpos(i);
        rpos_bz(ceil(ch(i)/2))=rpos(i);
    elseif rem(ch(i),2)==0
        bt(:,ch(i)/2)=b(:,i);
        bt2(:,ch(i)/2)=b2(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);

        zpos_bt(ceil(ch(i)/2))=zpos(i);
        rpos_bt(ceil(ch(i)/2))=rpos(i);
    end
end

[bz, ok_bz, ok_bz_plot] = ng_replace(bz, ok_bz, sheet_date);
[bz2, ok_bz, ok_bz_plot] = ng_replace(bz2, ok_bz, sheet_date);

ok_bt([3 4 5 6 7 8 9 10 15 21 27 30 42 43 49 53 69 84 87 92 94 95 96 97 98 99 100]) = false;

[zq,rq]=meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));

grid2D=struct('zq',zq,'rq',rq);
clear zq rq

%data2Dcalc.m
r_EF   = 0.5 ;
n_EF   = 234. ;

if date<221119
    z1_EF   = 0.875;%0.68;
    z2_EF   = -0.830;%-0.68;
else
    z1_EF   = 0.78;
    z2_EF   = -0.78;
end
[Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
clear EF r_EF n_EF i_EF z_EF

data2D=struct('psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'psi2',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Bt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'trange',trange);

for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間(線形fit)
    vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
    B_z = -Bz_EF+vq;
    B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);

    vq2 = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz2, ok_bz, t);
    B_z2 = -Bz_EF+vq2;
    B_t2 = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt2, ok_bt, t);
    
    %%PSI計算
    data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
    data2D.psi2(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z2.*grid2D.rq(:,1),1);

    %このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bt(:,:,i)=B_t;
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);

    % 磁場データと静電データ間の補間
    new_grid_Br = griddata(grid2D.zq(1,:), grid2D.rq(:,1), data2D.Br(:,:,i), phi_mesh_z, phi_mesh_r);
    new_grid_Bz = griddata(grid2D.zq(1,:), grid2D.rq(:,1), data2D.Bz(:,:,i), phi_mesh_z, phi_mesh_r);
    new_grid_Bt = griddata(grid2D.zq(1,:), grid2D.rq(:,1), data2D.Bt(:,:,i), phi_mesh_z, phi_mesh_r);
end
data2D.Et=diff(data2D.psi,1,3).*1e+6;
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D.Et=data2D.Et./(2.*pi.*grid2D.rq);
new_grid_Et = griddata(grid2D.zq(1,:), grid2D.rq(:,1), data2D.Et(:,:,start-430), phi_mesh_z, phi_mesh_r);

absB_grid = new_grid_Bt.^2 + new_grid_Br.^2 + new_grid_Bz.^2;

ok_z = zpos_bz(ok_bz_plot); %z方向の生きているチャンネル
ok_r = rpos_bz(ok_bz_plot); %r方向の生きているチャンネル

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end


figure('Position', [0 0 1500 1500],'visible','on');
%  t_start=470+start;
for m=1:Nplot %図示する時間
    i=(start-trange(1)+1)+(m-1)*dt; %end
    t=trange(i);
    if Nplot > 1
        if Nplot == 4
            subplot(2,2,m)
        elseif mod(Nplot,4) == 0
            subplot(Nplot/4,4,m)
        elseif Nplot < 4
            subplot(1,Nplot,m)
        end
    end

    vz = (Er_grid.*new_grid_Bt - new_grid_Br.*new_grid_Et)./absB_grid;  % Vz
    vr = -(Ez_grid.*new_grid_Bt - new_grid_Bz.*new_grid_Et)./absB_grid; % Vr 
    vd = sqrt(vz.^2 + vr.^2)*1e-3;                                      %% |V|

    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi2(:,:,i)),[-20e-3:0.2e-3:40e-3],'black','LineWidth',1)

hold on
    
    title(string(t)+' us','FontSize',16)
    xlabel('z [m]','FontSize',14)
    ylabel('r [m]','FontSize',14)

    set(gca,'FontSize',14)
    xlim([-0.15 0.15])
    ylim([0.11 0.29])
    axis square

%============== fig drift vector==============================
    q = quiver(phi_mesh_z,phi_mesh_r,vz,vr,2);
    q.LineWidth = 2;
    q.MaxHeadSize = 3;
    q.Color = 'r';

%%% fig potential & overlay with psi  %%%
    contour_phi = linspace(min(phi_grid,[],'all'),max(phi_grid,[],'all'),100);
    [c,h] = contourf(phi_mesh_z,phi_mesh_r,phi_grid,contour_phi,'Fill','on','edgecolor','none');    
    colormap(redblue(100));
    clb = colorbar;
    clim([-150 150])
    clb.Label.String = '\phi_{f} [V]';
    clb.Label.FontSize = 16;
    h.FaceAlpha = 0.6;


hold off

end
end
