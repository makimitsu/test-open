close all

%各PCのパスを定義
run define_path.m

date = 230524;%【input】実験日
dtacq.num = 39;
dtacq.shot = 1462;
dtacq.tfshot = 1444;
time = 475;
trange = 430:590;%【input】磁気プローブ計算時間範囲(430:590)
n_CH = 28;
n_z = 1;
factor = 0.15;
cmap = 'Bt';
cbar = true;
plot_list_1 = [2 5 6 8];
plot_list_2 = [2 4 6 9];

unit = 1e2;

[grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
[~,data2D_tfoffset,~,~] = load_pcb200ch_with_tfoffset(date,dtacq,pathname);

% figure('Position', [0 0 1500 600],'visible','on');
figure('Position', [0 0 600 600],'visible','on');
i=(time-trange(1)+1);
t=trange(i);
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
contour(grid2D.zq(1,:)*unit,grid2D.rq(:,1)*unit,squeeze(data2D.psi(:,:,i)),-20e-3:0.2e-3:40e-3,'black','LineWidth',1)
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
hold on
% plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
hold on

switch cmap
    case 'Bz'
        contourf(grid2D.zq(1,:)*unit,grid2D.rq(:,1)*unit,data2D.Bz(:,:,i)*1e3,100,'LineStyle','none','FaceAlpha',0.5)%Bz
        clim([-30,30])%Bz
        if cbar
            colorbar('Location','eastoutside')
            c = colorbar;
            c.Label.String = 'Bz [mT]';
        end
        % colormap(jet)
    case 'Br'
        contourf(grid2D.zq(1,:)*unit,grid2D.rq(:,1)*unit,data2D.Br(:,:,i)*1e3,100,'LineStyle','none','FaceAlpha',0.5)%Br
        clim([-50,50])%Br
        if cbar
            colorbar('Location','eastoutside')
            c = colorbar;
            c.Label.String = 'Br [mT]';
        end
        % colormap(jet)
    case 'psi'
        contourf(grid2D.zq(1,:)*unit,grid2D.rq(:,1)*unit,data2D.psi(:,:,i)*1e3,100,'LineStyle','none','FaceAlpha',0.5)%psi
        clim([0,4])%psi
        if cbar
            colorbar('Location','eastoutside')
            c = colorbar;
            c.Label.String = 'Psi [Wb]';
        end
        colormap(jet)
    case 'Bt'
        contourf(grid2D.zq(1,:)*unit,grid2D.rq(:,1)*unit,data2D_tfoffset.Bt(:,:,i)*1e3,100,'LineStyle','none','FaceAlpha',0.5)%Bt
        clim([100,400])%Bt
        if cbar
            colorbar('Location','eastoutside')
            c = colorbar;
            c.Label.String = 'Bt [mT]';
        end
        % colormap(jet)
    case 'Jt'
        contourf(grid2D.zq(1,:)*unit,grid2D.rq(:,1)*unit,-1.*data2D.Jt(:,:,i),100,'LineStyle','none','FaceAlpha',0.5)%jt
        clim([-0.8*1e+6,0.8*1e+6]) %jt
        if cbar
            colorbar('Location','eastoutside')
            c = colorbar;
            c.Label.String = 'Jt [A/m^{2}]';
        end
        colormap(jet)
    case 'Et'
        contourf(grid2D.zq(1,:)*unit,grid2D.rq(:,1)*unit,-1.*data2D.Et(:,:,i),100,'LineStyle','none','FaceAlpha',0.5)%Et
        clim([-300,0])%Et
        if cbar
            colorbar('Location','eastoutside')
            c = colorbar;
            c.Label.String = 'Et [V/m]';
        end
        colormap(jet)
    otherwise
        close(fig)
        warning('Input error in cmap.')%cmapの入力エラー
        return;
end
axis image
axis tight manual
hold on

%ドップラープローブ計測点配列を生成

load([pathname.mat,'/ionflow_magpres/230524/shot15-24.mat'],'mpoints','V_i_all','T_i_all');
mpoints_1 = mpoints;
V_i_all_1 = V_i_all(:,:,plot_list_1);
T_i_all_1 = T_i_all(:,:,plot_list_1);

load([pathname.mat,'/ionflow_magpres/230526/shot2-12.mat'],'mpoints','V_i_all','T_i_all');
mpoints_2 = mpoints;
V_i_all_2 = V_i_all(:,:,plot_list_2);
T_i_all_2 = T_i_all(:,:,plot_list_2);

mpoints_1.r = [mpoints_1.r(1,1); mpoints_2.r];
mpoints_1.z= [mpoints_1.z(1,1); mpoints_2.z];
V_i_all_2(2,:,:) = V_i_all_1(3,:,:);
V_i_all_2 = [V_i_all_1(1,:,:);V_i_all_2];
V_i_all_1 = V_i_all_2;
T_i_all_2(2,:,:) = T_i_all_1(3,:,:);
T_i_all_2 = [T_i_all_1(1,:,:);T_i_all_2];
T_i_all_1 = T_i_all_2;

V_i_mean_1 = round(mean(V_i_all_1(:,:,:),3),1);
V_i_sigma_1 = std(V_i_all_1(:,:,:),0,3);
T_i_mean_1 = mean(T_i_all_1(:,1,:),3);
T_i_sigma_1 = std(T_i_all_1(:,1,:),0,3);

mpoints_1.z = mpoints_1.z;
mpoints_1.r = mpoints_1.r;
factor = factor;

% %温度カラーマップ
% T_icon = repmat(T_i(:,1),1,5);%等高線図用平均イオン温度
% zcon = mpoints.z(1,1);
% s = pcolor([zcon-2*unit zcon-1*unit zcon zcon+1*unit zcon+2*unit],mpoints.r,T_icon);
% s.FaceColor = 'interp';
% s.EdgeAlpha = 0;
% colormap('jet')
% colorbar
% clim([0 150])
% hold on

r_start = 1;
r_end = 6;
% plot(mpoints.z,mpoints.r,'xr');%ドップラープローブ計測点を表示
plot(mpoints_1.z(r_start:r_end,1),mpoints_1.r(r_start:r_end,1),'xk');%ドップラープローブ計測点を表示
hold on
if mpoints_1.n_z == 1
    q = quiver(mpoints_1.z(r_start:r_end,1),mpoints_1.r(r_start:r_end,1),V_i_mean_1(r_start:r_end,1)*factor,V_i_mean_1(r_start:r_end,2)*factor);
end
q.LineWidth = 2;
q.MaxHeadSize = 10;
q.AutoScale = 'off';
q.Color = 'r';
% for i = r_start:r_end
%     txt = text(mpoints_1.z(i,1)+0.5,mpoints_1.r(i,1)-0.2,['(',num2str(V_i_mean_1(i,1)),', ',num2str(V_i_mean_1(i,2)),')']);
%     txt.FontSize = 18;
%     txt.Color = 'r';
%     txt.FontWeight = 'bold';
% end
% xlim([-17 17])
% ylim([0.06 0.33])
xlim([-3.0 7.0])
ylim([6 25])
xticks(-3.0:2.0:7.0)
yticks(7.5:2.5:25)
daspect([1 1 1])
ax = gca;
ax.FontSize = 18;
xlabel('Z [cm]')
ylabel('R [cm]')