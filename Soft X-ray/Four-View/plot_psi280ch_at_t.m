function plot_psi280ch_at_t(PCB,pathname)
% shot = PCB.shot;
trange = PCB.trange;
% start = PCB.start;
time = PCB.time;

[grid2D,data2D] = process_PCBdata_280ch(PCB,pathname);
% [grid2D,data2D] = process_PCBdata_200ch(PCB,pathname);

% ***********************************************

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

[magAxisList,xPointList] = get_axis_x_multi(grid2D,data2D); %時間ごとの磁気軸、X点を検索

% プロット部分
% figure('Position', [0 0 1500 1500],'visible','on');

figure;
i = find(trange==time);
% contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Br(:,:,i),40,'LineStyle','none')
contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none');
colorbar;
colormap(jet)
axis image
axis tight manual
hold on
contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),10,'black','LineWidth',3)
% contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.Br(:,:,i)),20,'black')
% plot(magAxisList.z(:,i),magAxisList.r(:,i),'ko');
% plot(xPointList.z(i),xPointList.r(i),'kx');
hold off
colorbar
xlabel('z [m]');
ylabel('r [m]');
title(string(time)+' us')
ax = gca;
ax.FontSize = 18;

% figure('Position', [0 0 1500 1500],'visible','on');
% start=30;
% dt = 4;
% %  t_start=470+start;
%  for m=1:16 %図示する時間
%      i=start+m.*dt; %end
%      t=trange(i);
%      subplot(4,4,m)
%     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none')
%     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),40,'LineStyle','none')
%     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
%     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')
% %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')
%     colormap(jet)
%     axis image
%     axis tight manual
% %     caxis([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
% %     caxis([-0.01,0.01])%Bz
%      % clim([-0.1,0.1])%Bt
%     % clim([-5e-3,5e-3])%psi
% %     caxis([-500,400])%Et
% %     colorbar('Location','eastoutside')
%     %カラーバーのラベル付け
% %     c = colorbar;
% %     c.Label.String = 'Jt [A/m^{2}]';
%     hold on
% %     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
% %     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%     % contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black','LineWidth',1)
% %     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
% %     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
%      % plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
% 
%     % timeIndex = find(trange==t);
%     % [magaxis,xpoint] = get_axis_x(grid2D,data2D,t);
%     % plot(magaxis.z,magaxis.r,'ro');
%     % plot(xpoint.z,xpoint.r,'rx');
%     plot(magAxisList.z(:,i),magAxisList.r(:,i),'ko');
%     plot(xPointList.z(i),xPointList.r(i),'kx');
% 
%     hold off
%     title(string(t)+' us')
% %     xlabel('z [m]')
% %     ylabel('r [m]')
%  end

 % sgtitle(strcat('shot',num2str(shot)));

 % drawnow

end