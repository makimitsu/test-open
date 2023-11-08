%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%プロット時間を指定して磁気面を1枚プロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_psi200ch_at_t(time,trange,grid2D,data2D,ok_z,ok_r)

figure('Position', [0 0 1500 1500],'visible','on');
i=(time-trange(1)+1);
t=trange(i);
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),-20e-3:0.2e-3:40e-3,'black','LineWidth',1)
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
hold on
% plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
hold on
title(string(t)+' us')
xlabel('z [m]')
ylabel('r [m]')
ax = gca;
ax.FontSize = 14;

% figure('Position', [0 0 1500 1500],'visible','on');
% for m=1:16 %図示する時間
%     i=(start-trange(1)+1)+(m-1)*dt; %end
%     t=trange(i);
%     subplot(4,4,m)
%     plot(grid2D.zq(1,:),data2D.Bt(25,:,i))
%     title(string(t)+' us')
%     ylim([-0.03 0.03])
% end
% saveas(gcf,strcat('/Users/yunhancai/Downloads/files/a039_',num2str(shot)))
% saveas(gcf,strcat('/Users/yunhancai/Downloads/a039_',num2str(shot),'.png'))
% close
end
