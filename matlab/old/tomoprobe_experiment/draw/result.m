% save('shot15.mat','Vx','Vy','drReF')

%再構成速度分布関数を描画
% fprintf('Summit is (Vz = %.1f, Vr = %.1f)[km/s].\n',Vx(Summitx),Vy(Summity));
levels = [drReF(Summity,Summitx)*0.1, drReF(Summity,Summitx)*0.2, ...
    drReF(Summity,Summitx)*0.3, drReF(Summity,Summitx)*0.4, ...
    drReF(Summity,Summitx)*0.5, drReF(Summity,Summitx)*0.6, ...
    drReF(Summity,Summitx)*0.7, drReF(Summity,Summitx)*0.8, ...
    drReF(Summity,Summitx)*0.9];%等高線本数
figure('Position',[700 150 550 500])
% figure('Position',[700 150 550 500],'visible','off')
if draw_type == 1
    contourf(Vy, Vx, drReF, levels)%2次元等高線図
elseif draw_type == 2
    surf(Vy, Vx, drReF)%3次元表面プロット
end
colorbar('Ticks',levels,...
         'TickLabels',{'10%','20%','30%','40%','50%','60%','70%','80%','90%'})
hold on
for i = 1:numTheta*numLambda
    plot(plot_point(i,1), plot_point(i,2), 'r+')
    hold on
end
xticks(-50:10:50)
yticks(-50:10:50)
ax = gca;
ax.FontSize = 18;
title(['R = ', num2str(r_measured),' mm (',num2str(time),' us)'],'FontSize',22);
xlabel('V_{Z} [km/s]','FontSize',22);
ylabel('V_{R} [km/s]','FontSize',22);
% saveas(gcf,[num2str(time),'us_',num2str(r_measured),'mm.png'])
hold off

