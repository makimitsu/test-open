function [] = plot_temp_at_t(t,NR,Tz,Tr)
%温度の一次元分布をプロット
x = linspace(0.125,0.275,NR);
xx = 0.125:0.00125:0.275;
err = 15*ones(1,NR);
%異常な温度を0に仮置き
% for i = 1:NZ
%     for j = 1:NR
%         if Ti(i,j)>300
%            Ti(i,j) = 0;
%         end
%     end
% end
h = zeros(1,4);
figure('Position',[600 150 800 600])
yy = pchip(x,Tz(2,:),xx);
h(1) = errorbar(x, Tz(2,:),err,'b*','MarkerSize', 15,'LineWidth', 3,'DisplayName','z=+2.1cm');
hold on
h(3) = plot(xx,yy,'--b','LineWidth', 3);
hold on

yy = pchip(x,Tz(1,:),xx);
h(2) = errorbar(x, Tz(1,:),err,'ro','MarkerSize', 15,'LineWidth', 3,'DisplayName','z=-2.1cm');
hold on
h(4) = plot(xx,yy,'-r','LineWidth', 3);

title(['Tz/',num2str(t),' \mus'])
xlabel('r [m]')
ylabel('Ti [eV]')
xlim([0.125 0.275])
xticks(0.125:0.025:0.275)
yticks(0:50:200)
ax = gca;
ax.FontSize = 20;
legend(h(1:2));
legend('off')
hold off

figure('Position',[300 50 800 600])
yy = pchip(x,Tr(2,:),xx);
h(1) = errorbar(x, Tr(2,:),err,'b*','MarkerSize', 15,'LineWidth', 3,'DisplayName','z=+2.1cm');
hold on
h(3) = plot(xx,yy,'--b','LineWidth', 3);
hold on
yy = pchip(x,Tr(1,:),xx);
h(2) = errorbar(x, Tr(1,:),err,'ro','MarkerSize', 15,'LineWidth', 3,'DisplayName','z=-2.1cm');
hold on
h(4) = plot(xx,yy,'-r','LineWidth', 3);
title(['Tr/',num2str(t),' \mus'])
xlabel('r [m]')
ylabel('Ti [eV]')
xlim([0.125 0.275])
xticks(0.125:0.025:0.275)
yticks(0:50:500)
ax = gca;
ax.FontSize = 20;
legend(h(1:2));
legend('off')
hold off
