function [] = plot_velocity(NT,NR,Vmean_z2,Vmean_r2,Sz2,Sr2)
NT_sub = 3;%プロット数
if NT_sub == 6
    presetmark = ['o' 'd' 'x' '*' '^' 's'];
    presetcolor =  [...
        1    0    0; % 赤
        1.            0.65,         0.        ; % 橙
        0    0.69    0.48; % 緑
        0.3    0.77    1; % 水色
        0    0.4470    0.7410;    % 青
        0.5    0    0.5; % 紫
        0.6350    0.0780    0.1840; % chacole
        0             0             0]; %black
elseif NT_sub == 3
    presetmark = ['o' 'x' '^'];
    presetcolor =  [...
        1    0    0; % 赤
%         1.            0.65,         0.        ; % 橙
        0    0.69    0.48; % 緑
%         0.3    0.77    1; % 水色
        0    0.4470    0.7410;    % 青
        0.5    0    0.5; % 紫
        0.6350    0.0780    0.1840; % chacole
        0             0             0]; %black
end

% presetcolor =  [0    0.4470    0.7410;    % 青
%     0.8500    0.3250    0.0980; % 赤
%     0.9290    0.6940    0.1250; % 山吹色
%     0.4940    0.1840    0.5560; % 紫
%     0.4660    0.6740    0.1880; % 緑
%     1.            0.5,         0.        ; % Orange
%     0.6350    0.0780    0.1840; % chacole
%     0             0             0]; %black

%一次元分布をプロット
x = linspace(0.125,0.100+0.025*NR,NR);
xx = 0.125:0.00125:0.100+0.025*NR;

%Vmean_z
h1 = zeros(2*NT_sub,1);
leg1 = zeros(NT_sub,1);
figure('Position',[600 150 800 600])
for t = 1:NT_sub
    leg1(t) = 2*t-1;
    if NT_sub == 6
        h1(2*t-1) = errorbar(Vmean_z2(t,:),x,Sz2(t,:),'horizontal',presetmark(t),'Color',presetcolor(t,:),'MarkerSize', 15,'LineWidth', 3,'DisplayName',['t=',num2str(463+2*t),'\mus']');
        hold on
        yy = pchip(x,Vmean_z2(t,:),xx);
        h1(2*t) = plot(yy,xx,'--','Color',presetcolor(t,:),'LineWidth', 3);
    elseif NT_sub == 3
        h1(2*t-1) = errorbar(Vmean_z2(2*t-1,:),x,Sz2(2*t-1,:),'horizontal',presetmark(t),'Color',presetcolor(t,:),'MarkerSize', 15,'LineWidth', 3,'DisplayName',['t=',num2str(461+4*t),'\mus']');
        hold on
        yy = pchip(x,Vmean_z2(2*t-1,:),xx);
        h1(2*t) = plot(yy,xx,'--','Color',presetcolor(t,:),'LineWidth', 3);
        hold on
    end
end

xlabel('m_{i}V^{2}_{0Z}/2k_{B} [eV]')
ylabel('R [m]')
xlim([0 2])
ylim([0.125 0.100+0.025*NR])
yticks(0.125:0.025:0.100+0.025*NR)
ax = gca;
ax.FontSize = 20;
legend(h1([leg1]));
% legend('off')
hold off

%Vmean_r
h2 = zeros(2*NT_sub,1);
leg2 = zeros(NT_sub,1);
figure('Position',[300 50 800 600])
for t = 1:NT_sub
    leg2(t) = 2*t-1;
    if NT_sub == 6
        h2(2*t-1) = errorbar(Vmean_r2(t,:),x,Sr2(t,:),'horizontal',presetmark(t),'Color',presetcolor(t,:),'MarkerSize', 15,'LineWidth', 3,'DisplayName',['t=',num2str(463+2*t),'\mus']');
        hold on
        yy = pchip(x,Vmean_r2(t,:),xx);
        h2(2*t) = plot(yy,xx,'--','Color',presetcolor(t,:),'LineWidth', 3);
    elseif NT_sub == 3
        h2(2*t-1) = errorbar(Vmean_r2(2*t-1,:),x,Sr2(2*t-1,:),'horizontal',presetmark(t),'Color',presetcolor(t,:),'MarkerSize', 15,'LineWidth', 3,'DisplayName',['t=',num2str(461+4*t),'\mus']');
        hold on
        yy = pchip(x,Vmean_r2(2*t-1,:),xx);
        h2(2*t) = plot(yy,xx,'--','Color',presetcolor(t,:),'LineWidth', 3);
        hold on
    end
end

xlabel('m_{i}V^{2}_{0R}/2k_{B} [eV]')
ylabel('R [m]')
xlim([0 60])
ylim([0.125 0.100+0.025*NR])
yticks(0.125:0.025:0.100+0.025*NR)
ax = gca;
ax.FontSize = 30;
legend(h2([leg2]));
% legend('off')
hold off

%Total flow energy
h3 = zeros(2*NT_sub,1);
leg3 = zeros(NT_sub,1);
figure('Position',[300 50 800 600])
for t = 1:NT_sub
    leg3(t) = 2*t-1;
    if NT_sub == 6
        h3(2*t-1) = errorbar(Vmean_z2(t,:)+Vmean_r2(t,:),x,Sz2(t,:)+Sr2(t,:),'horizontal',presetmark(t),'Color',presetcolor(t,:),'MarkerSize', 15,'LineWidth', 3,'DisplayName',['t=',num2str(463+2*t),'\mus']');
        hold on
        yy = pchip(x,Vmean_z2(t,:)+Vmean_r2(t,:),xx);
        h3(2*t) = plot(yy,xx,'--','Color',presetcolor(t,:),'LineWidth', 3);
    elseif NT_sub == 3
        h3(2*t-1) = errorbar(Vmean_z2(2*t-1,:)+Vmean_r2(2*t-1,:),x,Sz2(2*t-1,:)+Sr2(2*t-1,:),'horizontal',presetmark(t),'Color',presetcolor(t,:),'MarkerSize', 15,'LineWidth', 3,'DisplayName',['t=',num2str(461+4*t),'\mus']');
        hold on
        yy = pchip(x,Vmean_z2(2*t-1,:)+Vmean_r2(2*t-1,:),xx);
        h3(2*t) = plot(yy,xx,'--','Color',presetcolor(t,:),'LineWidth', 3);
        hold on
    end
end

xlabel('m_{i}V^{2}_{0}/2k_{B} [eV]')
ylabel('R [m]')
xlim([0 60])
ylim([0.125 0.100+0.025*NR])
yticks(0.125:0.025:0.100+0.025*NR)
ax = gca;
ax.FontSize = 30;
legend(h3([leg3]));
% legend('off')
hold off
