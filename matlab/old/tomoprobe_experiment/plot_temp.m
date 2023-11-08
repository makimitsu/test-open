function [] = plot_temp(NT,NR,Tz,Tr,Sz,Sr)
presetmark = ['o' 'x' '^'];
presetcolor =  [...
%     0.6350    0.0780    0.1840; % chacole
    1    0    0; % ��
%     1.            0.65,         0.        ; % ��
    0    0.69    0.48; % ��
%     0.3    0.77    1; % ���F
    0    0.4470    0.7410;    % ��
    0.5    0    0.5; % ��

    0             0             0]; %black

% presetcolor =  [0    0.4470    0.7410;    % ��
%     0.8500    0.3250    0.0980; % ��
%     0.9290    0.6940    0.1250; % �R���F
%     0.4940    0.1840    0.5560; % ��
%     0.4660    0.6740    0.1880; % ��
%     1.            0.5,         0.        ; % Orange
%     0.6350    0.0780    0.1840; % chacole
%     0             0             0]; %black

%���x�̈ꎟ�����z���v���b�g
x = linspace(0.125,0.100+0.025*NR,NR);
xx = 0.125:0.00125:0.100+0.025*NR;

NT_sub = 3;
%Tz
T_Z0 = 22.78;%�I�t�Z�b�g
h1 = zeros(2*NT_sub,1);
leg1 = zeros(NT_sub,1);
figure('Position',[600 150 800 600])
for t = 1:NT_sub
    leg1(t) = 2*t-1;
    h1(2*t-1) = errorbar(Tz(2*t-1,:)-T_Z0,x,Sz(2*t-1,:),'horizontal',presetmark(t),'Color',presetcolor(t,:),'MarkerSize', 15,'LineWidth', 3,'DisplayName',['t=',num2str(461+4*t),'\mus']');
    hold on
    yy = pchip(x,Tz(2*t-1,:),xx);
    h1(2*t) = plot(yy-T_Z0,xx,'--','Color',presetcolor(t,:),'LineWidth', 3);
    hold on
end

% title('Tz')
xlabel('T_{Z} - T_{Z0} [eV]')
ylabel('R [m]')
ylim([0.125 0.100+0.025*NR])
yticks(0.125:0.025:0.100+0.025*NR)
% yticks(0:50:200)
ax = gca;
ax.FontSize = 30;
legend(h1([leg1]));
% legend('off')
hold off

%Tr
T_R0 = 70.87;%�I�t�Z�b�g
h2 = zeros(2*NT_sub,1);
leg2 = zeros(NT_sub,1);
figure('Position',[300 50 800 600])
for t = 1:NT_sub
    leg2(t) = 2*t-1;
    h2(2*t-1) = errorbar(Tr(2*t-1,:)-T_R0,x,Sr(2*t-1,:),'horizontal',presetmark(t),'Color',presetcolor(t,:),'MarkerSize', 15,'LineWidth', 3,'DisplayName',['t=',num2str(461+4*t),'\mus']');
    hold on
    yy = pchip(x,Tr(2*t-1,:),xx);
    h2(2*t) = plot(yy-T_R0,xx,'--','Color',presetcolor(t,:),'LineWidth', 3);
    hold on
end

% title('Tr')
xlabel('T_{R} - T_{R0} [eV]')
ylabel('R [m]')
ylim([0.125 0.100+0.025*NR])
yticks(0.125:0.025:0.100+0.025*NR)
% yticks(0:50:500)
ax = gca;
ax.FontSize = 30;
legend(h2([leg2]));
% legend('off')
hold off
