%%% plot raw ESP data %%%
% close all

run define_path.m

% %直接入力の場合
date = 230817;%計測日
n_ch = 21;%静電プローブCH数
ch = linspace(1,n_ch,n_ch);
trange = 450:0.1:500;%【input】計算時間範囲(0.1刻み)
shot = 16;%shot番号(ログから読み込む)

res_ratio = 50;%静電プローブ分圧比

ng_ch = [4 6 16];%死んだCH

phi = zeros(numel(trange),n_ch);
filename = sprintf("%s%03d%s",[pathname.ESP '/' num2str(date) '/ES_' num2str(date)], shot, '.csv');
ESPdata = readmatrix(filename,'Range',sprintf('B%d:V%d',trange(1)*10+2,trange(end)*10+2));
phi(:,:) = ESPdata;
if ng_ch
    phi(:,ng_ch) = [];%死んだCHを除去
    ch(ng_ch) = [];%死んだCHを除去
end
% phi = fliplr(phi);%(CH1のZ座標)>(CH2のZ座標)>...のため、列を反転
phi = phi.*res_ratio;%分圧比を掛ける
[~,raw_phi] = size(phi);

figure('Position', [0 300 700 500],'visible','on');
for i = 1:raw_phi
    plot(trange,phi(:,i))
    hold on
end
title('Raw ESP data')
xlabel('Time [us]')
ylabel('Floating Potential [V]')
% xlim([460 470])
ylim([-150 150])
legendStrings = "CH" + string(ch);
legend(legendStrings)

%Low pass filtered
figure('Position', [800 300 700 500],'visible','on');
for i = 1:raw_phi
    plot(trange,lowpass(phi(:,i),1E-4))
    hold on
end
title('Low pass filtered ESP data')
xlabel('Time [us]')
ylabel('Floating Potential [V]')
% xlim([460 470])
ylim([-150 150])
legendStrings = "CH" + string(ch);
legend(legendStrings)
