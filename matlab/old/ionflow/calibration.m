function [] = calibration(NX,gas)
%校正CH/ガス種(Ar:1,H:2)

%パラメータを定義
run define/parameter.m
NofCH = 1;

%校正ファイル読み込み
cal_filename = '/Volumes/experiment/results/Doppler/Andor/IDSP/221114/Xe_96120_29to32.asc';%ICCDファイル名
cal_data = importdata(cal_filename);
center_file = '221114_Xe_96120_calibration.txt';

%生データプロット用
% conX = cal_data(:,1);
% conY = cal_data(:,1);
% contour(conX,conY,cal_data(:,2:1025))

cal_X = cal_data(:,1);%X(ピクセル)軸を定義
[LofX,LofL]=size(cal_data);%X軸の長さを取得
cal_L = zeros(LofX,NofCH);%L(波長)軸を定義
% px2nm = zeros(NofCH,1);%nm/pixel
lambda = [lambda1 lambda2];

%縦位置特定
% cal_data_X = cal_data(500:700,2:1025);
% mean_cal = mean(cal_data_X(1:10,1:10),'all');
% cal_data_X = cal_data_X - mean_cal;
% spectrum_X = sum(cal_data(500:700,2:1025),1);%cal_dataの波長方向積分値 -> ch位置にピーク
% mean1 = mean(spectrum_X(800:1000,1))
% Y1 = spectrum_X(:,1)-mean1;
% f = fit(cal_X,Y1,'gauss2')
% plot(cal_X,Y1)

center = importdata(center_file);%中心座標を取得
centerY = center(:,2);%チャンネル対応中心Y座標

for i = 1:NofCH
    pixel = [center(i,3) center(i,4)];
    p = polyfit(pixel,lambda,1);
    px2nm(i,1) = p(1);
    cal_L(:,i) = polyval(p,cal_X);
end

spectrum_X=zeros(LofX,NofCH);%data1の分光結果を入れる

%第2ピーク検出用
% spectrum_X=zeros(400,NofCH);
% cal_data = cal_data(1:400,:);
% cal_X = cal_X(1:400,:);

spectrum_X(:,NX) = ...
    sum(cal_data(:,round(centerY(NX,1)-width):round(centerY(NX,1)+width)),2);
mean1 = mean(spectrum_X(100:200,NX));
Y1 = spectrum_X(:,NX)-mean1;
f = fit(cal_X,Y1,'gauss2')
plot(cal_X,Y1)
