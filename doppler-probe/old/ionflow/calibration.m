function [] = calibration(NX)
%NX=校正CH

%パラメータを定義
run define/parameter.m

center_file = 'Xe_96120_calibration_new.txt';%中心データファイル名
filename1 = '/Users/Ryo/Doppler/211017/ascii/Xe_96120.asc';%ICCDファイル名

center = importdata(center_file);%中心座標を取得
centerY = center(:,2);%チャンネル対応中心Y座標
data1 = importdata(filename1);
data1 = data1(100:600,:);

X1 = data1(:,1);%X1(ピクセル)軸を定義
[LofX1,n1]=size(data1);%X1軸の長さを取得
L1 = zeros(LofX1,NofCH);%L1(波長)軸を定義
px2nm = zeros(NofCH,1);%nm/pixel

lambda = [lambda1 lambda2];

for i = 1:NofCH
    pixel = [center(i,3) center(i,4)];
    p = polyfit(pixel,lambda,1);
    px2nm(i,1) = p(1);
    L1(:,i) = polyval(p,X1);
end

spectrum1=zeros(LofX1,NofCH);%data1の分光結果を入れる
spectrum1(:,NX) = ...
    sum(data1(:,round(centerY(NX,1)-width):round(centerY(NX,1)+width)),2);
mean1 = mean(spectrum1(200:300,NX));
Y1 = spectrum1(:,NX)-mean1;
f = fit(X1,Y1,'gauss2');
coef=coeffvalues(f);
plot(X1,Y1)
fprintf('%.1f\t%.1f\t%.3f\n',coef(2), coef(5), coef(3));
