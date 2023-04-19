function [] = calibration(NX)
%NX=�Z��CH

%�p�����[�^���`
run define/parameter.m

center_file = 'Xe_96120_calibration_new.txt';%���S�f�[�^�t�@�C����
filename1 = '/Users/Ryo/Doppler/211017/ascii/Xe_96120.asc';%ICCD�t�@�C����

center = importdata(center_file);%���S���W���擾
centerY = center(:,2);%�`�����l���Ή����SY���W
data1 = importdata(filename1);
data1 = data1(100:600,:);

X1 = data1(:,1);%X1(�s�N�Z��)�����`
[LofX1,n1]=size(data1);%X1���̒������擾
L1 = zeros(LofX1,NofCH);%L1(�g��)�����`
px2nm = zeros(NofCH,1);%nm/pixel

lambda = [lambda1 lambda2];

for i = 1:NofCH
    pixel = [center(i,3) center(i,4)];
    p = polyfit(pixel,lambda,1);
    px2nm(i,1) = p(1);
    L1(:,i) = polyval(p,X1);
end

spectrum1=zeros(LofX1,NofCH);%data1�̕������ʂ�����
spectrum1(:,NX) = ...
    sum(data1(:,round(centerY(NX,1)-width):round(centerY(NX,1)+width)),2);
mean1 = mean(spectrum1(200:300,NX));
Y1 = spectrum1(:,NX)-mean1;
f = fit(X1,Y1,'gauss2');
coef=coeffvalues(f);
plot(X1,Y1)
fprintf('%.1f\t%.1f\t%.3f\n',coef(2), coef(5), coef(3));
