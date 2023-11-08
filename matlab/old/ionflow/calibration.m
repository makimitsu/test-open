function [] = calibration(NX,gas)
%�Z��CH/�K�X��(Ar:1,H:2)

%�p�����[�^���`
run define/parameter.m
NofCH = 1;

%�Z���t�@�C���ǂݍ���
cal_filename = '/Volumes/experiment/results/Doppler/Andor/IDSP/221114/Xe_96120_29to32.asc';%ICCD�t�@�C����
cal_data = importdata(cal_filename);
center_file = '221114_Xe_96120_calibration.txt';

%���f�[�^�v���b�g�p
% conX = cal_data(:,1);
% conY = cal_data(:,1);
% contour(conX,conY,cal_data(:,2:1025))

cal_X = cal_data(:,1);%X(�s�N�Z��)�����`
[LofX,LofL]=size(cal_data);%X���̒������擾
cal_L = zeros(LofX,NofCH);%L(�g��)�����`
% px2nm = zeros(NofCH,1);%nm/pixel
lambda = [lambda1 lambda2];

%�c�ʒu����
% cal_data_X = cal_data(500:700,2:1025);
% mean_cal = mean(cal_data_X(1:10,1:10),'all');
% cal_data_X = cal_data_X - mean_cal;
% spectrum_X = sum(cal_data(500:700,2:1025),1);%cal_data�̔g�������ϕ��l -> ch�ʒu�Ƀs�[�N
% mean1 = mean(spectrum_X(800:1000,1))
% Y1 = spectrum_X(:,1)-mean1;
% f = fit(cal_X,Y1,'gauss2')
% plot(cal_X,Y1)

center = importdata(center_file);%���S���W���擾
centerY = center(:,2);%�`�����l���Ή����SY���W

for i = 1:NofCH
    pixel = [center(i,3) center(i,4)];
    p = polyfit(pixel,lambda,1);
    px2nm(i,1) = p(1);
    cal_L(:,i) = polyval(p,cal_X);
end

spectrum_X=zeros(LofX,NofCH);%data1�̕������ʂ�����

%��2�s�[�N���o�p
% spectrum_X=zeros(400,NofCH);
% cal_data = cal_data(1:400,:);
% cal_X = cal_X(1:400,:);

spectrum_X(:,NX) = ...
    sum(cal_data(:,round(centerY(NX,1)-width):round(centerY(NX,1)+width)),2);
mean1 = mean(spectrum_X(100:200,NX));
Y1 = spectrum_X(:,NX)-mean1;
f = fit(cal_X,Y1,'gauss2')
plot(cal_X,Y1)
