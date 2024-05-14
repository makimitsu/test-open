%マルチスリット(17CH、49CH、81CH)の横間隔を求める

clear all
close all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

date = 20240313;
ch = [17,49,81];
x_min = 1;
x_max = 1024;
y_min = 450;
y_max = 510;
idx_ax = linspace(1,1024,1024);
lambda_center = zeros(size(ch,2),1);%フィルタ中心波長
lambda_10_left = zeros(size(ch,2),1);%フィルタ10%左端
lambda_10_right = zeros(size(ch,2),1);%フィルタ10%右端
lambda_50_left = zeros(size(ch,2),1);%フィルタ半値左端
lambda_50_right = zeros(size(ch,2),1);%フィルタ半値右端
lambda_90_left = zeros(size(ch,2),1);%フィルタ90%左端
lambda_90_right = zeros(size(ch,2),1);%フィルタ90%右端
spectra = zeros(1024,size(ch,2));
figure
legStr = "";
hs = char.empty;
for i=1:size(ch,2)
    filename = [pathname.TE,'/Andor/',num2str(date),'/CVI_filter_', num2str(ch(i)), 'CH.asc'];
    data = importdata(filename);
    data = data(:,2:end);
    offset = mean(data(1:5,1:5),"all");
    data = data - offset;
    transmit = sum(data(:,y_min:y_max),2);
    [spectrum_max,~] = max(movmean(transmit,3));
    transmit = transmit/spectrum_max;%規格化
    spectra(:,i) = transmit;
    spectrum_movmean = movmean(transmit,5);
    lambda_10_left(i) = idx_ax(find(spectrum_movmean>0.1,1,'last'));
    lambda_10_right(i) = idx_ax(find(spectrum_movmean>0.1,1,'first'));
    lambda_50_left(i) = idx_ax(find(spectrum_movmean>0.5,1,'last'));
    lambda_50_right(i) = idx_ax(find(spectrum_movmean>0.5,1,'first'));
    lambda_90_left(i) = idx_ax(find(spectrum_movmean>0.9,1,'last'));
    lambda_90_right(i) = idx_ax(find(spectrum_movmean>0.9,1,'first'));
    % contour(data(min_x:max_x,min_y:max_y))
    h = plot(idx_ax,transmit,'LineWidth',2);
    [hs,legStr] = make_legend(hs,h,legStr,num2str(ch(i))+"°");
    hold on
    lambda_center(i) = (lambda_50_left(i) + lambda_50_right(i))/2;
end
newcolors = ["#005AFF" "#03AF7A" "#F6AA00" "#FF4B00"];
colororder(newcolors)
% xline(529.05,'LineWidth',3)%C VI 529.05nm
xlabel('Wavelength [nm]')
ylabel('Transmittance')
legend(hs,legStr)
ylim([0 1.2])
ax = gca;
ax.FontSize = 18;
hold off
