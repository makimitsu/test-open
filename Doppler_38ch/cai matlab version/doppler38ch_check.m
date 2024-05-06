clear all
% doppler_38ch_folder = '/Volumes/experiment/results/Doppler/Photron/230915/';
doppler_38ch_folder = '/Users/yunhancai/Google Drive/Data/Doppler/Photron/230914/';
doppler_38ch_name = 'H_dial48635_shot41_300kfps_3us_gain25';
doppler_38ch_path = strcat(doppler_38ch_folder,doppler_38ch_name,'/',doppler_38ch_name,'.tif');
% doppler_38ch_path = '/Users/yunhancai/Downloads/Doppler/Photron/shot60_Hbeta_11mTorr_gain25_400kfps_250frame_z=0.tif';
info = imfinfo(doppler_38ch_path);
x = 1:250;
max_intensity = zeros(size(x));
row = zeros(64,length(x));
for i = x
    thisPage = imread(doppler_38ch_path, i);
    max_intensity(i-x(1)+1) = max(thisPage,[],"all");
%     row(:,i-x(1)+1) = thisPage(:,24);
end

figure
subplot(3,1,1)
plot(x,max_intensity);
ylim([0 4095])
subplot(3,1,2)
contourf(imread(doppler_38ch_path, 152),50,'LineStyle','none');
% subplot(3,1,3)
% contourf(row,100,'LineStyle','none');