function [VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename,filter)
% date = 230119;
% number = 5;
% SXRfilename = '/Users/shinjirotakeda/Library/CloudStorage/OneDrive-TheUniversityofTokyo/Documents/SXR_Images/230119/shot024.tif';
% filter = false;

% % 画像を切り取る
N_projection = 80;
% VectorImages = CutImage(date,shot,N_projection/230,false);
% CheckFlag = true;
CheckFlag = false;
VectorImages = CutImage(date,N_projection,CheckFlag,SXRfilename,filter);
VectorImages1 = squeeze(VectorImages(1,:,:));
VectorImages2 = squeeze(VectorImages(2,:,:));
VectorImage1 = VectorImages1(number,:);
VectorImage2 = VectorImages2(number,:);
end

function VectorImages = CutImage(date,N_projection,check_flag,SXRfilename,filter)

RawImage = imread(SXRfilename);

if filter
    figure;imagesc(RawImage);
    [RawImage,~] = imnlmfilt(RawImage,'SearchWindowSize',91,'ComparisonWindowSize',15);
    figure;imagesc(RawImage);
% else
%     RawImage = RawImage;
end

% IW = 115; %Image Width
% IW = 120; %Image Width

PositionPath = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/',num2str(date),'/PositionCheck.tif');
CalibrationImage = imread(PositionPath);
if date == 230315
    centers = 1.0e3*[...
        0.2815,1.5525;0.4834,1.5468;0.8174,1.5699;1.0084,1.5579;...
        1.0215,1.1871;0.8132,1.1972;0.4748,1.1887;0.2770,1.1942;...
        0.2766,0.8066;0.4750,0.7979;0.8203,0.8339;1.0226,0.8280;...
        1.0169,0.4307;0.8187,0.4377;0.4819,0.4108;0.2794,0.4245];
    IW = 95;
elseif date == 230119
    centers = 1.0e3*[...
        0.2815,1.5525;0.5234,1.5168;0.8174,1.5699;1.0584,1.5279;...
        1.0665,1.1621;0.8532,1.1972;0.5248,1.1487;0.3370,1.1542;...
        0.3366,0.7766;0.5250,0.7679;0.8753,0.7939;1.0676,0.7830;...
        1.0669,0.4107;0.8787,0.4177;0.5319,0.3808;0.3394,0.3945];
    IW = 90;
elseif date == 230316
    [centers,radii]=FindFibe0rs(CalibrationImage,[90,130]);
    IW = 90;
    centers(:,1) = centers(:,1) + 40;
    centers(:,2) = centers(:,2) - 20;
else
    [centers,radii]=FindFibers(CalibrationImage,[70,130]);
    IW = 85;
end
% 上にずらすときは1つ目の数字を増やす
% 右にずらすときは2つめの数字を増やす
% IW = round(mean(radii));
centers = round(centers);


Center = zeros(2,8,2);
TimeRapsImage = zeros(2,8,2*IW,2*IW);

Center(1,:,:) = centers([1,3,6,8,9,11,14,16],:);
Center(2,:,:) = centers([2,4,5,7,10,12,13,15],:);
% 
% Center(2,1,:) = [530,1520];
% Center(2,2,:) = [1060,1530];
% Center(2,3,:) = [1070,1160];
% Center(2,4,:) = [530,1150];
% Center(2,5,:) = [530,770];
% Center(2,6,:) = [1070,780];
% Center(2,7,:) = [1070,410];
% Center(2,8,:) = [530,380];
% Center(1,:,1) = Center(2,:,1) - 190;
% Center(1,:,2) = Center(2,:,2) + 5;

% disp(Center);
% whos RawImage
% whos TimeRapsImage

for i = 1:8
    UpRangeV = Center(1,i,1)-IW+1:Center(1,i,1)+IW;
    UpRangeH = Center(1,i,2)-IW+1:Center(1,i,2)+IW;
    DownRangeV = Center(2,i,1)-IW+1:Center(2,i,1)+IW;
    DownRangeH = Center(2,i,2)-IW+1:Center(2,i,2)+IW;
    TimeRapsImage(1,i,:,:) = RawImage(UpRangeV,UpRangeH,1);
    TimeRapsImage(2,i,:,:) = RawImage(DownRangeV,DownRangeH,1);
    TimeRapsImage(1,i,:,:) = rot90(TimeRapsImage(1,i,:,:),2);
    TimeRapsImage(2,i,:,:) = rot90(TimeRapsImage(2,i,:,:),2);
end

BackGround = cast(RawImage(1050-IW+1:1050+IW,120-IW+1:120+IW,1),'double');
BackGround = ones(size(BackGround))*mean(BackGround,'all');

k = FindCircle(N_projection/2);
VectorImages = zeros(2,8,numel(k));

if check_flag
    f1=figure;
%     f.WindowState = 'maximized';
    f1.Position = [900,200,500,500];
%     f2=figure;
% %     f.WindowState = 'maximized';
%     f2.Position = [900,200,500,500];
end

CalibrationPath = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/',num2str(date),'/CalibrationFactor.mat');
if exist(CalibrationPath,'file')
    load(CalibrationPath,'CalibrationFactor');
else
    CalibrationFactor = clc_CalibrationFactor(date,N_projection);
end

resolution = N_projection/(IW*2);

for i=1:8
    SingleImage1 = squeeze(TimeRapsImage(1,i,:,:))-BackGround;
    SingleImage2 = squeeze(TimeRapsImage(2,i,:,:))-BackGround;
    GrayImage1 = fliplr(SingleImage1(:,:));
    GrayImage2 = fliplr(SingleImage2(:,:));
    GrayImage1(GrayImage1<0) = 0;
    GrayImage2(GrayImage2<0) = 0;
    RoughImage1 = imresize(GrayImage1, resolution, 'nearest');
    RoughImage1 = cast(RoughImage1,'double');
    RoughImage2 = imresize(GrayImage2, resolution, 'nearest');
    RoughImage2 = cast(RoughImage2,'double');
    
%     whos RoughImage1
%     whos SD
%     
%     RoughImage1 = RoughImage1./fliplr(SD./100);
%     RoughImage2 = RoughImage2./fliplr(SD./100);
    
    if check_flag
        figure(f1);
        i_str = num2str(i);
        title1 = strcat('1,',i_str);
        title2 = strcat('2,',i_str);
        subplot(4,4,2*(i-1)+1);imagesc(RoughImage1);title(title1);
%         caxis([50,60]);
        subplot(4,4,2*(i-1)+2);imagesc(RoughImage2);title(title2);
%         caxis([50,60]);
%         if i == 8
%             colorbar;
%         end
        % figure(f2);
        % RoughCalibrated1 = RoughImage1;
        % RoughCalibrated1(k) = RoughCalibrated1(k).*squeeze(CalibrationFactor(1,i,:));
        % subplot(4,4,2*(i-1)+1);imagesc(RoughCalibrated1);title(title1);
        % RoughCalibrated2 = RoughImage2;
        % RoughCalibrated2(k) = RoughCalibrated2(k).*squeeze(CalibrationFactor(1,i,:));
        % subplot(4,4,2*(i-1)+2);imagesc(RoughCalibrated2);title(title2);
    end
    
    % figure;imagesc(RoughImage1);
    % RoughImage1(k) = RoughImage1(k).*squeeze(CalibrationFactor(1,i,:));
    % figure;imagesc(RoughImage1);

    VectorImages(1,i,:) = RoughImage1(k);
    VectorImages(2,i,:) = RoughImage2(k);
    % VectorImages(1,i,:) = RoughImage1(k).*squeeze(CalibrationFactor(1,i,:));
    % VectorImages(2,i,:) = RoughImage2(k).*squeeze(CalibrationFactor(2,i,:));
end

end

function k = FindCircle(L)
R = zeros(2*L);
for i = 1:2*L
    for j = 1:2*L
        R(i,j) = sqrt((L-i+0.5)^2+(j-L-0.5)^2);
    end
end
% figure;imagesc(R)
k = find(R<L);
end