function [imageVector1,imageVector2] = get_sxr_image_new(date,number,SXRfilename, filter)

rawImage = imread(SXRfilename);
projectionNumber = 80;
doCheck = false;

if filter
    % figure;imagesc(rawImage);
    [rawImage,~] = imnlmfilt(rawImage,'SearchWindowSize',91,'ComparisonWindowSize',15);
    % figure;imagesc(rawImage);
end

% 位置情報ファイルからファイバーの位置（＋半径）を取得
positionPath = 'fiberPositions_2view.xlsx';
positionData = readmatrix(positionPath,'Sheet',num2str(date),'Range','C2:E17');
Center = zeros(2,8,2);
for i = 1:2
    Center(i,:,:) = positionData(1+8*(i-1):8+8*(i-1),1:2);
end
IW = positionData(1,3);

Center = round(Center);

% 切り取った画像を格納するための配列
timeSeries = zeros(2,8,2*IW,2*IW);

% バックグラウンドノイズのデータを取得
backgroundImage = cast(rawImage(1:2*IW,1:2*IW,1),'double');
backgroundNoise = ones(size(backgroundImage))*mean(backgroundImage,'all');

% 切り取った画像のうち実際に使う部分（ファイバー部分）を切り出し
k = find_circle(projectionNumber/2);
imageVectors = zeros(2,8,numel(k));

if doCheck
    f1=figure;
    % f1.Position = [900,200,500,500];
    f1.Position = [200,250,1060,500];
    % f2=figure;
    % f2.Position = [900,200,500,500];
end

% % 校正用データ（matファイル）が存在する場合はそれを取得、しなければ計算
% calibrationPath = strcat(getenv('SXR_IMAGE_DIR'),'/',num2str(date),'/calibrationFactor.mat');
% if exist(calibrationPath,'file')
%     load(calibrationPath,'calibrationFactor');
%     if numel(calibrationFactor(1,1,:)) ~= numel(k)
%         calibrationFactor = get_calibration_factor(date,projectionNumber);
%     end
% else
%     calibrationFactor = get_calibration_factor(date,projectionNumber);
% end

% 再構成計算に使用可能なサイズに変換するための解像度を取得
resolution = projectionNumber/(IW*2);

% 画像データの行列化
for i=1:8
    % 画像切り取りの縦方向・横方向範囲を指定
    horizontalRange1 = Center(1,i,1)-IW+1:Center(1,i,1)+IW;
    verticalRange1 = Center(1,i,2)-IW+1:Center(1,i,2)+IW;
    horizontalRange2 = Center(2,i,1)-IW+1:Center(2,i,1)+IW;
    verticalRange2 = Center(2,i,2)-IW+1:Center(2,i,2)+IW;
    % 生画像を切り取ってファイバー数×フレーム数×IW×IWの配列に格納
    timeSeries(1,i,:,:) = rawImage(verticalRange1,horizontalRange1,1);
    timeSeries(2,i,:,:) = rawImage(verticalRange2,horizontalRange2,1);

    % バックグラウンドノイズ成分を差し引く
    singleImage1 = squeeze(timeSeries(1,i,:,:))-backgroundNoise;
    singleImage2 = squeeze(timeSeries(2,i,:,:))-backgroundNoise;
    grayImage1 = flipud(singleImage1(:,:));
    grayImage2 = flipud(singleImage2(:,:));
    grayImage1(grayImage1<0) = 0;
    grayImage2(grayImage2<0) = 0;
    % 再構成計算に使用可能なサイズに変換
    roughImage1 = imresize(grayImage1, resolution, 'nearest');
    roughImage1 = cast(roughImage1,'double');
    roughImage2 = imresize(grayImage2, resolution, 'nearest');
    roughImage2 = cast(roughImage2,'double');

    % % 校正データを用いた明るさの修正（カメラの影の補正とか）や入射角度の補正をここでやりたい
    % % figure;imagesc(roughImage1);
    % roughImage1(k) = roughImage1(k).*squeeze(calibrationFactor(1,i,:));
    % roughImage2(k) = roughImage2(k).*squeeze(calibrationFactor(2,i,:));
    % roughImage3(k) = roughImage3(k).*squeeze(calibrationFactor(3,i,:));
    % roughImage4(k) = roughImage4(k).*squeeze(calibrationFactor(4,i,:));
    % % figure;imagesc(roughImage1);

    % 切り出した画像を表示
    if doCheck
        figure(f1);
        i_str = num2str(i);
        title1 = strcat('1,',i_str);
        title2 = strcat('2,',i_str);
        subplot(2,8,2*(i-1)+1);imagesc(roughImage1);title(title1);
%         caxis([50,60]);
        subplot(2,8,2*(i-1)+2);imagesc(roughImage2);title(title2);
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

    % ベクトル化
    imageVectors(1,i,:) = roughImage1(k);
    imageVectors(2,i,:) = roughImage2(k);
end

imageVectors1 = squeeze(imageVectors(1,:,:));
imageVectors2 = squeeze(imageVectors(2,:,:));

imageVector1 = imageVectors1(number,:);
imageVector2 = imageVectors2(number,:);

end


function k = find_circle(L)
R = zeros(2*L);
for i = 1:2*L
    for j = 1:2*L
        R(i,j) = sqrt((L-i+0.5)^2+(j-L-0.5)^2);
    end
end
% figure;imagesc(R)
k = find(R<L);
end