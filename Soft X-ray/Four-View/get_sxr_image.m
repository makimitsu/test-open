function [VectorImage1,VectorImage2] = get_sxr_image(date,number,ProjectionNumber,TIFImagePath,doFilter)
% この関数は，一枚のTIF画像から再構成直前のベクトル画像を生成します．
% 事前フィルタ処理や，時系列処理はここで行われます

% --以下竹田さんメモ
% 画像を切り取ってベクトル形式で出力
% 全時間帯の画像を全て出力するかどうか
% 使う分だけでいい？
% 4つのベクトルを出力するか
% 画像は校正用イメージをもとに切り取る感じで
% 入射角度補正をかけるかどうか
% 校正用データをもとに半径決めてその半径内で強度の補正をつけると小さい円の外周部が過剰に大きくなる
% 半径は固定？それとも可変？
% 4つのファイバーそれぞれの半径が同じくらいであれば4つの異なる半径を採用
% 大体全部同じようなら小さめで採用
% そもそも今回は断面全体は使わない予定なので問題なし？

doCheck = true;

% 生画像の取得
TIFImage = imread(TIFImagePath);

% 非線形フィルターをかける（必要があれば）
if doFilter
    figure;imagesc(TIFImage);
    [TIFImage,~] = imnlmfilt(TIFImage,'SearchWindowSize',91,'ComparisonWindowSize',15);
    figure;imagesc(TIFImage);
end

% ファイバーの中心を特定し，時系列順に並べます．
% 校正用画像PositionCheck.tifから円検出を行い，時系列データをcentersに格納します
PositionCheckImagePath = strcat('G:\My Drive\X-ray\Data\TIF\',num2str(date),'/PositionCheck.tif');
PositionCheckImage = imread(PositionCheckImagePath);
[centers,radius]=find_fibers2(PositionCheckImage);
RadiusGeneral = round(mean(radius));
centers = round(centers);

% それぞれのファイバーの中心を格納するための配列を定義（ここまでFindFibersでやれば良くない？）
Center = zeros(2,8,2);
Center(1,:,:) = centers([1,3,6,8,9,11,14,16],:);
Center(2,:,:) = centers([2,4,5,7,10,12,13,15],:);

% 切り取った画像を格納するための配列
timeRapsImage = zeros(2,4,2*RadiusGeneral,2*RadiusGeneral);

% % 各時間帯ごとに画像を切り取って格納
% for i = 1:8
%     UpRangeV = Center(1,i,1)-IW+1:Center(1,i,1)+IW;
%     UpRangeH = Center(1,i,2)-IW+1:Center(1,i,2)+IW;
%     DownRangeV = Center(2,i,1)-IW+1:Center(2,i,1)+IW;
%     DownRangeH = Center(2,i,2)-IW+1:Center(2,i,2)+IW;
%     TimeRapsImage(1,i,:,:) = rawImage(UpRangeV,UpRangeH,1);
%     TimeRapsImage(2,i,:,:) = rawImage(DownRangeV,DownRangeH,1);
%     TimeRapsImage(1,i,:,:) = rot90(TimeRapsImage(1,i,:,:),2);
%     TimeRapsImage(2,i,:,:) = rot90(TimeRapsImage(2,i,:,:),2);
% end

% バックグラウンドノイズのデータを取得
BackGround = cast(TIFImage(650-RadiusGeneral+1:650+RadiusGeneral,120-RadiusGeneral+1:120+RadiusGeneral,1),'double');
BackGround = ones(size(BackGround))*mean(BackGround,'all');

% 切り取った画像のうち実際に使う部分（ファイバー部分）を切り出し
k = FindCircle(ProjectionNumber/2);
Image_vectors = zeros(2,8,numel(k));

if doCheck
    f1=figure;
    f1.Position = [900,200,500,500];
    % f2=figure;
    % f2.Position = [900,200,500,500];
end

% 明るさに関する校正係数（matファイル）を取得します．とりあえず今はお休み
% CalibrationPath = strcat('G:\My Drive\X-ray\Data\TIF\',num2str(date),'/CalibrationFactor.mat');
% if exist(CalibrationPath,'file')
%     load(CalibrationPath,'CalibrationFactor');
% else
%     CalibrationFactor = get_calibration_factor(date,projectionNumber);
% end

% 再構成計算に使用可能なサイズに変換するための比を取得
DownscaleRate = ProjectionNumber/(RadiusGeneral*2);

% 画像データの行列化
for i=1:4
    % 画像切り取りの縦方向・横方向範囲を指定
    verticalRange1 = Center(1,i,1)-RadiusGeneral+1:Center(1,i,1)+RadiusGeneral;
    horizontalRange1 = Center(1,i,2)-RadiusGeneral+1:Center(1,i,2)+RadiusGeneral;
    verticalRange2 = Center(2,i,1)-RadiusGeneral+1:Center(2,i,1)+RadiusGeneral;
    horizontalRange2 = Center(2,i,2)-RadiusGeneral+1:Center(2,i,2)+RadiusGeneral;
    % 生画像を切り取ってファイバー数×フレーム数×IW×IWの配列に格納
    timeRapsImage(1,i,:,:) = TIFImage(verticalRange1,horizontalRange1,1);
    timeRapsImage(2,i,:,:) = TIFImage(verticalRange2,horizontalRange2,1);
    % timeRapsImage(1,i,:,:) = rot90(timeRapsImage(1,i,:,:),2);
    % timeRapsImage(2,i,:,:) = rot90(timeRapsImage(2,i,:,:),2);

    % バックグラウンドノイズ成分を差し引く
    singleImage1 = squeeze(timeRapsImage(1,i,:,:))-BackGround;
    singleImage2 = squeeze(timeRapsImage(2,i,:,:))-BackGround;
    % GrayImage1 = fliplr(SingleImage1(:,:));
    % GrayImage2 = fliplr(SingleImage2(:,:));
    grayImage1 = flipud(singleImage1(:,:));
    grayImage2 = flipud(singleImage2(:,:));
    grayImage1(grayImage1<0) = 0;
    grayImage2(grayImage2<0) = 0;
    % 再構成計算に使用可能なサイズに変換
    roughImage1 = imresize(grayImage1, DownscaleRate, 'nearest');
    roughImage1 = cast(roughImage1,'double');
    roughImage2 = imresize(grayImage2, DownscaleRate, 'nearest');
    roughImage2 = cast(roughImage2,'double');

    % 校正データを用いた明るさの修正（カメラの影の補正とか）や入射角度の補正をここでやりたい
    % figure;imagesc(RoughImage1);
    % RoughImage1(k) = RoughImage1(k).*squeeze(CalibrationFactor(1,i,:));
    % figure;imagesc(RoughImage1);

    % 切り出した画像を表示
    if doCheck
        figure(f1);
        i_str = num2str(i);
        title1 = strcat('1,',i_str);
        title2 = strcat('2,',i_str);
        subplot(4,4,2*(i-1)+1);imagesc(roughImage1);title(title1);
%         caxis([50,60]);
        subplot(4,4,2*(i-1)+2);imagesc(roughImage2);title(title2);
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
    Image_vectors(1,i,:) = roughImage1(k);
    Image_vectors(2,i,:) = roughImage2(k);
    % VectorImages(1,i,:) = RoughImage1(k).*squeeze(CalibrationFactor(1,i,:));
    % VectorImages(2,i,:) = RoughImage2(k).*squeeze(CalibrationFactor(2,i,:));
end

VectorImages1 = squeeze(Image_vectors(1,:,:));
VectorImages2 = squeeze(Image_vectors(2,:,:));
VectorImage1 = VectorImages1(number,:);
VectorImage2 = VectorImages2(number,:);
end


function InCircleIndex = FindCircle(Radius)
Square = zeros(2*Radius);
for i = 1:2*Radius
    for j = 1:2*Radius
        Square(i,j) = sqrt((Radius-i+0.5)^2+(j-Radius-0.5)^2);
    end
end
InCircleIndex = find(Square<Radius);
end