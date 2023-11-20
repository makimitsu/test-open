function [imageVector1,imageVector2,imageVector3,imageVector4] = get_sxr_image(Date,TimeFrameNumber,N_before,TIFImagePath,ApplyFilter)
% この関数では，一枚のTIF画像からファイバの断面のみを切り出し，再構成直前のベクトル画像を生成します．
% 上下整理，フィルタ処理，時系列整理などは，ここに含まれます．

% ***********ここの数字*と*findfibers2の数字両方いじること************
NumberOfFrame = 8;

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

% 画像を切り取る
doCheck = false;

% 生画像の取得
% TIFImage = imread(TIFImagePath);
TIFImage = TIFImagePath;
% imagesc(TIFImage);

% 非線形フィルターをかける（必要があれば）これより上の階層の関数で，既にフィルタを適用しているのでここは不要
% if ApplyFilter
%     figure;imagesc(TIFImage);
%     [TIFImage,~] = imnlmfilt(TIFImage,'SearchWindowSize',91,'ComparisonWindowSize',15);
%     figure;imagesc(TIFImage);
% end

% ファイバーの位置を検索するための校正用画像を取得
PositionCheckImagePath = strcat('G:\My Drive\X-ray\Data\TIF\',num2str(Date),'/PositionCheck.tif');
PositionCheckImage = imread(PositionCheckImagePath);

% % 校正用画像からファイバーの位置（＋半径）を取得
% [centers,radii]=find_fibers(calibrationImage,[70,130]);
% IW = 85;
% % IW = round(mean(radii));
% centers = round(centers);
% % それぞれのファイバーの中心を格納するための配列を定義（ここまでFindFibersでやれば良くない？）
% Center = zeros(2,8,2);
% Center(1,:,:) = centers([1,3,6,8,9,11,14,16],:);
% Center(2,:,:) = centers([2,4,5,7,10,12,13,15],:);

% 校正用画像からファイバーの位置（＋半径）を取得
[CenterPositions,radius] = find_fibers2(PositionCheckImage);
CenterPositions = round(CenterPositions);
radii = radius*2;

% 切り取った画像を格納するための配列
Seriesed_TIF = zeros(4,NumberOfFrame,radii,radii);

% バックグラウンドノイズのデータを取得
backgroundImage = cast(TIFImage(1:radii,1:radii,1),'double'); % ここの範囲の取り方reservations
backgroundNoise = ones(size(backgroundImage))*mean(backgroundImage,'all');

% 切り取った画像のうち実際に使う部分（ファイバー部分）を切り出す準備
k = find_circle(N_before/2);
imageVectors = zeros(4,NumberOfFrame,numel(k));

if doCheck
    f1=figure;
    % f1.Position = [900,200,500,500];
    f1.Position = [200,250,1060,500];
    % f2=figure;
    % f2.Position = [900,200,500,500];
end

% 明るさ校正データ.matを取得，あるいは新規に計算します 今はお休み
% CalibrationPath = strcat(getenv('SXR_IMAGE_DIR'),'/',num2str(Date),'/CalibrationFactor.mat');
% if exist(CalibrationPath,'file')
%     load(CalibrationPath,'CalibrationFactor');
% else
%     CalibrationFactor = get_calibration_factor(Date,N_before);
% end

% 再構成計算に使用可能なサイズに変換するための解像度を取得
DownscaleRate = N_before/radii;

% 画像データの行列化
for i=1:NumberOfFrame
    % 画像切り取りの縦方向・横方向範囲を指定
    X_Range1 = CenterPositions(1,i,1)-radius+1:CenterPositions(1,i,1)+radius;
    Y_Range1 = CenterPositions(1,i,2)-radius+1:CenterPositions(1,i,2)+radius;
    X_Range2 = CenterPositions(2,i,1)-radius+1:CenterPositions(2,i,1)+radius;
    Y_Range2 = CenterPositions(2,i,2)-radius+1:CenterPositions(2,i,2)+radius;
    X_Range3 = CenterPositions(3,i,1)-radius+1:CenterPositions(3,i,1)+radius;
    Y_Range3 = CenterPositions(3,i,2)-radius+1:CenterPositions(3,i,2)+radius;
    X_Range4 = CenterPositions(4,i,1)-radius+1:CenterPositions(4,i,1)+radius;
    Y_Range4 = CenterPositions(4,i,2)-radius+1:CenterPositions(4,i,2)+radius;
    % 生画像を切り取ってファイバー数 × 時間 × radius × radiusの配列に格納
    Seriesed_TIF(1,i,:,:) = TIFImage(Y_Range1,X_Range1,1);
    Seriesed_TIF(2,i,:,:) = TIFImage(Y_Range2,X_Range2,1);
    Seriesed_TIF(3,i,:,:) = TIFImage(Y_Range3,X_Range3,1);
    Seriesed_TIF(4,i,:,:) = TIFImage(Y_Range4,X_Range4,1);

    % バックグラウンドノイズ成分を差し引く
    singleTIF1 = squeeze(Seriesed_TIF(1,i,:,:))-backgroundNoise;
    singleTIF2 = squeeze(Seriesed_TIF(2,i,:,:))-backgroundNoise;
    singleTIF3 = squeeze(Seriesed_TIF(3,i,:,:))-backgroundNoise;
    singleTIF4 = squeeze(Seriesed_TIF(4,i,:,:))-backgroundNoise;
    grayImage1 = flipud(singleTIF1(:,:));
    grayImage2 = flipud(singleTIF2(:,:));
    grayImage3 = flipud(singleTIF3(:,:));
    grayImage4 = flipud(singleTIF4(:,:));
    grayImage1(grayImage1<0) = 0;
    grayImage2(grayImage2<0) = 0;
    grayImage3(grayImage3<0) = 0;
    grayImage4(grayImage4<0) = 0;

    % 再構成計算に使用するサイズまでダウンコンバートします．
    roughImage1 = imresize(grayImage1, DownscaleRate, 'nearest');
    roughImage1 = cast(roughImage1,'double');
    roughImage2 = imresize(grayImage2, DownscaleRate, 'nearest');
    roughImage2 = cast(roughImage2,'double');
    roughImage3 = imresize(grayImage3, DownscaleRate, 'nearest');
    roughImage3 = cast(roughImage3,'double');
    roughImage4 = imresize(grayImage4, DownscaleRate, 'nearest');
    roughImage4 = cast(roughImage4,'double');

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
        title3 = strcat('3,',i_str);
        title4 = strcat('4,',i_str);
        subplot(4,8,4*(i-1)+1);imagesc(roughImage1);title(title1);
%         caxis([50,60]);
        subplot(4,8,4*(i-1)+2);imagesc(roughImage2);title(title2);
%         caxis([50,60]);
        subplot(4,8,4*(i-1)+3);imagesc(roughImage3);title(title3);
        subplot(4,8,4*(i-1)+4);imagesc(roughImage4);title(title4);
        set(gcf,'Name','FiberTrimCheck','NumberTitle','off');
    end

    % これまでの二次元行列を，円形部分のみ抽出し一次元化します
    imageVectors(1,i,:) = roughImage1(k);
    imageVectors(2,i,:) = roughImage2(k);
    imageVectors(3,i,:) = roughImage3(k);
    imageVectors(4,i,:) = roughImage4(k);
end

% 時系列ベクトルをファイバごとに分けて整理します
imageVectors1 = squeeze(imageVectors(1,:,:));
imageVectors2 = squeeze(imageVectors(2,:,:));
imageVectors3 = squeeze(imageVectors(3,:,:));
imageVectors4 = squeeze(imageVectors(4,:,:));

% 目的の時間のみを取り出します
disp(num2str(TimeFrameNumber));
imageVector1 = imageVectors1(TimeFrameNumber,:);
imageVector2 = imageVectors2(TimeFrameNumber,:);
imageVector3 = imageVectors3(TimeFrameNumber,:);
imageVector4 = imageVectors4(TimeFrameNumber,:);
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