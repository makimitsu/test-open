function check_position(date,sxrFilename)

% ファイバー位置の校正用コード
% エクセルからファイバー位置を読み込んで画像を切り取り、切り取った画像をsubplot等で表示
% calibrationImage.tiffあたりを読み込んで切り取る機能を入れたい
% ついでに読み込んだ情報を元に中心位置、半径あたりも記録
% calibrationImage.tiffから中心位置、半径を読み取り、記録するコードは別？
% こっちはExcelファイルからデータを読み込んで表示するためだけのコードにする

projectionNumber = 80;


% 生画像の取得
rawImage = imread(sxrFilename);


% 校正用画像からファイバーの位置（＋半径）を取得
positionPath = '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/Four-View/fiberPositions.xlsx';
centers = readmatrix(positionPath,'Sheet',num2str(date),'Range','C2:D33');
Center = zeros(4,8,2);
for i = 1:4
    Center(i,:,:) = centers(1+8*(i-1):8+8*(i-1),:);
end
IW = 60;

Center = round(Center);

% 切り取った画像を格納するための配列
timeSeries = zeros(4,8,2*IW,2*IW);

% バックグラウンドノイズのデータを取得
backgroundImage = cast(rawImage(1:2*IW,1:2*IW,1),'double');
backgroundNoise = ones(size(backgroundImage))*mean(backgroundImage,'all');

% 切り取った画像のうち実際に使う部分（ファイバー部分）を切り出し
k = find_circle(projectionNumber/2);
imageVectors = zeros(4,8,numel(k));

f1=figure;
% f1.Position = [900,200,500,500];
f1.Position = [200,250,1060,500];
% f2=figure;
% f2.Position = [900,200,500,500];

% 再構成計算に使用可能なサイズに変換するための解像度を取得
resolution = projectionNumber/(IW*2);

% 画像データの行列化
for i=1:8
    % 画像切り取りの縦方向・横方向範囲を指定
    horizontalRange1 = Center(1,i,1)-IW+1:Center(1,i,1)+IW;
    verticalRange1 = Center(1,i,2)-IW+1:Center(1,i,2)+IW;
    horizontalRange2 = Center(2,i,1)-IW+1:Center(2,i,1)+IW;
    verticalRange2 = Center(2,i,2)-IW+1:Center(2,i,2)+IW;
    horizontalRange3 = Center(3,i,1)-IW+1:Center(3,i,1)+IW;
    verticalRange3 = Center(3,i,2)-IW+1:Center(3,i,2)+IW;
    horizontalRange4 = Center(4,i,1)-IW+1:Center(4,i,1)+IW;
    verticalRange4 = Center(4,i,2)-IW+1:Center(4,i,2)+IW;
    % 生画像を切り取ってファイバー数×フレーム数×IW×IWの配列に格納
    timeSeries(1,i,:,:) = rawImage(verticalRange1,horizontalRange1,1);
    timeSeries(2,i,:,:) = rawImage(verticalRange2,horizontalRange2,1);
    timeSeries(3,i,:,:) = rawImage(verticalRange3,horizontalRange3,1);
    timeSeries(4,i,:,:) = rawImage(verticalRange4,horizontalRange4,1);

    % バックグラウンドノイズ成分を差し引く
    singleImage1 = squeeze(timeSeries(1,i,:,:))-backgroundNoise;
    singleImage2 = squeeze(timeSeries(2,i,:,:))-backgroundNoise;
    singleImage3 = squeeze(timeSeries(3,i,:,:))-backgroundNoise;
    singleImage4 = squeeze(timeSeries(4,i,:,:))-backgroundNoise;
    grayImage1 = flipud(singleImage1(:,:));
    grayImage2 = flipud(singleImage2(:,:));
    grayImage3 = flipud(singleImage3(:,:));
    grayImage4 = flipud(singleImage4(:,:));
    grayImage1(grayImage1<0) = 0;
    grayImage2(grayImage2<0) = 0;
    grayImage3(grayImage3<0) = 0;
    grayImage4(grayImage4<0) = 0;
    % 再構成計算に使用可能なサイズに変換
    roughImage1 = imresize(grayImage1, resolution, 'nearest');
    roughImage1 = cast(roughImage1,'double');
    roughImage2 = imresize(grayImage2, resolution, 'nearest');
    roughImage2 = cast(roughImage2,'double');
    roughImage3 = imresize(grayImage3, resolution, 'nearest');
    roughImage3 = cast(roughImage3,'double');
    roughImage4 = imresize(grayImage4, resolution, 'nearest');
    roughImage4 = cast(roughImage4,'double');

% 切り出した画像を表示
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

    % ベクトル化
    imageVectors(1,i,:) = roughImage1(k);
    imageVectors(2,i,:) = roughImage2(k);
    imageVectors(3,i,:) = roughImage3(k);
    imageVectors(4,i,:) = roughImage4(k);
end

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