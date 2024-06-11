% function MyFindCircleApp
clearvars -except date
% close all
global OFFSET_fiber

%%%%実験オペレーションの取得
prompt = {'Date:'};
definput = {'','',''};
if exist('date','var')
    definput{1} = num2str(date);
end
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    return
end
date = str2double(cell2mat(answer));

up=70;

% フォルダを指定
sxrFilePath = getenv('SXR_IMAGE_DIR');
folderPath = fullfile(sxrFilePath,num2str(date));

% Excelシートにその日付のシートが存在するかどうかで条件分岐
% Excelファイルのパスとシート名を指定
excelFilePath = 'fiberPositions_2view.xlsx';
sheetName = num2str(date);

% 中心位置のインデックス
centerIdx1 = [1,3,7,5,9,11,15,13];
centerIdx2 = [2,4,8,6,10,12,16,14];

centerIdx = [centerIdx1;centerIdx2];

% Excelファイルが存在するか確認
if exist(excelFilePath, 'file') == 2
    % Excelファイルを読み込み
    sheetNames = sheetnames(excelFilePath);
    
    % 指定のシートが存在するか確認
    if any(strcmp(sheetNames, sheetName))
        % シートを読み込み
        calibrationTable = readtable(excelFilePath, 'Sheet', sheetName);
        disp(['Successfully loaded the content of sheet "', sheetName, '".']);
        centerX = reshape(calibrationTable.center_X,8,2).';
        centerY = reshape(calibrationTable.center_Y,8,2).';
        centerPositions = cat(3,centerX,centerY);
        radii = calibrationTable.radii;
        disp('Select one image file to test the position calibration.');
        % uigetfileダイアログを開いて.tifファイルを選択
        [fileName, folderPath] = uigetfile(fullfile(folderPath, '*.tif'));
        
        % 選択がキャンセルされた場合の処理
        if isequal(fileName, 0)
            disp('ファイルの選択がキャンセルされました。');
            return
        else
            ShotImagePath = fullfile(folderPath, fileName);
            disp(['テスト用画像: ', ShotImagePath]);
        end
        doCalculation = false;
    else
        disp(['Sheet "', sheetName, '" does not exist.']);
        doCalculation = true;
    end
else
    disp(['Excel file "', excelFilePath, '" does not exist.']);
    doCalculation = true;
end

if doCalculation
    disp('Select two image files to test the position calibration.')
    % uigetfileダイアログを開いて.tifファイルを選択
    [fileNames, folderPath] = uigetfile(fullfile(folderPath, '*.tif'), 'MultiSelect', 'on');
    
    % 選択がキャンセルされた場合の処理
    if isequal(fileNames, 0)
        disp('ファイルの選択がキャンセルされました。');
        return
    else
        % 2つのファイルに対する処理
        if numel(fileNames) == 2
            % 1つ目のファイルに対する処理
            FiberCalibrationImagePath = fullfile(folderPath, fileNames{1});
            disp(['校正画像: ', FiberCalibrationImagePath]);
            
            % 2つ目のファイルに対する処理
            ShotImagePath = fullfile(folderPath, fileNames{2});
            disp(['テスト用画像: ', ShotImagePath]);
        else
            disp('2つのファイルを選択してください。');
        end
    end
    
    % このブロックでは，ファイバーの較正画像から32個の円を漏れなく検出することを目指します
    % 個別のTIF画像で調整すべき変数は，画像path, imfindcircles(RadiusRenge, Sensitivity), numberです．
    % 画像を読み込み，前処理を施します．
    calibrationImage = imread(FiberCalibrationImagePath);% imagesc(CalibrationImage);title(CalibrationImage,'RawImage');
        % 円を検出します．range,sensitivity,method がパラメータです．
    % この画像なら，range=[70,~], sensitivity=0.995がよさそう．[85,~]でも検知するけど，ちょっと大きめの円を取ってしまう
    if date == 230316
        radiusRange = [85,105];
    else    
        radiusRange = [65,75];
    end
    [centerPositions,radii] = find_fibers_16ch(calibrationImage,radiusRange);
   
end

% 次に補正をかけてゆく
ShotImage = imread(ShotImagePath);

OFFSET = zeros(size(centerPositions));
OFFSET_fiber = zeros([2,2]);

for i=1:2
    for j=1:2
        OFFSET(i,:,j)=OFFSET_fiber(i,j);
    end
end
FixedCenters = centerPositions+OFFSET;
RowPositions = reshape(FixedCenters(:,:,1),[16,1]);
ColumnPositions = reshape(FixedCenters(:,:,2),[16,1]);
FixedCentersForViscircles = [RowPositions ColumnPositions];

figure('Position',[200,250,800,600]);hold on;set(gcf,'Name','実際の放電画像に対する円検出（円を見つつ補正，点は補正前の中心）','NumberTitle','off');
% スライダーバー作成
handle_slider_x = uicontrol('Style','slider','Position',[250 30 300 40],'Min',-50,'Max',50,'Value',0);
% handle_slider_y = uicontrol('Style','slider','Position',[10 200 40 300],'Min',-10,'Max',90,'Value',0);
handle_slider_y = uicontrol('Style','slider','Position',[10 200 40 300],'Min',-90,'Max',10,'Value',0);
% ラジオボタン作成
fib1 = uicontrol('Style','radiobutton','String','Fiber1','Position',[100 60 60 20]);
fib2 = uicontrol('Style','radiobutton','String','Fiber2','Position',[100 40 60 20]);
% fib3 = uicontrol('Style','radiobutton','String','Fiber3','Position',[40 60 60 20]);
% fib4 = uicontrol('Style','radiobutton','String','Fiber4','Position',[40 40 60 20]);
% 校正データセーブ用のボタン作成
% ボタンを作成
btn = uicontrol('Style', 'pushbutton', 'String', 'save data', ...
    'Position', [650, 30, 80, 30], 'Callback', {@save_calibration_data, OFFSET,centerPositions,radii,date});
btn2 = uicontrol('Style', 'pushbutton', 'String', 'check', ...
    'Position', [550, 30, 80, 30], 'Callback', {@check_cropped_image, OFFSET,centerPositions,radii,date,ShotImage,centerIdx});

% 一旦描画する
imagesc(ShotImage,[50,up]);axis image;set(gca,'YDir','reverse');
% 円書く
viscircles(FixedCentersForViscircles,radii,'LineWidth',0.5,'EnhanceVisibility',0);
plot(reshape(centerPositions(:,:,1),[16,1]),reshape(centerPositions(:,:,2),[16,1]),'x','Color','red');
% イベントリスナーの作成
addlistener(handle_slider_x,'Value','PostSet',@(event,obj) UpdateImage(event,obj,OFFSET,centerPositions,handle_slider_x,handle_slider_y,ShotImage,radii,fib1,fib2,up));
addlistener(handle_slider_y,'Value','PostSet',@(event,obj) UpdateImage(event,obj,OFFSET,centerPositions,handle_slider_x,handle_slider_y,ShotImage,radii,fib1,fib2,up));

hold off;
% end

function UpdateImage(~,~,OFFSET,CenterPositions,handle_slider_x,handle_slider_y,ShotImage,radii,fib1,fib2,up)
global OFFSET_fiber
% 大事な部分を再描画する
imagesc(ShotImage,[50,up]);axis image;%set(gca,'YDir','normal');
% 値の取得
value_x = round(handle_slider_x.Value);
% value_y = round(handle_slider_y.Value);
value_y = -1*round(handle_slider_y.Value);
% そんでいじる % [x(横);y(縦)]% ここ調整ポイント MyFindCircleとか，get_sxr_imageでdocheckを入れて手動調整してください ファイバごとで補正値違いそうなのでゆくゆく突き詰めたい

if fib1.Value == 1
    OFFSET_fiber(1,:)=[value_x,value_y];
end
if fib2.Value == 1
    OFFSET_fiber(2,:)=[value_x,value_y];
end

for i=1:2
    for j=1:2
        OFFSET(i,:,j)=OFFSET_fiber(i,j);
    end
end
FixedCenters = CenterPositions+OFFSET;
RowPositions = reshape(FixedCenters(:,:,1),[16,1]);
ColumnPositions = reshape(FixedCenters(:,:,2),[16,1]);
FixedCentersForViscircles = [RowPositions ColumnPositions];

viscircles(FixedCentersForViscircles,radii,'LineWidth',0.5,'EnhanceVisibility',0);hold on;
plot(reshape(CenterPositions(:,:,1),[16,1]),reshape(CenterPositions(:,:,2),[16,1]),'x','Color','red');
title(strcat('1=[',num2str(OFFSET_fiber(1,:)),'],2=[',num2str(OFFSET_fiber(2,:)),']',']'));

end

function save_calibration_data(~,~,OFFSET,CenterPositions,radii,date)
global OFFSET_fiber
for i=1:2
    for j=1:2
        OFFSET(i,:,j)=OFFSET_fiber(i,j);
    end
end

FixedCenters = CenterPositions+OFFSET;
centerX = reshape(FixedCenters(:,:,1).',[],1);
centerY = reshape(FixedCenters(:,:,2).',[],1);
idx = 1:numel(centerX);
number = ceil(idx./8);
timing = mod(idx,8);
timing(timing==0) = 8;
saveFile = 'fiberPositions_2view.xlsx';
% calibrationData = [FixedCenters(:,:,1),FixedCenters(:,:,2)];
calibrationData = [number.',timing.',centerX,centerY,radii];
calibrationTable = array2table(calibrationData,'VariableNames',{'number','timing','center_X','center_Y','radii'});
% writematrix(calibrationData,saveFile,'Sheet',num2str(date));
writetable(calibrationTable,saveFile,'Sheet',num2str(date));
disp('Saved succesfully!');

end

function check_cropped_image(~,~,OFFSET,CenterPositions,radii,date,ShotImage,centerIdx)
global OFFSET_fiber
for i=1:2
    for j=1:2
        OFFSET(i,:,j)=OFFSET_fiber(i,j);
    end
end
rawImage = ShotImage;
FixedCenters = CenterPositions+OFFSET;
crop_image(date,FixedCenters,radii,rawImage,centerIdx)
end

function [] = crop_image(date,FixedCenters,radii,rawImage,centerIdx)

projectionNumber = 50;

% % 位置情報ファイルからファイバーの位置（＋半径）を取得
% positionPath = 'fiberPositions.xlsx';
% positionData = readmatrix(positionPath,'Sheet',num2str(date),'Range','C2:E33');
% Center = zeros(4,8,2);
% for i = 1:4
%     Center(i,:,:) = positionData(1+8*(i-1):8+8*(i-1),1:2);
% end
Center = FixedCenters;
IW = radii(1);
% IW = positionData(1,3);

Center = round(Center);

% 切り取った画像を格納するための配列
timeSeries = zeros(2,8,2*IW,2*IW);


% バックグラウンドノイズのデータを取得
backgroundImage = cast(rawImage(1:2*IW,1:2*IW,1),'double');
backgroundNoise = ones(size(backgroundImage))*mean(backgroundImage,'all');

% 切り取った画像のうち実際に使う部分（ファイバー部分）を切り出し
k = find_circle(projectionNumber/2);

    f1=figure;
    f1.Position = [200,250,1060,500];


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

positionIdx = reshape(1:16,8,2).';
centerIdxMatrix = rot90(reshape(1:16,2,8).',3);

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
        figure(f1);
        i_str = num2str(i);
        title1 = strcat('1,',i_str);
        title2 = strcat('2,',i_str);

        % subplot(4,8,4*(i-1)+1);imagesc(roughImage1);title(title1);yticks([]);xticks([]);
        % subplot(4,8,4*(i-1)+2);imagesc(roughImage2);title(title2);yticks([]);xticks([]);
        % subplot(4,8,4*(i-1)+3);imagesc(roughImage3);title(title3);yticks([]);xticks([]);
        % subplot(4,8,4*(i-1)+4);imagesc(roughImage4);title(title4);yticks([]);xticks([]);
        subplot(2,8,positionIdx(centerIdxMatrix==centerIdx(1,i)));imagesc(roughImage1);title(title1);yticks([]);xticks([]);clim([0 20]);
        subplot(2,8,positionIdx(centerIdxMatrix==centerIdx(2,i)));imagesc(roughImage2);title(title2);yticks([]);xticks([]);clim([0 20]);


    % % ベクトル化
    % imageVectors(1,i,:) = roughImage1(k);
    % imageVectors(2,i,:) = roughImage2(k);
    % imageVectors(3,i,:) = roughImage3(k);
    % imageVectors(4,i,:) = roughImage4(k);
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