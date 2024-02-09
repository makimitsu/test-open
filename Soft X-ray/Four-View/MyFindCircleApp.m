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

up=54;

% フォルダを指定
sxrFilePath = getenv('SXR_IMAGE_DIR');
folderPath = fullfile(sxrFilePath,num2str(date));

% Excelシートにその日付のシートが存在するかどうかで条件分岐
% Excelファイルのパスとシート名を指定
excelFilePath = 'fiberPositions.xlsx';
sheetName = num2str(date);

% Excelファイルが存在するか確認
if exist(excelFilePath, 'file') == 2
    % Excelファイルを読み込み
    sheetNames = sheetnames(excelFilePath);
    
    % 指定のシートが存在するか確認
    if any(strcmp(sheetNames, sheetName))
        % シートを読み込み
        calibrationTable = readtable(excelFilePath, 'Sheet', sheetName);
        disp(['Successfully loaded the content of sheet "', sheetName, '".']);
        centerX = reshape(calibrationTable.center_X,4,8).';
        centerY = reshape(calibrationTable.center_Y,4,8).';
        CenterPositions = cat(3,centerX,centerY);
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
    CalibrationImage = imread(FiberCalibrationImagePath);% imagesc(CalibrationImage);title(CalibrationImage,'RawImage');
    CalibrationImage = wiener2(CalibrationImage,[10,10]);% imagesc(CalibrationImage);title(CalibrationImage,'WienerFiltered<RawImage');
    CalibrationImage = imadjust(CalibrationImage);% imagesc(CalibrationImage);title(CalibrationImage,'ContrastAdjusted<WienerFiltered<RawImage');
    % 円を検出します．range,sensitivity,method がパラメータです．
    % この画像なら，range=[70,~], sensitivity=0.995がよさそう．[85,~]でも検知するけど，ちょっと大きめの円を取ってしまう
    [centers,radii] = imfindcircles(CalibrationImage,[65,75],'Sensitivity',0.995,"Method","twostage");
    % 検出した円をnumber個だけ描画します．このとき，numberは32個の円が全て含まれる程度に大きく設定します．
    % numberに許される最大値はnumel(radii)です．優先度の低い円を考慮から外すためにnumberをある程度小さくする必要があります．
    % こちらを整理前のcentersとします．
    number = 60;
    
    % 検出した円を整理できないか試します．
    % 中心が近く同じ円を指すと判断されるものをgroupnumberでまとめます
    centers = [centers(1:number,:) radii(1:number) zeros([number,1])];
    groupnumber = 1;
    for i=1:number
        x1=centers(i,1);y1=centers(i,2);
        if centers(i,4) == 0
            for j=i:number
                x2=centers(j,1);y2=centers(j,2);
                r = sqrt((x1-x2)^2+(y1-y2)^2);
                if r < 80
                    centers(j,4) = groupnumber;
                end
            end
            groupnumber = groupnumber + 1;
        end
    end
    % ここでは，同じgroupnumberをもつ円情報を統合し，ひとつにします．
    % ただし同じ円を4回以上数えている場合は優先度の高い上位3つ以外は考慮から外します．
    for i = 1:groupnumber-1
        [row0,~,~] = find(centers == i);
        if numel(row0) >1
            row = row0(1);
        else 
            row = row0;
        end
        center_mean_x = sum(centers(row,1))/numel(row);
        center_mean_y = sum(centers(row,2))/numel(row);
        radius_mean = sum(centers(row,3))/numel(row);
        centers(row0,:) = [];centers(end+1,:) = [center_mean_x,center_mean_y,radius_mean,i];
    end
    radii = centers(:,3);centers = centers(:,1:2);
    % 整理後のcentersです．ただし時系列には整列されていません．
    figure;hold on;imagesc(CalibrationImage);viscircles(centers,radii,'LineWidth',1,'EnhanceVisibility',0);plot(centers(:,1),centers(:,2),'*','Color','red');set(gcf,'Name','ファイバーチェック画像に対する円検出','NumberTitle','off');hold off;
    % centersを時系列に整理します．
    % まず横軸x座標を基準に並び替えます
    centers = [centers radii];
    [~,I] = sort(centers(:,1),'descend');
    centers = centers(I,:);
    % 次に横軸が近しい4つのデータの中で，縦軸を基準に並び替えます．
    % これで，右上から左下に五十音表と同じ順に整理されます．
    for i = 1:8
        range = 1+4*(i-1):1:4+4*(i-1);
        [~,I] = sort(centers(range,2),'descend');
        I = I + 4*(i-1);
        centers(range,:) = centers(I,:);
    end
    % 次にX線フィルタ種類と時系列順に整理します．
    % X線フィルタは，右上を1，右下を2, 左下を4とします
    % CentersTimeRaps(3,:,:)はX線フィルタ3の時系列[8 2]行列です．
    % CentersTimeRaps(3,4,:)はX線フィルタ3の時系列4番目のxy座標です．
    % 蛇順のため原始的なコードが入りますが，次のブロックで上手くいっていることが確かめられます．
    CentersTimeRaps = zeros(4,8,3);
    for i = 1:4
        if i==1 || i==2
        j = i-1;
        CentersTimeRaps(i,:,:) = centers([1+j,3+j,11+j, ...
            9+j,17+j,19+j,27+j,25+j],:);
        else
        j=i-3;
        CentersTimeRaps(i,:,:) = centers([5+j,7+j,15+j, ...
            13+j,21+j,23+j,31+j,29+j],:);    
        end
    end
    CenterPositions = CentersTimeRaps(:,:,1:2);radii = repmat(min(CentersTimeRaps(:,:,3),[],'all'),[32,1]);
end

% 次に補正をかけてゆく
ShotImage = imread(ShotImagePath);

OFFSET = zeros(size(CenterPositions));
OFFSET_fiber = zeros([4,2]);
% OFFSET_fiber(1,:)=[-10,50];% [x(横);y(縦)]% ここ調整ポイント MyFindCircleとか，get_sxr_imageでdocheckを入れて手動調整してください ファイバごとで補正値違いそうなのでゆくゆく突き詰めたい
% OFFSET_fiber(2,:)=[-15,50];
% OFFSET_fiber(3,:)=[-30,45];
% OFFSET_fiber(4,:)=[-20,50];
for i=1:4
    for j=1:2
        OFFSET(i,:,j)=OFFSET_fiber(i,j);
    end
end
FixedCenters = CenterPositions+OFFSET;
RowPositions = reshape(FixedCenters(:,:,1),[32,1]);
ColumnPositions = reshape(FixedCenters(:,:,2),[32,1]);
FixedCentersForViscircles = [RowPositions ColumnPositions];

figure('Position',[200,250,800,600]);hold on;set(gcf,'Name','実際の放電画像に対する円検出（円を見つつ補正，点は補正前の中心）','NumberTitle','off');
% スライダーバー作成
handle_slider_x = uicontrol('Style','slider','Position',[250 30 300 40],'Min',-50,'Max',50,'Value',0);
handle_slider_y = uicontrol('Style','slider','Position',[10 200 40 300],'Min',-10,'Max',90,'Value',0);
% ラジオボタン作成
fib1 = uicontrol('Style','radiobutton','String','Fiber1','Position',[100 60 60 20]);
fib2 = uicontrol('Style','radiobutton','String','Fiber2','Position',[100 40 60 20]);
fib3 = uicontrol('Style','radiobutton','String','Fiber3','Position',[40 60 60 20]);
fib4 = uicontrol('Style','radiobutton','String','Fiber4','Position',[40 40 60 20]);
% 校正データセーブ用のボタン作成
% ボタンを作成
btn = uicontrol('Style', 'pushbutton', 'String', 'save data', ...
    'Position', [650, 30, 80, 30], 'Callback', {@save_calibration_data, OFFSET,CenterPositions,radii,date});

% 一旦描画する
imagesc(ShotImage,[50,up]);axis image;
% 円書く
viscircles(FixedCentersForViscircles,radii,'LineWidth',0.5,'EnhanceVisibility',0);
plot(reshape(CenterPositions(:,:,1),[32,1]),reshape(CenterPositions(:,:,2),[32,1]),'x','Color','red');
% イベントリスナーの作成
addlistener(handle_slider_x,'Value','PostSet',@(event,obj) UpdateImage(event,obj,OFFSET,CenterPositions,handle_slider_x,handle_slider_y,ShotImage,radii,fib1,fib2,fib3,fib4,up));
addlistener(handle_slider_y,'Value','PostSet',@(event,obj) UpdateImage(event,obj,OFFSET,CenterPositions,handle_slider_x,handle_slider_y,ShotImage,radii,fib1,fib2,fib3,fib4,up));

hold off;
% end

function UpdateImage(~,~,OFFSET,CenterPositions,handle_slider_x,handle_slider_y,ShotImage,radii,fib1,fib2,fib3,fib4,up)
global OFFSET_fiber
% 大事な部分を再描画する
imagesc(ShotImage,[50,up]);axis image;set(gca,'YDir','normal');
% 値の取得
value_x = round(handle_slider_x.Value);
value_y = round(handle_slider_y.Value);
% そんでいじる % [x(横);y(縦)]% ここ調整ポイント MyFindCircleとか，get_sxr_imageでdocheckを入れて手動調整してください ファイバごとで補正値違いそうなのでゆくゆく突き詰めたい

if fib1.Value == 1
    OFFSET_fiber(1,:)=[value_x,value_y];
end
if fib2.Value == 1
    OFFSET_fiber(2,:)=[value_x,value_y];
end
if fib3.Value == 1
    OFFSET_fiber(3,:)=[value_x,value_y];
end
if fib4.Value == 1
    OFFSET_fiber(4,:)=[value_x,value_y];
end

for i=1:4
    for j=1:2
        OFFSET(i,:,j)=OFFSET_fiber(i,j);
    end
end
FixedCenters = CenterPositions+OFFSET;
RowPositions = reshape(FixedCenters(:,:,1),[32,1]);
ColumnPositions = reshape(FixedCenters(:,:,2),[32,1]);
FixedCentersForViscircles = [RowPositions ColumnPositions];

viscircles(FixedCentersForViscircles,radii,'LineWidth',0.5,'EnhanceVisibility',0);hold on;
plot(reshape(CenterPositions(:,:,1),[32,1]),reshape(CenterPositions(:,:,2),[32,1]),'x','Color','red');
title(strcat('1=[',num2str(OFFSET_fiber(1,:)),'],2=[',num2str(OFFSET_fiber(2,:)),'],3=[',num2str(OFFSET_fiber(3,:)),'],4=[',num2str(OFFSET_fiber(4,:)),']'));

end

function save_calibration_data(~,~,OFFSET,CenterPositions,radii,date)
global OFFSET_fiber
for i=1:4
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
saveFile = 'fiberPositions.xlsx';
% calibrationData = [FixedCenters(:,:,1),FixedCenters(:,:,2)];
calibrationData = [number.',timing.',centerX,centerY,radii];
calibrationTable = array2table(calibrationData,'VariableNames',{'number','timing','center_X','center_Y','radii'});
% writematrix(calibrationData,saveFile,'Sheet',num2str(date));
writetable(calibrationTable,saveFile,'Sheet',num2str(date));
disp('Saved succesfully!');

end