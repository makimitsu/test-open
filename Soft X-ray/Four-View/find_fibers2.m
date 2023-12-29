function [CenterPositions,radius]=find_fibers2(imageFile)

radiusRange = [65, 75];
NumberOfFrame = 8;

% このブロックでは，ファイバーの較正画像から32個の円を漏れなく検出することを目指します.
% 実行環境によって調整すべき変数は，較正画像pathであるFiberCalibrationImagePathです．
% すべての円を漏れなく検出できないとき調整すべき変数はimfindcircles(RadiusRenge, Sensitivity), numberです．
% 画像を読み込み，前処理を施します．
% FiberCalibrationImagePath = "/Volumes/experiment/results/X-ray/230704/Calibration/shot006.tif";
% CalibrationImage = imread(FiberCalibrationImagePath);
% imageFile = imread(imageFilePath);
imageFile = wiener2(imageFile,[10,10]);
imageFile = imadjust(imageFile);
% 円を検出します．range,sensitivity,method がパラメータです．
% この画像なら，range=[70,~], sensitivity=0.995がよさそう．[85,~]でも検知するけど，ちょっと大きめの円を取ってしまう
[centers,radii] = imfindcircles(imageFile,radiusRange,'Sensitivity',0.995,"Method","twostage"); % 半径は必ず整数
% 検出した円をnumber個だけ描画します．このとき，numberは32個の円が全て含まれる程度に大きく設定します．
% numberに許される最大値はnumel(radii)で，TIF画像によりますがおよそ200程度ぽいです．優先度の低い円を考慮から外すためにnumberをある程度小さくする必要があります．
% こちらを整理前のcentersとします．numberを調整したので検出に漏れはありませんがダブりがあります．

% figure;hold on;
number = 60;
% imagesc(imageFile);viscircles(centers(1:number,:),radii(1:number,:));
% plot(centers(1:number,1),centers(1:number,2),'*','Color','red');
% hold off;

% ダブりで検出した円を整理します．
% 中心が近く同じ円を指すと判断されるものをgroupnumberでまとめます
% ここで，その判断基準であるcenter同士の距離の閾値(if r<80)は上手くいくようにTIF画像ごとに変える必要があります．
% 1,2行目が中心座標、3行目が半径、4行目がgroupnumbaerであるような配列circleInformationを作成
circleInformation = [centers(1:number,:) radii(1:number) zeros([number,1])];
groupnumber = 1;
for i=1:number
    x1=circleInformation(i,1);y1=circleInformation(i,2);
    if circleInformation(i,4) == 0
        for j=i:number
            x2=circleInformation(j,1);y2=circleInformation(j,2);
            r = sqrt((x1-x2)^2+(y1-y2)^2);
            if r < 80
                circleInformation(j,4) = groupnumber;
            end
        end
        groupnumber = groupnumber + 1;
    end
end
% 同じgroupnumberをもつ円情報を統合し，ひとつにします．
% 同じ円を2回以上数えている場合は確度の最も高いものを採用します．
for i = 1:groupnumber-1
    [row0,~,~] = find(circleInformation == i);
    rows = find(circleInformation(:,4)==i);
    if numel(rows) >1
        row = rows(1);
    else 
        row = rows;
    end
    % if numel(row0) >1
    %     row = row0(1);
    % else 
    %     row = row0;
    % end
    center_mean_x = sum(circleInformation(row,1))/numel(row);
    center_mean_y = sum(circleInformation(row,2))/numel(row);
    radius_mean = sum(circleInformation(row,3))/numel(row);
    circleInformation(row0,:) = [];circleInformation(end+1,:) = [center_mean_x,center_mean_y,radius_mean,i];
end
radii = circleInformation(:,3);circleInformation = circleInformation(:,1:2);

% % 整理後のcentersです．ただし時系列には整列されていません．
% figure;hold on;
% imagesc(imageFile);
% viscircles(circleInformation,radii);
% plot(circleInformation(:,1),circleInformation(:,2),'*','Color','red');
% hold off;

% centersを時系列に整理します．
% まず横軸x座標を基準に並び替えます
circleInformation = [circleInformation radii];
[~,I] = sort(circleInformation(:,1),'descend');
circleInformation = circleInformation(I,:);
% 次に横軸x座標が近しい4つのデータの中で，縦軸を基準に並び替えます．
% これで，右上から左下に五十音表と同じ順に整理されます．
for i = 1:NumberOfFrame
    range = 1+4*(i-1):1:4+4*(i-1);
    [~,I] = sort(circleInformation(range,2),'descend');
    I = I + 4*(i-1);
    circleInformation(range,:) = circleInformation(I,:);
end
% 次にX線フィルタ種類と時系列順に整理します．
% X線フィルタは，右上を1，右下を2, 左上を3とします
% CentersTimeRaps(3,:,:)はX線フィルタ3の時系列[4 2]行列です．
% CentersTimeRaps(3,4,:)はX線フィルタ3の時系列4番目のxy座標です．
% 蛇順のため原始的なコードが入りますが，次のブロックで上手くいっていることが確かめられます．
CenterPositions = zeros(4,NumberOfFrame,3);
if NumberOfFrame == 4
    CenterPositions(1,:,:) = circleInformation([4,2,10,12],:);% ここの順番reservation
    CenterPositions(2,:,:) = circleInformation([3,1,9,11],:);
    CenterPositions(3,:,:) = circleInformation([8,6,14,16],:);
    CenterPositions(4,:,:) = circleInformation([7,5,13,15],:);
else
    CenterPositions(1,:,:) = circleInformation([4,2,10,12,20,18,26,28],:);
    CenterPositions(2,:,:) = circleInformation([3,1,9,11,19,17,25,27],:);
    CenterPositions(3,:,:) = circleInformation([8,6,14,16,24,22,30,32],:);
    CenterPositions(4,:,:) = circleInformation([7,5,13,15,23,21,29,31],:);
end
% disp(num2str(CenterPositions(:,:,3)));
radius = round(mean(CenterPositions(:,:,3),'all'));
% radius = min(CenterPositions(:,:,3),[],'all');
CenterPositions(:,:,3) = [];% これはradius分をすてたということです。

%::::中心位置の補正を入れます::::
OFFSET = zeros(size(CenterPositions));
OFFSET_fiber = zeros([4,2]);
% 230929
% OFFSET_fiber(1,:)=[-11,51];% [x(横);y(縦)]% ここ調整ポイント MyFindCircleとか，get_sxr_imageでdocheckを入れて手動調整してください ファイバごとで補正値違いそうなのでゆくゆく突き詰めたい
% OFFSET_fiber(2,:)=[-12,51];
% OFFSET_fiber(3,:)=[-36,44];
% OFFSET_fiber(4,:)=[-32,44];

% 230920/shot27
% OFFSET_fiber(1,:)=[-7,46];% [x(横);y(縦)]% ここ調整ポイント MyFindCircleとか，get_sxr_imageでdocheckを入れて手動調整してください ファイバごとで補正値違いそうなのでゆくゆく突き詰めたい
% OFFSET_fiber(2,:)=[-7,40];
% OFFSET_fiber(3,:)=[-27,44];
% OFFSET_fiber(4,:)=[-12,42];

% 230920/shot15
% OFFSET_fiber(1,:)=[-7,39];% [x(横);y(縦)]% ここ調整ポイント MyFindCircleとか，get_sxr_imageでdocheckを入れて手動調整してください ファイバごとで補正値違いそうなのでゆくゆく突き詰めたい
% OFFSET_fiber(2,:)=[9,45];
% OFFSET_fiber(3,:)=[-17,39];
% OFFSET_fiber(4,:)=[-17,39];

% 231215/
OFFSET_fiber(1,:)=[-12,33];
OFFSET_fiber(2,:)=[-12,33];
OFFSET_fiber(3,:)=[-21,33];
OFFSET_fiber(4,:)=[-21,33];

% 231216/
OFFSET_fiber(1,:)=[-7,47];
OFFSET_fiber(2,:)=[-2,47];
OFFSET_fiber(3,:)=[-21,38];
OFFSET_fiber(4,:)=[-20,38];

for i=1:4
    for j=1:2
        OFFSET(i,:,j)=OFFSET_fiber(i,j);
    end
end
CenterPositions = CenterPositions+OFFSET;
%::::::::::::::::::::::::::

end