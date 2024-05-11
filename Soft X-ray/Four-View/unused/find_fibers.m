function [centers,radii]=find_fibers(IM,radiusRange)

figure;imagesc(IM);
[centers,radii] = imfindcircles(IM,radiusRange,'Sensitivity',0.95);
viscircles(centers,radii);
% viscircles(centers(1:16,:),radii(1:16));
hold on
plot(centers(:,1),centers(:,2),'*');
% disp(centers);
% disp(radii);

% 中心の近すぎる円があればより確度の高いものを選択したい
D = pdist(centers); % ユークリッド距離を計算
Z = squareform(D); % 計算結果をもとに距離行列を生成
% 最近接点との距離が閾値を下回る点を抽出
Z(Z==0) = NaN;
[M,I] = min(Z,[],2,'omitnan');
Idx = find(M<distanceThreshold);
% 最近接点が対角要素より右にあれば最近接点を消去、左にあればその点を消去
eliminateList = zeros(1,numel(I));
for i = 1:numel(I)
    if ~ismember(i,Idx)
        continue
    end
    if i < I(i) % 対角要素より右
        % 最近接点の番号(I(i))を消去リストに入れる
        if ~ismember(I(i),eliminateList)
            eliminateList(i) = I(i);
        end
    else % 対角要素より左
        % iを消去リストに入れる
        if ~ismember(i,eliminateList)
            eliminateList(i) = i;
        end
    end
end
% この動作を順番にやるかどうか
% 消去する点（行）のリストを作ってそこにぶち込む？
% 最後にその行を消去
centers(eliminateList,:) = [];
radii(eliminateList) = [];

figure;imagesc(IM);
viscircles(centers,radii);
% viscircles(centers(1:16,:),radii(1:16));
hold on
plot(centers(:,1),centers(:,2),'*');
% 問題は32個の円全てを正しく認識できていない場合（現在はこっち）


circleData = [centers,radii];
[~,I] = sort(circleData(:,1),'descend');
A = circleData(I,:);
first = A(1:4,:);
second = A(5:8,:);
third = A(9:12,:);
fourth = A(13:16,:);
[~,I1] = sort(first(:,2));
first = first(I1,:);
[~,I2] = sort(second(:,2),'descend');
second = second(I2,:);
[~,I3] = sort(third(:,2));
third = third(I3,:);
[~,I4] = sort(fourth(:,2),'descend');
fourth = fourth(I4,:);

CircleDataNew = [first;second;third;fourth];
% hold on
% plot(CircleData_New(:,1),CircleData_New(:,2),'r');

centers = fliplr(CircleDataNew(:,1:2));
radii = CircleDataNew(:,3);

end