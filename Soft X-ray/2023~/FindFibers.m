function [centers,radii]=FindFibers(IM,radiusRange)

% figure;imagesc(IM);
[centers,radii] = imfindcircles(IM,radiusRange,'Sensitivity',0.98);
% viscircles(centers,radii);
% hold on 
% plot(centers(:,1),centers(:,2),'*');
% disp(centers);
% disp(radii);

CircleData = [centers,radii];
[~,I] = sort(CircleData(:,1),'descend');
A = CircleData(I,:);
First = A(1:4,:);
Second = A(5:8,:);
Third = A(9:12,:);
Fourth = A(13:16,:);
[~,I1] = sort(First(:,2));
First = First(I1,:);
[~,I2] = sort(Second(:,2),'descend');
Second = Second(I2,:);
[~,I3] = sort(Third(:,2));
Third = Third(I3,:);
[~,I4] = sort(Fourth(:,2),'descend');
Fourth = Fourth(I4,:);

CircleData_New = [First;Second;Third;Fourth];
% hold on
% plot(CircleData_New(:,1),CircleData_New(:,2),'r');

centers = fliplr(CircleData_New(:,1:2));
radii = CircleData_New(:,3);

end