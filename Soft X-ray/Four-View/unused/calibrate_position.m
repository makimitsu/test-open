function calibrate_position(date)

fiberPositionFile = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/',num2str(date),'/PositionCheck.tif');
calibrationImage = imread(fiberPositionFile);

[Center,IW] = find_fibers2(calibrationImage,[65,75]);

centers = zeros(32,2);
for i = 1:4
    centers(1+8*(i-1):8+8*(i-1),:) = Center(i,:,:);
end

positionData = [centers IW*ones(32,1)];
positionPath = '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/Four-View/fiberPositions.xlsx';
writematrix(positionData,positionPath,'Sheet',num2str(date),'Range','C2:E33');
% disp(positionData);

end