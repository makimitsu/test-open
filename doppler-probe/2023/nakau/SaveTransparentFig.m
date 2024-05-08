function SaveTransparentFig(fig)
% Save figure with transparent background
set(fig, 'Color', 'white');
image = getframe(fig);
image = image.cdata;
[h,w,~] = size(image);
Alpha = ones(h,w);
white = [255 255 255]; % RGB value

for i=1:h
    for j=1:w
        if nnz(squeeze(image(i,j,:)) == white') == 3
            Alpha(i,j) = 0;
        end
    end
end

% savefile
filter = {'*.png';'*.*'};
[file, ~] = uiputfile(filter);
imwrite(image,string(file),'Alpha',Alpha);
end