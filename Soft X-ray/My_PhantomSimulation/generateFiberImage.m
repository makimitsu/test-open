function [FiberImage,Image] = generateFiberImage(gm2d,EE)

PhantomEE = reshape(EE,1,[]);

FiberImageArray = gm2d*(PhantomEE)';
FiberImage=awgn(FiberImageArray,5,'measured'); % 5 related to 20%; 10 related to 10%;
% FiberImage = MyAddWhiteGaussianNoise(FiberImageArray,3);
% FiberImage = FiberImageArray;
% FiberImage(FiberImage<0)=0;
FiberImage = FiberImage.*1.2;
k = find_circle(50/2);
Image = zeros(50);
Image(k) = FiberImage;

figure(3);imagesc(Image);set(gcf,'Name','巻き戻しファイバー画像','NumberTitle','off');

Image = imnlmfilt(Image,"ComparisonWindowSize",21,"SearchWindowSize",41);
FiberImage = Image(k);
% figure(4);imagesc(Image);set(gcf,'Name','フィルタ後巻き戻しファイバー画像','NumberTitle','off');

end

function k = find_circle(L)
R = zeros(2*L);
for i = 1:2*L
    for j = 1:2*L
        R(i,j) = sqrt((L-i+0.5)^2+(j-L-0.5)^2);
    end
end
k = find(R<L);
end