load(['230807_Xe4835_calibation.mat'],'calib_Ar');
calib_Ar(1,:) = [];
writematrix(calib_Ar,'230807_Xe4835_calibation.txt','Delimiter','tab')
