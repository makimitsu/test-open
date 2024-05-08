calib_Ar = ["CH_number" "CH_position" "lambda1_position" "nm/pixel" "instrument" "strength"];

%--【Input】----
date = 230807;%校正実験日

for i = 1:96
    if  exist([num2str(date),'/calibation',num2str(i),'.mat'],'file') ~= 0
        load([num2str(date),'/calibation',num2str(i),'.mat'],'cal_result')
        calib_Ar = cat(1,calib_Ar,cal_result);
    end
end
save([num2str(date),'_Xe4835_calibation','.mat'],'calib_Ar')
