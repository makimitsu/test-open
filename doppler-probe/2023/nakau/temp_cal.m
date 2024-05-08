NofCH = 32;%チャンネル数
NT = 3;%計測トリガ時刻数
NData = 2;%同一計測トリガでの計測ショット数
nz = 1;

T = zeros(NT,1);
err = zeros(NT,1);
T_multi = zeros(NofCH,NData);
T_mean = zeros(NT,1);

for t = 1:NT
    for ndata = 1:NData
        load(['magflow/mat/temp_',num2str(462+4*t),'us_',num2str(ndata),'.mat'],'temp')
        T_multi(:,ndata) = temp(1:NofCH,1);
    end
    T = mean(T_multi,2);
    T_mean(t) = mean(T)%平均温度
    err(t) = sqrt(var(T_mean))/sqrt(NT*NData)%標準誤差
end

errorbar(T_mean,err)

