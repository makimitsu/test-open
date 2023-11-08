function [] = main_multi_temp(inversion_method,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare)
%逆変換(1~5),スペクトル描画,結果描画,等高線(1)/三次元(2),分析描画,比較描画

nz = 1;%計測点Z方向番号(nz = 1:-2.1cm, nz = 2:+2.1cm)
NT = 6;%計測トリガ時間数(1~6)
NR = 7;%計測点数(1~7)
NData = 5;%データ数(1~5)
Tz_t = zeros(NT,NR);%時刻tの平均Tz
Tr_t = zeros(NT,NR);%時刻tの平均Tr
Sz_t = zeros(NT,NR);%時刻tのTzの標準偏差
Sr_t = zeros(NT,NR);%時刻tのTrの標準偏差

% Vmean_z = zeros(NData,NR);
% Vmean_r = zeros(NData,NR);
% Tz = zeros(NData,NR);
% Tr = zeros(NData,NR);

% for t = 1:NT
%     for j = 1:NR
%         for k = 1:NData
%             [Vmean_z(k,j),Vmean_r(k,j),Tz(k,j),Tr(k,j)] = main(463+2*t,nz,j,k,inversion_method,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare);
%         end
%     end
%     %保存
%     save(['mat/save_fwhm_plus_',num2str(463+2*t),'us.mat'],'Vmean_z','Vmean_r','Tz','Tr')
%     Tz_t(t,:) = mean(Tz);
%     Tr_t(t,:) = mean(Tr);
%     Sz_t(t,:) = std(Tz);
%     Sr_t(t,:) = std(Tr);
% end

% for t = 1:NT
%     %読み込む
%     load(['mat/save_fwhm_minus_',num2str(463+2*t),'us.mat'],'Tz','Tr')
%     Tz_t(t,:) = mean(Tz);
%     Tr_t(t,:) = mean(Tr);
%     Sz_t(t,:) = std(Tz);
%     Sr_t(t,:) = std(Tr);
% end
% %-------------------擬似温度の一次元分布をプロット-------------------
% plot_temp(NT,NR,Tz_t,Tr_t,Sz_t,Sr_t)

for t = 1:NT
    %読み込む
    load(['mat/save_fwhm_minus_',num2str(463+2*t),'us.mat'],'Vmean_z','Vmean_r')
    Tz_t(t,:) = mean(Vmean_z);
    Tr_t(t,:) = mean(Vmean_r);
    Sz_t(t,:) = std(Vmean_z);
    Sr_t(t,:) = std(Vmean_r);
end
%-------------------速度成分の一次元分布をプロット-------------------
plot_velociy(NT,NR,Tz_t,Tr_t,Sz_t,Sr_t)


end
