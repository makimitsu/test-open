function [] = main_multi_velocity(inversion_method,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare)
%逆変換(1~5),スペクトル描画,結果描画,等高線(1)/三次元(2),分析描画,比較描画

mp = 1.67e-27;%陽子質量(kg)
kB = 1.60e-19;%ボルツマン定数(J/eV)
A = 40;%質量数

nz = 1;%計測点Z方向番号(nz = 1:-2.1cm, nz = 2:+2.1cm)
NT = 6;%計測トリガ時間数(1~6)
NR = 7;%計測点数(1~7)
NData = 5;%データ数(1~5)
Vmean_z_t = zeros(NT,NR);%時刻tの平均Vmean_z
Vmean_r_t = zeros(NT,NR);%時刻tの平均Vmean_r
Sz_t = zeros(NT,NR);%時刻tのVmean_zの標準偏差
Sr_t = zeros(NT,NR);%時刻tのVmean_rの標準偏差
Vmean_z2 = zeros(NData,NR);
Vmean_r2 = zeros(NData,NR);
Vmean_z2_t = zeros(NT,NR);%時刻tの平均V^2mean_z
Vmean_r2_t = zeros(NT,NR);%時刻tの平均V^2mean_r
Sz2_t = zeros(NT,NR);%時刻tのV^2mean_zの標準偏差
Sr2_t = zeros(NT,NR);%時刻tのV^2mean_rの標準偏差
% Vmean_t = zeros(NT,NR);%時刻tの平均Vmean

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

for t = 1:NT
    %読み込む
    load(['mat/save_fwhm_minus_',num2str(463+2*t),'us.mat'],'Vmean_z','Vmean_r')
    for j = 1:NR
        for k = 1:NData
            Vmean_z2(k,j) = Vmean_z(k,j)^2*10^6*A*mp/(2*kB);
            Vmean_r2(k,j) = Vmean_r(k,j)^2*10^6*A*mp/(2*kB);
        end
    end
    Vmean_z_t(t,:) = mean(Vmean_z);
    Vmean_r_t(t,:) = mean(Vmean_r);
    Vmean_z2_t(t,:) = mean(Vmean_z2);
    Vmean_r2_t(t,:) = mean(Vmean_r2);
    Sz_t(t,:) = std(Vmean_z);
    Sr_t(t,:) = std(Vmean_r);
    Sz2_t(t,:) = std(Vmean_z2);
    Sr2_t(t,:) = std(Vmean_r2);
end
%-------------------擬似温度の一次元分布をプロット-------------------
plot_velocity(NT,NR,Vmean_z2_t,Vmean_r2_t,Sz2_t,Sr2_t)

end
