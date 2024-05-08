function [] = main_multi_at_t(t,NZ,factor,inversion_method,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare)
%時間,Z方向データ数,矢印サイズ,逆変換(1~5),スペクトル描画,結果描画,等高線(1)/三次元(2),分析描画,比較描画

NR = 7;%計測点数
ndata = 1;%データ番号
z_measured = zeros(NZ,NR);
r_measured = zeros(NZ,NR);
Vmean_z = zeros(NZ,NR);
Vmean_r = zeros(NZ,NR);
Tz = zeros(NZ,NR);
Tr = zeros(NZ,NR);
absV = zeros(NZ,NR);

for i = 1:NZ
    for j = 1:NR
        [Vmean_z(i,j),Vmean_r(i,j),Tz(i,j),Tr(i,j)] = main(t,i,j,ndata,inversion_method,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare);
        if i == 1
            z_measured(i,j) = -4.25/2/100;
        end
        if i == 2
            z_measured(i,j) = 4.25/2/100;
        end
        r_measured(i,j) = (10 + 2.5*j)/100;
        absV(i,j) = sqrt(Vmean_z(i,j)^2 + Vmean_r(i,j)^2);
    end
end
absV = round(absV,1);
Tz = round(Tz,1);
Tr = round(Tr,1);

%-------------------磁気面、フローの2次元分布をプロット-------------------
run magprobe/main_magflow.m %磁気面

% figure('Position',[600 150 230 600])
plot(z_measured,r_measured,'xr');
hold on

for i = 1:NZ
    q = quiver(z_measured(i,:),r_measured(i,:),Vmean_z(i,:)*factor,Vmean_r(i,:)*factor);
    q.LineWidth = 2;
    q.MaxHeadSize = 10;
    q.AutoScale = 'off';
    if i == 1
        q.Color = 'r';
    end
    if i == 2
        q.Color = 'b';
    end
    for j = 1:NR
        if i == 1
            txt1 = text(z_measured(i,j)+0.02,r_measured(i,j)-0.002,[num2str(absV(i,j))]);
            txt1.FontSize = 40;
            txt1.Color = 'r';
            txt1.FontWeight = 'bold';
        end
        if i == 2
            txt2 = text(z_measured(i,j)+0.005,r_measured(i,j)-0.002,[num2str(absV(i,j))]);
            txt2.FontSize = 40;
            txt2.Color = 'b';
            txt2.FontWeight = 'bold';
        end
    end
    hold on
end
% title('Flow Vector Map [km/s]')
% xlabel('z [cm]')
% ylabel('r [cm]')
% xlim([-5 5])
% ylim([5 35])
% ax = gca;
% ax.FontSize = 12;
hold off

% %-------------------擬似温度の一次元分布をプロット-------------------
% plot_temp_at_t(t,NR,Tz,Tr)
end
