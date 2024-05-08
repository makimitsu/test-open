function [Vmean_z,Vmean_r,Ti_z,Ti_r] = main(t,nz,nr,ndata,inversion_method,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare)
%時間,計測点Z方向番号,計測点R方向番号,データ番号,逆変換(1~5),スペクトル描画,結果描画,等高線(1)/三次元(2),分析描画,比較描画
tic
time = t;
r_measured = (100 + 25*nr);
%パラメータを定義
run define/parameter.m

%スペクトルデータPを取り込む
run data/dataP.m

%観測スペクトルを描画
if draw_spectra
    run draw/spectra.m
end

%W(F->P)を作成
run make/funcW.m

%---------------Invertion theory----------------

%pinvを使ってW^(-1)を求める。
if inversion_method == 1
    run method/pseudo.m
    run method/goto0.m
end

%Tikhonov 0th
if inversion_method == 2
    run method/Tikhonov0.m
    run method/goto0.m
end

%Tikhonov 1st
if inversion_method == 3
    run method/Tikhonov1.m
    run method/goto0.m
end

%Tikhonov 2nd
if inversion_method == 4
    run method/Tikhonov2.m
    run method/goto0.m
end

%非負制約SIRT
if inversion_method == 5
    run method/SIRT.m
end

%minimum Fischer information
if inversion_method == 6
    run method/MFI.m
    run method/goto0.m
end

%結果の表示
if inversion_method ~= false
    drReF = reshape(ReF, [numVx,numVy]);%描画用に整形
    drReF = flipud(drReF);
    [Summity,Summitx] = find(drReF == max(drReF(:)));
    
    %再構成結果を描画
    if draw_result
        run draw/result.m
    end
    
    %再構成結果から得られるスペクトルを計算、スペクトルデータと比較
    if draw_compare
        run draw/compare.m
    end
    %流速を計算
    Fsum = 0;
    Vmean_z = 0;
    Vmean_r = 0;
    for i = 1:numVy
        for j = 1:numVx
            if drReF(i,j)> drReF(Summity,Summitx)*0.5
                Vmean_z = Vmean_z + drReF(i,j)*Vx(j);
                Vmean_r = Vmean_r + drReF(i,j)*Vy(i);
                Fsum = Fsum + drReF(i,j);
            end
        end
    end
    Vmean_z = Vmean_z/Fsum;
    Vmean_r = Vmean_r/Fsum;
    Ti_z = 0;
    Ti_r = 0;
    %     fprintf('(Vz = %.1f, Vr = %.1f)[km/s].\n',Vmean_z,Vmean_r);
    %     %擬似温度を計算(分散から計算)
    %     Vdisp = zeros(2,1);
    %     for i = 1:numVy
    %         for j = 1:numVx
    %             if drReF(i,j)> drReF(Summity,Summitx)*0
    %                 Vdisp(1,1) = Vdisp(1,1) + drReF(i,j)*(Vx(j)-Vmean_z)^2;
    %                 Vdisp(2,1) = Vdisp(2,1) + drReF(i,j)*(Vy(i)-Vmean_r)^2;
    %             end
    %         end
    %     end
    %
    %     Vdisp(1,1) = Vdisp(1,1)/Fsum - (mean(inst(1+4*(nr-1):4+4*(nr-1)),'all')*deltaLambda/lambda0*Vc)^2;
    %     Vdisp(2,1) = Vdisp(2,1)/Fsum - (mean(inst(1+4*(nr-1):4+4*(nr-1)),'all')*deltaLambda/lambda0*Vc)^2;
    % %     Vdisp(1,1) = Vdisp(1,1)/Fsum;%装置関数なし
    % %     Vdisp(2,1) = Vdisp(2,1)/Fsum;%装置関数なし
    %     Ti_z = Vdisp(1,1)*10^6*A*mp/(2*kB);
    %     Ti_r = Vdisp(2,1)*10^6*A*mp/(2*kB);
    %     fprintf('(Tz = %.1f, Tr = %.1f)[eV].\n',Ti_z,Ti_r);
    
%     %擬似温度を計算(半値幅から計算)
%     W_0 = zeros(numLambda,numVx*numVy);%0度スペクトルを計算
%     for i = 1:numLambda
%         Unit_0 = [cos(0) sin(0)];
%         for j = 1:numVx*numVy
%             %W(:,j)の対応速度番号x,yを取得：対応速度はVx(x),Vy(y)
%             x = 1 + idivide(j-1, int8(numVy), 'floor');
%             y = int8(mod(j,numVy));
%             if y == 0
%                 y = numLambda;
%             end
%             V = [Vx(x), Vy(y)];
%             D = dot(V, Unit_0);%視線方向速度
%             %      W(i,j) = 1 - abs(Lambda(r)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc));%一次関数フィルタ
%             W_0(i,j) = 1 - ((ShiftLambda(i,1)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc)))^2;%二次関数フィルタ
%             if W_0(i,j) < 0
%                 W_0(i,j) = 0;
%             end
%         end
%     end
%     P_0 = W_0*ReF;
%     SummitP_0 = find(P_0 == max(P_0));
%     %プロット
% %     figure('Position',[300 550 600 400])
% %     plot(Vx,P_0/max(P_0), 'b','LineWidth',2)
% %     hold on
% %     yline(0.5, '--r','LineWidth',2);
% %     xlabel('V_{Z} [km/s]')
% %     ylabel('Relative intensity')
% %     xlim([-60 60])
% %     xticks(-60:20:60)
% %     yticks(0:0.1:1)
% %     ax = gca;
% %     ax.FontSize = 18;
% %     hold off
%     for i = 1:numLambda - SummitP_0
%         if P_0(i+SummitP_0) < max(P_0)*0.5
%             I_xmax = i+SummitP_0;
%             break
%         end
%     end
%     fwhm_xmax = Vx(I_xmax);%半値x右端
%     for i = 1:SummitP_0
%         if P_0(SummitP_0-i+1) < max(P_0)*0.5
%             I_xmin = SummitP_0-i+1;
%             break
%         end
%     end
%     fwhm_xmin = Vx(I_xmin);%半値x左端
%     fwhm_x = fwhm_xmax - fwhm_xmin;
%     sigma_x = fwhm_x/2.35;
% %     S_x = sigma_x^2 - (mean(inst(1+4*(nr-1):4+4*(nr-1)),'all')*deltaLambda/lambda0*Vc)^2;
%     S_x = sigma_x^2;%装置関数なし
%     Ti_z = S_x*10^6*A*mp/(2*kB);
%     
%     W_90 = zeros(numLambda,numVx*numVy);%0度スペクトルを計算
%     for i = 1:numLambda
%         Unit_90 = [cos(pi/2) sin(pi/2)];
%         for j = 1:numVx*numVy
%             %W(:,j)の対応速度番号x,yを取得：対応速度はVx(x),Vy(y)
%             x = 1 + idivide(j-1, int8(numVy), 'floor');
%             y = int8(mod(j,numVy));
%             if y == 0
%                 y = numLambda;
%             end
%             V = [Vx(x), Vy(y)];
%             D = dot(V, Unit_90);%視線方向速度
%             %      W(i,j) = 1 - abs(Lambda(r)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc));%一次関数フィルタ
%             W_90(i,j) = 1 - ((ShiftLambda(i,1)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc)))^2;%二次関数フィルタ
%             if W_90(i,j) < 0
%                 W_90(i,j) = 0;
%             end
%         end
%     end
%     P_90 = W_90*ReF;
%     SummitP_90 = find(P_90 == max(P_90));
%     %プロット
% %     figure('Position',[300 550 600 400])
% %     plot(Vx,P_90/max(P_90), 'b','LineWidth',2)
% %     hold on
% %     yline(0.5, '--r','LineWidth',2);
% %     xlabel('V_{R} [km/s]')
% %     ylabel('Relative intensity')
% %     xlim([-60 60])
% %     xticks(-60:20:60)
% %     yticks(0:0.1:1)
% %     ax = gca;
% %     ax.FontSize = 18;
% %     hold off
%     for i = 1:numLambda - SummitP_90
%         if P_90(i+SummitP_90) < max(P_90)*0.5
%             I_ymax = i+SummitP_90;
%             break
%         end
%     end
%     fwhm_ymax = Vy(I_ymax);%半値y右端
%     for i = 1:SummitP_90
%         if P_0(SummitP_90-i+1) < max(P_90)*0.5
%             I_ymin = SummitP_90-i+1;
%             break
%         end
%     end
%     fwhm_ymin = Vy(I_ymin);%半値x左端
%     fwhm_y = fwhm_ymax - fwhm_ymin;
%     sigma_y = fwhm_y/2.35;
% %     S_y = sigma_y^2 - (mean(inst(1+4*(nr-1):4+4*(nr-1)),'all')*deltaLambda/lambda0*Vc)^2;
%     S_y = sigma_y^2;%装置関数なし
%     Ti_r = S_y*10^6*A*mp/(2*kB);
%         fprintf('Reconstructd: (Tz = %.1f, Tr = %.1f)[eV].\n',Ti_z,Ti_r);
end
toc
end
