function [] = main(t,ndata,gas,nz,plot_fitting,cal_flow,save_flow,plot_flow,factor)
%時刻/データ番号/ガス種(Ar:1,H:2)/z方向データ数(数値)/フィッティングを表示(TorF)/流速を計算(TorF)/流速を保存(TorF)/流速を表示(TorF)/矢印サイズ(数値)

%パラメータを定義
run define/parameter.m
r_measured = zeros(NofCH/4,nz);%ベクトルプロットr座標1列目data1、2列目data2
z_measured = zeros(NofCH/4,nz);%ベクトルプロットz座標1列目data1、2列目data2

%------------------------ファイル名と対応する計測点を入力--------------------------
run data_221113.m
z_measured(:,1) = 0;%data1のz座標
r_measured(:,1) = [10 10 15 15 20 20 25 25];%data1のr座標
% filename1 = '/Users/Ryo/Doppler/211111/ascii/shot8_466us_w=2_gain=3800.asc';%data1ファイル名
% filename1 = '/Volumes/experiment/results/Doppler/Andor/IDSP/221113/shot2_474us_w=2_gain=3000.asc';%data1ファイル名
% filename1 = '/Volumes/experiment/results/Doppler/Andor/IDSP/221101/shot35_474_w=2_3200.asc';%data1ファイル名
% z_measured(:,1) = -2.1;%data1のz座標
% r_measured(:,1) = [12.5, 15, 17.5, 20, 22.5, 25, 27.5];%data1のr座標
% filename2 = '/Users/Ryo/Doppler/211110/ascii/shot50_464us_w=2_gain=3800.asc';%data2ファイル名
% filename2 = '/Volumes/experiment/results/Doppler/Andor/IDSP/221101/shot31_474_w=2_3100.asc';%data2ファイル名
% z_measured(:,2) = 2.1;%data2のz座標
% r_measured(:,1) = [22.5, 25, 27.5, 30, 32.5, 35, 37.5];%data2のr座標
%----------------------------------------------------------------------------

%較正データからプローブ各チャンネル対応箇所を取得
if gas == 1
    %アルゴンの時
    center_file = '221114_Xe_96120_calibration.txt';%中心データファイル名
elseif gas == 2
    %水素の時
    center_file = 'Hbeta_calibration.txt';%中心データファイル名
end
center = importdata(center_file);%中心座標を取得
% centerX = repmat(center(:,3),1,nz);%チャンネル対応中心相対X座標
if nz == 1
    centerY = center(:,2);%チャンネル対応中心Y座標
elseif nz == 2
    centerY = repmat(center(:,2),1,nz);%チャンネル対応中心Y座標
end

%data1からスペクトル群を取得
data1 = importdata(filename1);
X1 = data1(:,1);%X1(ピクセル)軸を定義
[LofX1,n1]=size(data1);%X1軸の長さを取得
L1 = zeros(LofX1,NofCH);%L1(波長)軸を定義
px2nm = zeros(NofCH,1);%nm/pixel
if gas == 1
    %アルゴンの時
    lambda = [lambda1 lambda2];
    for i = 1:NofCH
        pixel = [center(i,3) center(i,4)];
        p = polyfit(pixel,lambda,1);
        px2nm(i,1) = p(1);
        L1(:,i) = polyval(p,X1);
    end
elseif gas == 2
    %水素の時
    for i = 1:NofCH
        px2nm(i,1) = 0.00536;
        L1(:,i) = px2nm(i,1)*(X1-center(i,3))+lambda0;
    end
end
spectrum1 = zeros(LofX1,NofCH);%data1の分光結果を入れる
for i = 1:NofCH
    spectrum1(:,i) = ...
        sum(data1(:,round(centerY(i,1)-Width):round(centerY(i,1)+Width)),2);
end

%data2からスペクトル群を取得
if nz == 2
    data2 = importdata(filename2);
    X2 = data2(:,1);%X2(ピクセル)軸を定義
    [LofX2,n2]=size(data2);%X2軸の長さを取得
    spectrum2=zeros(LofX2,NofCH);%data2の分光結果を入れる
    for i = 1:NofCH
        spectrum2(:,i) = ...
            sum(data2(:,round(centerY(i,2)-Width):round(centerY(i,2)+Width)),2);
    end
end

%波長方向の移動平均から離れた外れ値を移動平均に変更(ノイズ除去)
M1 = movmean(spectrum1,l);
% M1 = (M1*l - spectrum1) / (l-1);
if nz == 2
    M2 = movmean(spectrum2,l);
    % M2 = (M2*l - spectrum2) / (l-1);
end
for i = 1:NofCH
    for j = 1:LofX1
        if spectrum1(j,i) > M1(j,i)*1.1
            spectrum1(j,i) = M1(j,i);
        elseif spectrum1(j,i) < M1(j,i)*0.9
            spectrum1(j,i) = M1(j,i);
        end
        if nz == 2
            if spectrum2(j,i) > M2(j,i)*1.1
                spectrum2(j,i) = M2(j,i);
            elseif spectrum2(j,i) < M2(j,i)*0.9
                spectrum2(j,i) = M2(j,i);
            end
        end
    end
end


%バックグラウンド引いてないとき用のノイズ除去
for j = 1:NofCH
    goukei = 0;
    for i = 1:200
        goukei = goukei + spectrum1(i,j);
    end
    heikin = goukei/200;
    for i = 1:1024
        spectrum1(i,j) = spectrum1(i,j) - heikin;
    end
end
% spectrum1(1024,8) 


%計算結果取得用の配列を用意
amp = zeros(NofCH,nz);%振幅1列目data1、2列目data2
shift = zeros(NofCH,nz);%中心1列目data1、2列目data2
sigma = zeros(NofCH,nz);%シグマ1列目data1、2列目data2
V = zeros(NofCH/4,nz*2);%1列目data1・Vz(km/s)、2列目data1・Vr(km/s)、3列目data2・Vz(km/s)、4列目data2・Vr(km/s)
absV = zeros(NofCH/4,nz);%1列目data1・V(km/s)、2列目data2・V(km/s)
T = zeros(NofCH,nz);%1列目data1・温度(eV)、2列目data1・温度(eV)
Ti = zeros(NofCH/4,nz);%1列目data1・温度(eV)、2列目data1・平均温度(eV)
offset = zeros(NofCH/4,nz);

%data1のガウスフィッティング
if plot_fitting
    figure('Position',[300 50 1000 1000])
end
for k = 1:NofCH/4
    %0度ペアスペクトルからオフセットを検出
    for i = 1:2
        S1 = [L1(:,(k-1)*4+i) spectrum1(:,(k-1)*4+i)]; %(k-1)*4+i番目のチャンネルの[波長,強度]
        s1 = size(S1);%S1の[行数,列数]
        MAX1 = max(spectrum1(:,(k-1)*4+i)); %スペクトルの最大値
        j = 1;
        while j < s1(1)+1 %SNの悪いデータを除く
            if S1(j,2) < MAX1*0.7
                S1(j,:) = [];
            else
                j = j+1;
            end
            s1 = size(S1);
        end
        f = fit(S1(:,1),S1(:,2),'gauss1');
        coef=coeffvalues(f);
        amp((k-1)*4+i,1) = coef(1);
        shift((k-1)*4+i,1) = coef(2)-lambda0;
        sigma((k-1)*4+i,1) = sqrt(coef(3)^2-(center((k-1)*4+i,5)*px2nm((k-1)*4+i,1))^2);
        T((k-1)*4+i,1) = 1.69e8*A*(2*sigma((k-1)*4+i,1)*sqrt(2*log(2))/lambda0)^2;
    end
    offset(k,1) = (shift((k-1)*4+1,1) + shift((k-1)*4+2,1))/2%対向視線から得られたオフセット[nm]
    %オフセットを引いてガウスフィッティング
    for i = 1:4
        S1 = [L1(:,(k-1)*4+i)-offset(k,1) spectrum1(:,(k-1)*4+i)]; %(k-1)*4+i番目のチャンネルの[波長,強度]
        s1 = size(S1);%S1の[行数,列数]
        MAX1 = max(spectrum1(:,(k-1)*4+i)); %スペクトルの最大値
        j = 1;
        while j < s1(1)+1 %SNの悪いデータを除く
            if S1(j,2) < MAX1*0.5
                S1(j,:) = [];
            else
                j = j+1;
            end
            s1 = size(S1);
        end
        f = fit(S1(:,1),S1(:,2),'gauss1');
        coef=coeffvalues(f);
        amp((k-1)*4+i,1) = coef(1);
        shift((k-1)*4+i,1) = coef(2)-lambda0;
        sigma((k-1)*4+i,1) = sqrt(coef(3)^2-(center((k-1)*4+i,5)*px2nm((k-1)*4+i,1))^2);
        T((k-1)*4+i,1) = 1.69e8*A*(2*sigma((k-1)*4+i,1)*sqrt(2*log(2))/lambda0)^2;
        if plot_fitting
            subplot(NofCH/4,4,(k-1)*4+i);
            plot(f,S1(:,1),S1(:,2));
            xline(lambda0);
            title([num2str(T((k-1)*4+i,1)),' eV'])
            legend('off')
            xlabel('Wavelength [nm]')
            ylabel('Intensity [cnt]')
        end
    end
end
if plot_fitting
    sgtitle('Fitting data1 (Horizontal：Channel Vertical：Position)')
end

%data1の流速、温度を計算
if cal_flow
    for i = 1:NofCH/4
        set = (i-1)*4;
        Va = -shift(set+2,1)/lambda0*Vc;
        Vb = -shift(set+1,1)/lambda0*Vc;
        V(i,1) = Va;%Vz
        V(i,2) = 0;%Vr
        absV(i,1) = sqrt(V(i,1)^2 + V(i,2)^2);
        %        fprintf('(z, r) = (%.2f, %.2f) cm\n (Vz, Vr) = (%.2f, %.2f) km/s\n'...
        %             ,z(i,1),r(i,1),V(i,1),V(i,2));
        %             fprintf('Ti: %.2f, %.2f, %.2f, %.2f eV\n'...
        %                     ,T(set+1,1),T(set+2,1),T(set+3,1),T(set+4,1));
        %         Ti(i,1) = (T(set+1,1)+T(set+2,1)+T(set+3,1)+T(set+4,1))/4;
        Ti(i,1) = trimmean(T(set+1:set+4,1),50);
%         fprintf('%.2f\n',Ti(i,1));
    end
end

%data1の計算結果を保存
if save_flow
    save(['magflow/mat/save_',num2str(t),'us_',num2str(ndata),'.mat'],'V','Ti')
end

%data2のガウスフィッティング
if nz == 2
    if plot_fitting
        figure('Position',[300 50 1000 1000])
    end
    for k = 1:NofCH/4
        %0度ペアスペクトルからオフセットを検出
        for i = 1:2
            S2 = [L1(:,(k-1)*4+i) spectrum2(:,(k-1)*4+i)]; %(k-1)*4+i番目のチャンネルの[波長,強度]
            s2 = size(S1);%S1の[行数,列数]
            MAX2 = max(spectrum2(:,(k-1)*4+i)); %スペクトルの最大値
            j = 1;
            while j < s2(1)+1 %SNの悪いデータを除く
                if S2(j,2) < MAX2*0.5
                    S2(j,:) = [];
                else
                    j = j+1;
                end
                s2 = size(S2);
            end
            f = fit(S2(:,1),S2(:,2),'gauss1');
            coef=coeffvalues(f);
            amp((k-1)*4+i,2) = coef(1);
            shift((k-1)*4+i,2) = coef(2)-lambda0;
            sigma((k-1)*4+i,2) = sqrt(coef(3)^2-(center((k-1)*4+i,5)*px2nm((k-1)*4+i,1))^2);
            T((k-1)*4+i,2) = 1.69e8*A*(2*sigma((k-1)*4+i,2)*sqrt(2*log(2))/lambda0)^2;
        end
        offset(k,2) = (shift((k-1)*4+1,2) + shift((k-1)*4+2,2))/2;%対向視線から得られたオフセット[nm]
        %オフセットを引いてガウスフィッティング
        for i = 1:4
            S2 = [L1(:,(k-1)*4+i) spectrum2(:,(k-1)*4+i)]; %(k-1)*4+i番目のチャンネルの[波長,強度]
            s2 = size(S1);%S1の[行数,列数]
            MAX2 = max(spectrum2(:,(k-1)*4+i)); %スペクトルの最大値
            j = 1;
            while j < s2(1)+1 %SNの悪いデータを除く
                if S2(j,2) < MAX2*0.5
                    S2(j,:) = [];
                else
                    j = j+1;
                end
                s2 = size(S2);
            end
            f = fit(S2(:,1),S2(:,2),'gauss1');
            coef=coeffvalues(f);
            amp((k-1)*4+i,2) = coef(1);
            shift((k-1)*4+i,2) = coef(2)-lambda0;
            sigma((k-1)*4+i,2) = sqrt(coef(3)^2-(center((k-1)*4+i,5)*px2nm((k-1)*4+i,1))^2);
            T((k-1)*4+i,2) = 1.69e8*A*(2*sigma((k-1)*4+i,2)*sqrt(2*log(2))/lambda0)^2;
            if plot_fitting
                subplot(NofCH/4,4,(k-1)*4+i);
                plot(f,S2(:,1),S2(:,2));
                xline(lambda0);
                title([num2str(T((k-1)*4+i,2)),' eV'])
                legend('off')
                xlabel('Wavelength [nm]')
                ylabel('Intensity [cnt]')
            end
        end
    end
    if plot_fitting
        sgtitle('Fitting data1 (Horizontal：Channel Vertical：Position)')
    end
    
    %data2の流速、温度を計算
    if cal_flow
        for i = 1:NofCH/4
            set = (i-1)*4;
            Va = -shift(set+4,2)/lambda0*Vc;
            Vb = -shift(set+3,2)/lambda0*Vc;
            V(i,3) = (Va-Vb)/(2*cos(Angle*pi/180));%Vz
            V(i,4) = -(Va+Vb)/(2*sin(Angle*pi/180));%Vr
            absV(i,2) = sqrt(V(i,3)^2 + V(i,4)^2);
            %        fprintf('(z, r) = (%.2f, %.2f) cm\n (Vz, Vr) = (%.2f, %.2f) km/s\n'...
            %             ,z(i,2),r(i,2),V(i,3),V(i,4));
            %             fprintf('Ti: %.2f, %.2f, %.2f, %.2f eV\n'...
            %                     ,T(set+1,2),T(set+2,2),T(set+3,2),T(set+4,2));
            Ti1=(T(set+1,2)+T(set+2,2)+T(set+3,2)+T(set+4,2))/4;
            fprintf('%.2f\n',Ti1);
        end
    end
end

%流速、温度をプロット
if cal_flow
    if plot_flow
        figure('Position',[600 150 300 600])
        Ticon = repmat(Ti(:,1),1,5);%等高線図用平均イオン温度
        zcon = z_measured(1,1);
        s = pcolor([zcon-2 zcon-1 zcon zcon+1 zcon+2],r_measured,Ticon);
        s.FaceColor = 'interp';
        s.EdgeAlpha = 0;
        colormap('jet')
        colorbar
        hold on
        plot(z_measured,r_measured,'xr');
        hold on
        if nz == 1
            q = quiver(z_measured,r_measured,V(:,1)*factor,V(:,2)*factor);
        end
        if nz == 2
            q = quiver(z_measured,r_measured,[V(:,1),V(:,3)]*factor,[V(:,2),V(:,4)]*factor);
        end
        q.LineWidth = 2;
        q.MaxHeadSize = 10;
        q.AutoScale = 'off';
        q.Color = 'k';
        xlim([min(z_measured,[],'all')-2.5 max(z_measured,[],'all')+2.5])
        ylim([min(r_measured,[],'all')-2.5 max(r_measured,[],'all')+2.5])
        title('Ion Flow [km/s]','Color','black','FontWeight','bold')
        absV = round(absV,1);
        for j = 1:nz
            for i = 1:NofCH/4
                txt = text(z_measured(i,j)+0.5,r_measured(i,j)-0.2,num2str(absV(i,j)));
                txt.FontSize = 18;
                txt.Color = 'k';
                txt.FontWeight = 'bold';
            end
        end
        xlabel('Z [cm]')
        ylabel('R [cm]')
        ax = gca;
        ax.FontSize = 14;
        grid on
        daspect([1 1 1])
        hold off
    end
end
