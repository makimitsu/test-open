function [V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_fit,save_flow)
%実験日/ICCD変数/計測点配列/pathname/offsetを表示/ガウスフィッティングをプロット/流速を保存

%------フィッティング調整用【input】-------
w_CH = 4;%【input】チャンネル方向(Y方向)足し合わせ幅
l_mov = 7;%【input】波長方向移動平均長さ
th_ratio = 0.8;%【input】フィッティング時閾値
l_L1 = 41;%【input】波長軸の切り取り長さ
d_L1 = 25;%【input】波長方向位置ずれ調整
d_CH = 2;%【input】CH方向位置ずれ調整

%物理定数
Vc = 299792.458;%光速(km/s)
Angle = 30;%入射角(度)
switch ICCD.line
    case 'Ar'%アルゴンの時
        A = 40;%原子量
        lambda0 = 480.602;%使用スペクトル(nm)
        lambda1 = 480.7019;%校正ランプスペクトル(nm)
        lambda2 = 479.2619;%校正ランプスペクトル(nm)
    case 'H'%水素の時
        A = 1;%原子量
        lambda0 = 486.135;%使用スペクトル(nm)
        warning('Sorry, not ready for H experiment.')%ICCD.lineの入力エラー
        return;
    otherwise
        warning('Input error in ICCD.line.')%ICCD.lineの入力エラー
        return;
end
center = load_calibration(date,ICCD);

dir1 = [pathname.IDSP '/' num2str(date)];%ディレクトリ1
if mpoints.n_z == 1
    filename1 = [dir1 '/shot' num2str(ICCD.shot) '_' num2str(ICCD.trg) 'us_w=' num2str(ICCD.exp_w) '_gain=' num2str(ICCD.gain) '.asc'];%ICCDファイル名
    if not(exist(filename1,"file"))
        warning([filename1,' does not exist.']);
        V_i = char.empty;
        absV = char.empty;
        T_i = char.empty;
        return
    end
    %data1からスペクトル群を取得
    data1 = importdata(filename1);
end

%----------------------------------------------------------------------------
% centerX = repmat(center(:,3),1,mpoints.n_z);%チャンネル対応中心相対X座標
switch mpoints.n_z
    case 1
        centerY = center(:,2);%チャンネル対応中心Y座標
    case 2
        warning('Sorry, not ready for n_z = 2.')%ICCD.lineの入力エラー
        return
        centerY = repmat(center(:,2),1,mpoints.n_z);%チャンネル対応中心Y座標
end

X1 = data1(:,1);%X1(ピクセル)軸を定義
[l_X1,~]=size(data1);%X1軸の長さを取得
L1 = zeros(l_X1,mpoints.n_CH);%L1(波長)軸を定義
L1_shaped = zeros(l_L1,mpoints.n_CH);%切り取った波長軸
px2nm = zeros(mpoints.n_CH,1);%nm/pixel
switch ICCD.line
    case 'Ar'%アルゴンの時
        %     lambda = [lambda1 lambda2];
        for i = 1:mpoints.n_CH
            %         pixel = [center(i,3) center(i,4)];
            %         p = polyfit(pixel,lambda,1);
            %         px2nm(i,1) = p(1);
            %         L1(:,i) = polyval(p,X1);
            px2nm(i,1) = center(i,4);
            L1(:,i) = lambda1 - px2nm(i,1)*(X1(:,1)-center(i,3));
            L1_shaped(:,i) = L1(round(center(i,3))-(l_L1-1)/2+d_L1:round(center(i,3))+(l_L1-1)/2+d_L1,i);
        end
    case 'H'%水素の時
        for i = 1:mpoints.n_CH
            px2nm(i,1) = 0.00536;
            L1(:,i) = px2nm(i,1)*(X1-center(i,3))+lambda0;
        end
end
spectrum1 = zeros(l_L1,mpoints.n_CH);%data1の分光結果を入れる
for i = 1:mpoints.n_CH
    spectrum1(:,i) = ...
        sum(data1(round(center(i,3))-(l_L1-1)/2+d_L1:round(center(i,3))+(l_L1-1)/2+d_L1,round(centerY(i,1)+d_CH-w_CH):round(centerY(i,1)+d_CH+w_CH)),2);
end

% %data2からスペクトル群を取得
% if mpoints.n_z > 1
%     data2 = importdata(filename2);
%     X2 = data2(:,1);%X2(ピクセル)軸を定義
%     [LofX2,~]=size(data2);%X2軸の長さを取得
%     spectrum2=zeros(LofX2,mpoints.n_CH);%data2の分光結果を入れる
%     for i = 1:mpoints.n_CH
%         spectrum2(:,i) = ...
%             sum(data2(:,round(centerY(i,2)-width):round(centerY(i,2)+width)),2);
%     end
% end

%波長方向の移動平均から離れた外れ値を移動平均に変更(ノイズ除去)
M1 = movmean(spectrum1,l_mov);
% M1 = (M1*l - spectrum1) / (l-1);
% if mpoints.n_z == 2
%     M2 = movmean(spectrum2,mov_l);
%     % M2 = (M2*l - spectrum2) / (l-1);
% end
for i = 1:mpoints.n_CH
    for j = 1:l_L1
        if spectrum1(j,i) > M1(j,i)*1.1
            spectrum1(j,i) = M1(j,i);
        elseif spectrum1(j,i) < M1(j,i)*0.9
            spectrum1(j,i) = M1(j,i);
        end
        % if mpoints.n_z == 2
        %     if spectrum2(j,i) > M2(j,i)*1.1
        %         spectrum2(j,i) = M2(j,i);
        %     elseif spectrum2(j,i) < M2(j,i)*0.9
        %         spectrum2(j,i) = M2(j,i);
        %     end
        % end
    end
end

%計算結果取得用の配列を用意
amp = zeros(mpoints.n_CH,mpoints.n_z);%振幅1列目data1、2列目data2
shift = zeros(mpoints.n_CH,mpoints.n_z);%中心1列目data1、2列目data2
sigma = zeros(mpoints.n_CH,mpoints.n_z);%シグマ1列目data1、2列目data2
V_i = zeros(mpoints.n_r,mpoints.n_z*2);%1列目data1・Vz(km/s)、2列目data1・Vr(km/s)、3列目data2・Vz(km/s)、4列目data2・Vr(km/s)
absV = zeros(mpoints.n_r,mpoints.n_z);%1列目data1・V(km/s)、2列目data2・V(km/s)
T_CH = zeros(mpoints.n_CH,mpoints.n_z);%1列目data1・温度(eV)、2列目data1・温度(eV)
T_i = zeros(mpoints.n_r,mpoints.n_z);%1列目data1・温度(eV)、2列目data1・平均温度(eV)
offset = zeros(mpoints.n_r,mpoints.n_z);
error_CH = zeros(mpoints.n_CH,1);

for k = 1:mpoints.n_r
    %0度ペアスペクトルからオフセットを検出
    for i = 1:2
        i_CH = (k-1)*4+i;%CH番号
        S1 = [L1_shaped(:,i_CH) spectrum1(:,i_CH)]; %[波長,強度]
        s1 = size(S1);%S1の[行数,列数]
        MAX1 = max(spectrum1(:,i_CH)); %スペクトルの最大値
        j = 1;
        while j < s1(1)+1 %SNの悪いデータを除く
            if S1(j,2) < MAX1*th_ratio
                S1(j,:) = [];
            else
                j = j+1;
            end
            s1 = size(S1);
        end
        try
            f = fit(S1(:,1),S1(:,2),'gauss1');
        catch ME
            warning(['Fitting failed in CH ',num2str(i_CH),'.']);
            error_CH(i_CH,1) = 1;
            break
        end
        coef=coeffvalues(f);
        shift(i_CH,1) = coef(2)-lambda0;
    end
    offset(k,1) = (shift((k-1)*4+1,1) + shift((k-1)*4+2,1))/2;%対向視線から得られたオフセット[nm]
end

%異常値を持つoffsetを内挿して置き換え
ori_offset = offset;%変更前を保管
buf = offset(:,1);%offsetをコピー
ind_0 = find(buf == 0);%0を検出
ind_non0 = find(buf ~= 0);%non0を検出
buf(ind_0) = [];%0を削除
buf = filloutliers(buf,"linear");%外れ値を線型内挿で置き換え
offset(ind_non0,1) = buf;%offsetを更新
offset(ind_0,1) = mean(offset(ind_non0,1));%offsetを更新
for k = 1:mpoints.n_r
    if offset(k,1) ~= ori_offset(k,1)
        disp(['Offset in Row ',num2str(k),' is replaced from ',num2str(ori_offset(k,1)),' to ',num2str(offset(k,1))'.'])
    end
end

%data1のガウスフィッティング
if plot_fit
    figure('Position',[300 50 1000 1000],'visible','on')
    sgtitle(['Fitting data (Horizontal：View Line, Vertical：Measured Position)',newline, ...
        'shot',num2str(ICCD.shot),'-',num2str(ICCD.trg),'us-w=',num2str(ICCD.exp_w),'-gain=',num2str(ICCD.gain),'.asc'])
end
for k = 1:mpoints.n_r
    % if (error_CH((k-1)*4+1,1) == 0) && (error_CH((k-1)*4+2,1) == 0)%対抗視線のフィッティングに成功
        if show_offset
            disp(['Offset = ',num2str(offset(k,1)),' in Row ',num2str(k),'.'])
        end
        %オフセットを引いてガウスフィッティング
        for i = 1:4
            i_CH = (k-1)*4+i;%CH番号
            S1 = [L1_shaped(:,i_CH)-offset(k,1) spectrum1(:,i_CH)]; %[波長,強度]
            s1 = size(S1);%S1の[行数,列数]
            MAX1 = max(spectrum1(:,i_CH)); %スペクトルの最大値
            j = 1;
            while j < s1(1)+1 %SNの悪いデータを除く
                if S1(j,2) < MAX1*0.5
                    S1(j,:) = [];
                else
                    j = j+1;
                end
                s1 = size(S1);
            end
            try
                f = fit(S1(:,1),S1(:,2),'gauss1');
                coef=coeffvalues(f);
                amp(i_CH,1) = coef(1);
                shift(i_CH,1) = coef(2)-lambda0;
                sigma(i_CH,1) = sqrt(coef(3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
                T_CH(i_CH,1) = 1.69e8*A*(2*sigma(i_CH,1)*sqrt(log(2))/lambda0)^2;
                if plot_fit
                    subplot(mpoints.n_r,4,i_CH);
                    plot(f,S1(:,1),S1(:,2));
                    % plot(S1(:,1),S1(:,2));
                    xline(lambda0);
                    title([num2str(T_CH(i_CH,1)),' eV'])
                    legend('off')
                    xlabel('Wavelength [nm]')
                    ylabel('Intensity [cnt]')
                end
            catch ME
                warning(['Fitting failed in CH ',num2str(i_CH),'.']);
                break
            end
        end
    % else
    %     warning(['Offset failed in Row ',num2str(k),'.']);
    % end
end
if plot_fit
    if save_fit
        time = round(ICCD.trg+ICCD.exp_w/2);%計測時刻
        if not(exist([pathname.fig,'/fit/',num2str(date)],'dir'))
            mkdir(sprintf("%s/fit", pathname.fig), sprintf("%s", num2str(date)));
        end
        saveas(gcf,[pathname.fig,'/fit/',num2str(date),'/','shot', num2str(ICCD.shot),'_',num2str(time),'us_fit.png'])
        hold off
        close
    else
        hold off
    end
end

%data1の流速、温度を計算
for i = 1:mpoints.n_r
    set = (i-1)*4;
    Va = -shift(set+4,1)/lambda0*Vc;
    Vb = -shift(set+3,1)/lambda0*Vc;
    V_i(i,1) = (Va-Vb)/(2*cos(Angle*pi/180));%Vz
    V_i(i,2) = -(Va+Vb)/(2*sin(Angle*pi/180));%Vr
    absV(i,1) = sqrt(V_i(i,1)^2 + V_i(i,2)^2);
    %        fprintf('(z, r) = (%.2f, %.2f) cm\n (Vz, Vr) = (%.2f, %.2f) km/s\n'...
    %             ,z(i,1),r(i,1),V(i,1),V(i,2));
    %             fprintf('T_i: %.2f, %.2f, %.2f, %.2f eV\n'...
    %                     ,T(set+1,1),T(set+2,1),T(set+3,1),T(set+4,1));
    %         T_i(i,1) = (T(set+1,1)+T(set+2,1)+T(set+3,1)+T(set+4,1))/4;
    T_i(i,1) = trimmean(T_CH(set+1:set+4,1),50);
    %         fprintf('%.2f\n',T_i(i,1));
end

% %data2のガウスフィッティング
% if mpoints.n_z == 2
%     if plot_fit
%         figure('Position',[300 50 1000 1000])
%     end
%     for k = 1:mpoints.n_r
%         %0度ペアスペクトルからオフセットを検出
%         for i = 1:2
%             i_CH = (k-1)*4+i;%CH番号
%             S2 = [L1_shaped(:,i_CH) spectrum2(:,i_CH)]; %[波長,強度]
%             s2 = size(S1);%S1の[行数,列数]
%             MAX2 = max(spectrum2(:,i_CH)); %スペクトルの最大値
%             j = 1;
%             while j < s2(1)+1 %SNの悪いデータを除く
%                 if S2(j,2) < MAX2*0.5
%                     S2(j,:) = [];
%                 else
%                     j = j+1;
%                 end
%                 s2 = size(S2);
%             end
%             f = fit(S2(:,1),S2(:,2),'gauss1');
%             coef=coeffvalues(f);
%             amp(i_CH,2) = coef(1);
%             shift(i_CH,2) = coef(2)-lambda0;
%             sigma(i_CH,2) = sqrt(coef(3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
%             T_CH(i_CH,2) = 1.69e8*A*(2*sigma(i_CH,2)*sqrt(log(2))/lambda0)^2;
%         end
%         offset(k,2) = (shift((k-1)*4+1,2) + shift((k-1)*4+2,2))/2;%対向視線から得られたオフセット[nm]
%         if show_offset
%             disp(offset(k,2))
%         end
%         %オフセットを引いてガウスフィッティング
%         for i = 1:4
%             i_CH = (k-1)*4+i;%CH番号
%             S2 = [L1_shaped(:,i_CH) spectrum2(:,i_CH)]; %[波長,強度]
%             s2 = size(S1);%S1の[行数,列数]
%             MAX2 = max(spectrum2(:,i_CH)); %スペクトルの最大値
%             j = 1;
%             while j < s2(1)+1 %SNの悪いデータを除く
%                 if S2(j,2) < MAX2*0.5
%                     S2(j,:) = [];
%                 else
%                     j = j+1;
%                 end
%                 s2 = size(S2);
%             end
%             f = fit(S2(:,1),S2(:,2),'gauss1');
%             coef=coeffvalues(f);
%             amp(i_CH,2) = coef(1);
%             shift(i_CH,2) = coef(2)-lambda0;
%             sigma(i_CH,2) = sqrt(coef(3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
%             T_CH(i_CH,2) = 1.69e8*A*(2*sigma(i_CH,2)*sqrt(log(2))/lambda0)^2;
%             if plot_fit
%                 subplot(mpoints.n_r,4,i_CH);
%                 plot(f,S2(:,1),S2(:,2));
%                 xline(lambda0);
%                 title([num2str(T_CH(i_CH,2)),' eV'])
%                 legend('off')
%                 xlabel('Wavelength [nm]')
%                 ylabel('Intensity [cnt]')
%             end
%         end
%     end
%     if plot_fit
%         sgtitle('Fitting data1 (Horizontal：Channel Vertical：Position)')
%     end
%
%     %data2の流速、温度を計算
%     for i = 1:mpoints.n_r
%         set = (i-1)*4;
%         Va = -shift(set+4,2)/lambda0*Vc;
%         Vb = -shift(set+3,2)/lambda0*Vc;
%         V_i(i,3) = (Va-Vb)/(2*cos(Angle*pi/180));%Vz
%         V_i(i,4) = -(Va+Vb)/(2*sin(Angle*pi/180));%Vr
%         absV(i,2) = sqrt(V_i(i,3)^2 + V_i(i,4)^2);
%         %        fprintf('(z, r) = (%.2f, %.2f) cm\n (Vz, Vr) = (%.2f, %.2f) km/s\n'...
%         %             ,z(i,2),r(i,2),V(i,3),V(i,4));
%         %             fprintf('T_i: %.2f, %.2f, %.2f, %.2f eV\n'...
%         %                     ,T(set+1,2),T(set+2,2),T(set+3,2),T(set+4,2));
%         T_i1=(T_CH(set+1,2)+T_CH(set+2,2)+T_CH(set+3,2)+T_CH(set+4,2))/4;
%         fprintf('%.2f\n',T_i1);
%     end
% end

%流速データを保存
if save_flow
    if not(exist([pathname.flowdata,'/',num2str(date)],'dir'))
        mkdir(sprintf("%s", pathname.flowdata), sprintf("%s", num2str(date)));
    end
    save([pathname.flowdata,'/',num2str(date),'/shot',num2str(ICCD.shot),'_',num2str(ICCD.trg),'us_w=',num2str(ICCD.exp_w),'_gain=',num2str(ICCD.gain),'.mat'],'V_i','absV','T_i')
end