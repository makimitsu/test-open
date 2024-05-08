xfunction [V_i,absV,T_i] = cal_ionvdist(date,expval,ICCD,mpoints,pathname,show_offset,plot_spectra,inversion_method,plot_analisis,plot_vdist,plot_type,save_fig,plot_compare,save_vdist,Ti_type)
%------フィッティング、再構成調整用【input】-------
l_mov = 7;%【input】波長方向移動平均長さ
th_ratio = 0.8;%【input】フィッティング時閾値
cut_l_L = 81;%【input】波長軸の切り取り長さ(奇数)
cut_d_L = 30;%【input】波長方向位置ずれ調整
cut_l_CH = 9;%【input】CH方向の切り取り長さ(奇数)
cut_d_CH = 2;%【input】CH方向位置ずれ調整
n_L = 51;%【input】再構成用波長軸配列の長さ(奇数, cut_l_L以下にする)

%物理定数
Vc = 299792.458;%光速(km/s)
mp = 1.67e-27;%陽子質量(kg)
kB = 1.60e-19;%ボルツマン定数(J/eV)
%装置変数
Angle = [0 30 150];%視線角度[度](0~180)
Theta = Angle*pi/180;%視線角度[rad]に変換
n_Theta = numel(Theta);%視線角度数
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

dir = [pathname.IDSP,'/',num2str(date)];%ディレクトリ1
filename = [dir,'/shot',num2str(ICCD.shot),'_',num2str(ICCD.trg),'us_w=',num2str(ICCD.exp_w),'_gain=',num2str(ICCD.gain),'.asc'];%ICCDファイル名
if not(exist(filename,"file"))
    warning([filename,' does not exist.']);
    V_i = char.empty;
    absV = char.empty;
    T_i = char.empty;
    return
end
%スペクトルを取得
data = importdata(filename);

%軸の定義
X = data(:,1);%X1(ピクセル)軸を定義(1~1024)
[l_X,~]=size(data);%X1軸の長さを取得
L = zeros(l_X,mpoints.n_CH);%λ軸(チャンネルにより変わる)
L_shaped = zeros(cut_l_L,mpoints.n_CH);%λ軸(チャンネルにより変わる)
px2nm = zeros(mpoints.n_CH,1);%nm/pixel

switch ICCD.line
    case 'Ar'%アルゴンの時
        for i = 1:mpoints.n_CH
            px2nm(i,1) = center(i,4);
            L(:,i) = lambda1 - px2nm(i,1)*(X(:,1)-center(i,3));
            L_shaped(:,i) = L(round(center(i,3))-(cut_l_L-1)/2+cut_d_L:round(center(i,3))+(cut_l_L-1)/2+cut_d_L,i);
        end
        d_L = mean(px2nm, 'all');%Δλ
    case 'H'%水素の時
        % for i = 1:mpoints.n_CH
        %     px2nm(i,1) = 0.00536;
        %     L1(:,i) = px2nm(i,1)*(X1-centerX(i))+lambda0;
        % end
        % d_L = 0.00536;%Δλ
    otherwise
        warning('Input error in ICCD.line.')%ICCD.lineの入力エラー
        return;
end

n_Vx = n_L;
d_Vx = d_L/lambda0*Vc;%ΔVx
Vx = linspace(-d_Vx*(n_Vx-1)/2,d_Vx*(n_Vx-1)/2,n_Vx);%Vx軸
n_Vy = n_L;
d_Vy = d_L/lambda0*Vc;%ΔVy
Vy = linspace(-d_Vy*(n_Vy-1)/2,d_Vy*(n_Vy-1)/2,n_Vy);%Vy軸

%分光データを整形
spectrum = zeros(cut_l_L,mpoints.n_CH);%data1の分光結果を入れる
for i = 1:mpoints.n_CH
    spectrum(:,i) = sum(data(round(center(i,3))-(cut_l_L-1)/2+cut_d_L:round(center(i,3))+(cut_l_L-1)/2+cut_d_L, ...
        round(center(i,2))-(cut_l_CH-1)/2+cut_d_CH:round(center(i,2))+(cut_l_CH-1)/2+cut_d_CH),2);
end

%波長方向の移動平均から離れた外れ値を移動平均に変更(ノイズ除去)
mov_spectrum = movmean(spectrum,l_mov);
for i = 1:mpoints.n_CH
    for j = 1:cut_l_L
        if spectrum(j,i) > mov_spectrum(j,i)*1.1
            spectrum(j,i) = mov_spectrum(j,i);
        elseif spectrum(j,i) < mov_spectrum(j,i)*0.9
            spectrum(j,i) = mov_spectrum(j,i);
        end
    end
end

shift = zeros(mpoints.n_CH,1);
offset = zeros(mpoints.n_r,1);
error_CH = zeros(mpoints.n_CH,1);

%-----対抗視線のフィッティングを行うことでオフセットを求め、λ軸を較正------
for k = 1:mpoints.n_r
    %0度ペアスペクトルからオフセットを検出
    for i = 1:2
        i_CH = (k-1)*4+i;%CH番号
        S = [L_shaped(:,i_CH) spectrum(:,i_CH)]; %[波長,強度]
        [S_rows,~] = size(S);
        S_max = max(spectrum(:,i_CH)); %スペクトルの最大値
        j = 1;
        while j < S_rows+1 %SNの悪いデータを除く
            if S(j,2) < S_max*th_ratio
                S(j,:) = [];
            else
                j = j+1;
            end
            [S_rows,~] = size(S);
        end
        try
            f = fit(S(:,1),S(:,2),'gauss1');
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
    if show_offset
        disp(['Offset = ',num2str(offset(k,1)),' in Row ',num2str(k),'.'])
    end
    for i = 1:4
        i_CH = (k-1)*4+i;%CH番号
        L_shaped(:,i_CH) = L_shaped(:,i_CH)-offset(k,1);
    end
end

%-----------------再構成できるデータ形式に変換---------------------
%V=0付近のデータを抽出
Lambda = zeros(n_L,mpoints.n_CH);%再構成用波長データ
spectrum_trans = zeros(n_L,mpoints.n_CH);%再構成用スペクトルデータ
for i = 1:mpoints.n_CH
    IDX = knnsearch(L_shaped(:,i),lambda0);%L1の中で最もlambda0に近いセル番号を取得
    Lambda(:,i) = L_shaped(IDX-(n_L-1)/2:IDX+(n_L-1)/2,i)-lambda0;
    spectrum_trans(:,i) = spectrum(IDX-(n_L-1)/2:IDX+(n_L-1)/2,i);
end
%データを規格化
f_st = 100;
for i = 1:mpoints.n_CH
    spectrum_trans(:,i) = spectrum_trans(:,i)/sum(spectrum_trans(:,i))*f_st;
end

%観測スペクトルを描画
if plot_spectra
    figure('Position',[300 50 1000 1000],'visible','on')
    sgtitle(['Fitting data (Horizontal：View Line, Vertical：Measured Position)',newline, ...
        'shot',num2str(ICCD.shot),'-',num2str(ICCD.trg),'us-w=',num2str(ICCD.exp_w),'-gain=',num2str(ICCD.gain),'.asc'])
    for k = 1:mpoints.n_r
        if show_offset
            disp(['Offset = ',num2str(offset(k,1)),' in Row ',num2str(k),'.'])
        end
        %オフセットを引いてガウスフィッティング
        for i = 1:4
            i_CH = (k-1)*4+i;%CH番号
            S = [Lambda(:,i_CH) spectrum_trans(:,i_CH)]; %[波長,強度]
            try
                f = fit(S(:,1),S(:,2),'gauss1');
                subplot(mpoints.n_r,4,i_CH);
                plot(f,S(:,1),S(:,2));
                plot(S(:,1),S(:,2));
                title(['CH ',num2str(i_CH)])
                legend('off')
                xlabel('Shift Wavelength [nm]')
                ylabel('Intensity [cnt]')
            catch ME
                warning(['Fitting failed in CH ',num2str(i_CH),'.']);
                break
            end
        end
    end
end

%1次元配列P,行列W,1次元配列Fを計測点ごとに用意・・・P(:,k) = W(:,:,k)×F(:,k)
P = zeros(n_Theta*n_L,mpoints.n_r);%P(:,k)がk番目の計測点の分光データ
W = zeros(n_Theta*n_L,n_Vx*n_Vy,mpoints.n_r);%W(:,:,k)がk番目のF->P変換行列
F = zeros(n_Vx*n_Vy,mpoints.n_r);%F(:,k)がk番目の計測点の再構成された速度分布
ppoints = zeros(n_Theta*n_L,2,mpoints.n_r);%ppoints(:,:,k)がk番目のプロット点
for k = 1:mpoints.n_r
    P(:,k) = cat(1, spectrum_trans(:,4*(k-1)+2),spectrum_trans(:,4*(k-1)+4),spectrum_trans(:,4*(k-1)+3));%0度,30度,150度を連結
    for i = 1:n_Theta*n_L
        %W(i,:,k)の対応シフト波長番号lを取得：対応シフト波長はL1_trans(l,4*(k-1)+t)
        l = int8(mod(i,n_L));
        if l == 0
            l = n_L;
        end
        %W(i,:,k)の対応視線角度番号tを取得：対応角度はTheta(t)
        t = 1 + idivide(i-1, int8(n_L), 'floor');
        %視線方向単位ベクトルv_thetaを計算
        v_theta = [cos(Theta(t)) sin(Theta(t))];
        %速度空間上の観測者座標を計算
        ppoints(i,1,k) = Lambda(l,4*(k-1)+t)/lambda0*Vc*v_theta(1);
        ppoints(i,2,k) = Lambda(l,4*(k-1)+t)/lambda0*Vc*v_theta(2);
        for j = 1:n_Vx*n_Vy
            %W(:,j)の対応速度番号x,yを取得：対応速度はVx(x),Vy(y)
            x = 1 + idivide(j-1, int8(n_Vy), 'floor');
            y = int8(mod(j,n_Vy));
            if y == 0
                y = n_Vy;
            end
            V = [Vx(x), Vy(y)];
            D = dot(V, v_theta);%視線方向速度
            %W(i,j,k) = 1 - abs(-Lambda(l,4*(k-1)+t)/lambda0*Vc-D)*(1/(d_L/lambda0*Vc));%一次関数フィルタ
            W(i,j,k) = 1 - ((-Lambda(l,4*(k-1)+t)/lambda0*Vc-D)*(1/(d_L/lambda0*Vc)))^2;%二次関数フィルタ
            if W(i,j,k) < 0
                W(i,j,k) = 0;
            end
        end
    end
    %---------------Invertion theory----------------
    switch inversion_method
        case 1 %pinvを使ってW^(-1)を求める。
            F(:,k) = method_pseudo(W(:,:,k),P(:,k));
        case 2 %Tikhonov 0th
            F(:,k) = method_Tikhonov0(W(:,:,k),P(:,k),plot_analisis);
        case 3 %Tikhonov 1st
            [F(:,k),~] = method_Tikhonov1(W(:,:,k),P(:,k),n_Vx,n_Vy,d_Vx,d_Vy,plot_analisis);
        case 4 %Tikhonov 2nd
            F(:,k) = method_Tikhonov2(W(:,:,k),P(:,k),n_Vx,n_Vy,d_Vx,d_Vy,plot_analisis);
        case 5 %非負制約SIRT
            F(:,k) = method_SIRT(W(:,:,k),P(:,k),plot_analisis);
        case 6 %minimum Fischer information
            F(:,k) = method_MFI(W(:,:,k),P(:,k),n_Vx,n_Vy,d_Vx,d_Vy,plot_analisis);
        otherwise
            warning('Input error in inversion_method.')%inversion_methodの入力エラー
            return;
    end
end

%Fを二次元画像形式draw_Fに変換
draw_F = zeros(n_Vx,n_Vy,mpoints.n_r);%draw_F(:,:,k)がk番目の二次元画像
summit_x = zeros(mpoints.n_r,1);%summit_x(k,1)がk番目の二次元画像の最大値のx座標
summit_y = zeros(mpoints.n_r,1);%summit_y(k,1)がk番目の二次元画像の最大値のy座標
for k = 1:mpoints.n_r
    draw_F(:,:,k) = reshape(F(:,k), [n_Vx,n_Vy]);
    draw_F(:,:,k) = flipud(draw_F(:,:,k));
    buf_F = draw_F(:,:,k);
    [summit_y(k,1),summit_x(k,1)] = find(buf_F == max(buf_F(:)));
    % fprintf('Summit is (Vz = %.1f, Vr = %.1f)[km/s].\n',Vx(summit_x(k,1)),Vy(summit_y(k,1)));
end

%流速を計算
V_i = zeros(mpoints.n_r,2);%1列目・Vz(km/s)、2列目・Vr(km/s)
absV = zeros(mpoints.n_r,1);%|V|(km/s)
for k = 1:mpoints.n_r
    sum_F_flow = 0;
    for i = 1:n_Vy
        for j = 1:n_Vx
            if draw_F(i,j,k)> draw_F(summit_y(k,1),summit_x(k,1),k)*0.5
                V_i(k,1) = V_i(k,1) + draw_F(i,j,k)*Vx(j);
                V_i(k,2) = V_i(k,2) + draw_F(i,j,k)*Vy(i);
                sum_F_flow = sum_F_flow + draw_F(i,j,k);
            end
        end
    end
    V_i(k,1) = V_i(k,1)/sum_F_flow;%Vz
    V_i(k,2) = V_i(k,2)/sum_F_flow;%Vr
    absV(k,1) = sqrt(V_i(k,1)^2 + V_i(k,2)^2);%|V|
    % fprintf('(z = %.1f, r = %.1f)[cm]: (Vz = %.1f, Vr = %.1f, |Vi| = %.1f)[km/s]\n', ...
    %     mpoints.z(k,1),mpoints.r(k,1),V_i(k,1),V_i(k,2),absV(k,1));
end

%温度を計算
T_i_z = zeros(mpoints.n_r,1);%r方向温度(eV)
T_i_r = zeros(mpoints.n_r,1);%z方向温度(eV)
switch Ti_type
    case 'dispersion'%擬似温度を計算(分散から計算)
        for k = 1:mpoints.n_r
            sum_F_temp = 0;
            disp_V = zeros(2,1);%1列目Vz分散,2列目Vr分散
            for i = 1:n_Vy
                for j = 1:n_Vx
                    if draw_F(i,j,k)> draw_F(summit_y(k,1),summit_x(k,1))*0
                        disp_V(1,1) = disp_V(1,1) + draw_F(i,j,k)*(Vx(j)-V_i(k,1))^2;
                        disp_V(2,1) = disp_V(2,1) + draw_F(i,j,k)*(Vy(i)-V_i(k,2))^2;
                        sum_F_temp = sum_F_temp + draw_F(i,j,k);
                    end
                end
            end
            % %装置関数を考慮
            % disp_V(1,1) = disp_V(1,1)/sum_F_temp - (mean(center(1+4*(k-1):4+4*(k-1),5),'all')*d_L/lambda0*Vc)^2;
            % disp_V(2,1) = disp_V(2,1)/sum_F_temp - (mean(center(1+4*(k-1):4+4*(k-1),5),'all')*d_L/lambda0*Vc)^2;
            %装置関数を考慮しない
            disp_V(1,1) = disp_V(1,1)/sum_F_temp;
            disp_V(2,1) = disp_V(2,1)/sum_F_temp;
            T_i_z(k,1) = disp_V(1,1)*10^6*A*mp/(2*kB);
            T_i_r(k,1) = disp_V(2,1)*10^6*A*mp/(2*kB);
            % fprintf('(z = %.1f, r = %.1f)[cm]: (Tz = %.1f, Tr = %.1f)[eV]\n', ...
            %     mpoints.z(k,1),mpoints.r(k,1),T_i_z(k,1),T_i_z(k,1));
        end
    case 'FWHM'%擬似温度を計算(半値幅から計算)
        figure('Position',[300 50 600 1200],'visible','on')
        sgtitle(['Full Width Half Maximum (Vertical：Measured Position)',newline, ...
            'shot',num2str(ICCD.shot),'-',num2str(ICCD.trg),'us-w=',num2str(ICCD.exp_w),'-gain=',num2str(ICCD.gain),'.asc'])
        for k = 1:mpoints.n_r
            %Z方向FWHMを計算
            W_0 = zeros(n_L,n_Vx*n_Vy);%0度スペクトルを計算
            for i = 1:n_L
                Unit_0 = [cos(0) sin(0)];
                for j = 1:n_Vx*n_Vy
                    %W(:,j)の対応速度番号x,yを取得：対応速度はVx(x),Vy(y)
                    x = 1 + idivide(j-1, int8(n_Vy), 'floor');
                    y = int8(mod(j,n_Vy));
                    if y == 0
                        y = n_L;
                    end
                    V = [Vx(x), Vy(y)];
                    D = dot(V, Unit_0);%視線方向速度
                    %      W_0(i,j) = 1 - abs(Lambda(i,4*(k-1)+1)/lambda0*Vc-D)*(1/(d_L/lambda0*Vc));%一次関数フィルタ
                    W_0(i,j) = 1 - ((Lambda(i,4*(k-1)+1)/lambda0*Vc-D)*(1/(d_L/lambda0*Vc)))^2;%二次関数フィルタ
                    if W_0(i,j) < 0
                        W_0(i,j) = 0;
                    end
                end
            end
            P_0 = W_0*F(:,k);
            summit_P_0 = find(P_0 == max(P_0));
            %プロット
            subplot(mpoints.n_r,2,2*(k-1)+1);
            plot(Vx,P_0/max(P_0), 'b','LineWidth',2)
            hold on
            yline(0.5, '--r','LineWidth',2);
            title(['r = ',num2str(mpoints.r(k)),'[cm]'])
            xlabel('V_{Z} [km/s]')
            ylabel('Intensity')
            % xlim([-60 60])
            % xticks(-60:30:60)
            yticks(0:0.5:1)
            ax = gca;
            ax.FontSize = 10;
            hold off
            for i = 1:n_L - summit_P_0
                if P_0(i+summit_P_0) < max(P_0)*0.5
                    I_right_x = i+summit_P_0;
                    break
                end
            end
            fwhm_right_x = Vx(I_right_x);%半値x右端
            for i = 1:summit_P_0
                if P_0(summit_P_0-i+1) < max(P_0)*0.5
                    I_left_x = summit_P_0-i+1;
                    break
                end
            end
            fwhm_left_x = Vx(I_left_x);%半値x左端
            fwhm_x = fwhm_right_x - fwhm_left_x;
            sigma_x = fwhm_x/(2*sqrt(2*log(2)));%FWHM = 2σ√2ln(2)
            %R方向FWHMを計算
            W_90 = zeros(n_L,n_Vx*n_Vy);%0度スペクトルを計算
            for i = 1:n_L
                Unit_90 = [cos(pi/2) sin(pi/2)];
                for j = 1:n_Vx*n_Vy
                    %W(:,j)の対応速度番号x,yを取得：対応速度はVx(x),Vy(y)
                    x = 1 + idivide(j-1, int8(n_Vy), 'floor');
                    y = int8(mod(j,n_Vy));
                    if y == 0
                        y = n_L;
                    end
                    V = [Vx(x), Vy(y)];
                    D = dot(V, Unit_90);%視線方向速度
                    %      W_90(i,j) = 1 - abs(Lambda(i,4*(k-1)+1)/lambda0*Vc-D)*(1/(d_L/lambda0*Vc));%一次関数フィルタ
                    W_90(i,j) = 1 - ((Lambda(i,4*(k-1)+1)/lambda0*Vc-D)*(1/(d_L/lambda0*Vc)))^2;%二次関数フィルタ
                    if W_90(i,j) < 0
                        W_90(i,j) = 0;
                    end
                end
            end
            P_90 = W_90*F(:,k);
            summit_P_90 = find(P_90 == max(P_90));
            %プロット
            subplot(mpoints.n_r,2,2*(k-1)+2);
            plot(Vy,P_90/max(P_90), 'b','LineWidth',2)
            hold on
            yline(0.5, '--r','LineWidth',2);
            title(['r = ',num2str(mpoints.r(k)),'[cm]'])
            xlabel('V_{R} [km/s]')
            ylabel('Intensity')
            % xlim([-60 60])
            % xticks(-60:20:60)
            % yticks(0:0.1:1)
            ax = gca;
            ax.FontSize = 10;
            for i = 1:n_L - summit_P_90
                if P_90(i+summit_P_90) < max(P_90)*0.5
                    I_right_y = i+summit_P_90;
                    break
                end
            end
            fwhm_right_y = Vy(I_right_y);%半値y右端
            for i = 1:summit_P_90
                if P_0(summit_P_90-i+1) < max(P_90)*0.5
                    I_left_y = summit_P_90-i+1;
                    break
                end
            end
            fwhm_left_y = Vy(I_left_y);%半値x左端
            fwhm_y = fwhm_right_y - fwhm_left_y;
            sigma_y = fwhm_y/(2*sqrt(2*log(2)));%FWHM = 2σ√2ln(2)
            % %装置関数を考慮
            % disp_V(1,1) = sigma_x^2 - (mean(center(1+4*(k-1):4+4*(k-1),5),'all')*d_L/lambda0*Vc)^2;
            % disp_V(2,1) = sigma_y^2 - (mean(center(1+4*(k-1):4+4*(k-1),5),'all')*d_L/lambda0*Vc)^2;
            %装置関数を考慮しない
            disp_V(1,1) = sigma_x^2;
            disp_V(2,1) = sigma_y^2;
            T_i_z(k,1) = disp_V(1,1)*10^6*A*mp/(2*kB);
            T_i_r(k,1) = disp_V(2,1)*10^6*A*mp/(2*kB);
            % fprintf('(z = %.1f, r = %.1f)[cm]: (Tz = %.1f, Tr = %.1f)[eV]\n', ...
            %     mpoints.z(k,1),mpoints.r(k,1),T_i_z(k,1),T_i_z(k,1));
        end
    otherwise
        warning('Input error in Ti_type.')%Ti_typeの入力エラー
        return;
end
T_i = T_i_z;

%速度分布をプロット
if plot_vdist
    plot_ionvdist(Vx,Vy,F,date,expval,ICCD,pathname,mpoints,ppoints,plot_type,save_fig)
end

%再構成結果から得られるスペクトルを計算、スペクトルデータと比較
if plot_compare
    plot_inversion_compare(F,W,P,Lambda,mpoints,Angle,ICCD)
end

%速度分布データを保存
if save_vdist
    if not(exist([pathname.vdistdata,'/',num2str(date)],'dir'))
        mkdir(sprintf("%s", pathname.vdistdata), sprintf("%s", num2str(date)));
    end
    save([pathname.vdistdata,'/',num2str(date),'/shot',num2str(ICCD.shot),'_', ...
        num2str(ICCD.trg),'us_w=',num2str(ICCD.exp_w),'_gain=',num2str(ICCD.gain),'.mat'], ...
        'V_i','absV','T_i','F','W','P','Lambda','Vx','Vy','ppoints','Angle')
end
end