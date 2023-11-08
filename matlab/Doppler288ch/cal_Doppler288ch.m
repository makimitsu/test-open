function [z_IDS,r_IDS,Ti_IDS,Ti_max_IDS,Ti_min_IDS,Em_IDS] = cal_Doppler288ch(date,IDS288ch,pathname,cal_LineInt,plot_2D_spectra,save_Ti)
%-----------------------解析オプション【input】----------------------------
plot_ICCD = false;%【input】ICCD画像をプロット
plot_CH_spectra = false;%【input】CHごとのスペクトルをプロット
plot_LineInt_interp = false;%【input】死んだCHの補間Ti,Emをプロット(cal_LineInt = trueが必要)
plot_LineInt_result = false;%【input】線積分イオン温度、発光強度分布をプロット(cal_LineInt = trueが必要)
plot_2D_interp = false;%【input】補間スペクトルをプロット

hw_lambda = 50;%【input】波長切り出し半幅
hw_ch = 5;%【input】CH方向切り出し半幅
hw_fit = 12;%【input】フィッティング波長切り出し半幅(< hw_lambda)
num_r = 30;%【input】r分割数(比例して計算時間が増える)
hw_lambdaA = 40;%【input】lambdaA半幅(< hw_lambda)
%------------------------------------------------------------------

%物理定数
switch IDS288ch.line
    case 'Ar'
        lambda0 = 480.6;%471.3;656.3;480.6;587.562;486.133;656.3;486.133;468.57;486.133;nm%線スペクトル波長
        mass = 39.95;%4.;12.01%イオン質量数
    otherwise
        warning('Input error in ICCD.line.')%ICCD.lineの入力エラー
        z_IDS = char.empty;
        r_IDS = char.empty;
        Ti_IDS = char.empty;
        Ti_max_IDS = char.empty;
        Ti_min_IDS = char.empty;
        Em_IDS = char.empty;
        return;
end

%------校正値(1つのファイルにまとめたい)---------------------------------
separation = [0 15 31 45 57 73 89 105 120 136 152 164 180 196 211 226 240 255];%CHをZ方向で切り分けるための値
resolution = -1.420387E-12*lambda0^3 - 2.156031E-09*lambda0^2 + 1.250038E-06*lambda0 + 3.830769E-03;0.0037714;...
    %-0.000000000001420387*lambda0^3 - 0.000000002156031*lambda0^2 + 0.000001250038*lambda0 + 0.003830769
z = importdata("z_negative.txt")*1e-3;%計測視線Z[m]
p = importdata("r.txt")*1e-3;%計測視線と中心軸の距離P[m]
edge = 0.33;%Rの最大値[m](2022/7　解析～　ポテンシャルの範囲より仮定)
calib = importdata("Ar_calibration.0916_remake.txt");%ICCD校正ファイル
%------------------------------------------------------------------

%データ読み込み(GUI)
dir = [pathname.IDS288ch,'/20',num2str(date)];%実験日ディレクトリ
[file,path] = uigetfile('*.asc','Select a file',dir);
if isequal(file,0)
    disp('User selected Cancel');
    z_IDS = char.empty;
    r_IDS = char.empty;
    Ti_IDS = char.empty;
    Ti_max_IDS = char.empty;
    Ti_min_IDS = char.empty;
    Em_IDS = char.empty;
    return
else
    disp(['User selected ', fullfile(path,file)]);
    filename = fullfile(path,file);
end
% %データ読み込み(手入力)
% filename = [dir,'/shot3_PF39_TF-2_EF120_gas0.62_delay468_width3.asc'];
%スペクトルを取得

data = importdata(filename);
data = data(:,2:1025);%1列目は分光データではないので削除

%校正ファイル読み込み
ch = calib.data(:,1);
center = calib.data(:,2);
smile = calib.data(:,3);
relative = calib.data(:,5);
instru = calib.data(:,6);
Ti_instru_CH = 1.69e8*mass*(2.*resolution*instru*sqrt(2.*log(2.))/lambda0).^2;

%各CHの波長軸を生成
lambda = zeros(2*hw_lambda+1,numel(ch));
idx_l0 = zeros(numel(ch),1);
px = transpose(linspace(1,1024,1024));
for i=1:numel(ch)
    lx=(px-1-smile(i))*resolution+lambda0;%+0.13
    idx_l0(i) = knnsearch(lx,lambda0);%lambdaの中で最もlambda0に近いセル番号を取得
    lambda(:,i) = lx(idx_l0(i)-hw_lambda:idx_l0(i)+hw_lambda);
end

%ICCD生画像を描画(目視で確認用)
if plot_ICCD
    figure('Position',[0 300 400 550])
    subplot(2,1,1)
    contour(px,px,data')
    title('ICCD Raw Image')
    xlabel('X　(lambda) [px]')
    ylabel('Y (position) [px]')
    hold on
    for i=1:numel(ch)
        hol_X = linspace(smile(i)-hw_lambda,smile(i)+hw_lambda);
        hol_Y = center(i)*ones(100,1);
        plot(hol_X,hol_Y,'r')
        hold on
        ver_X = smile(i)*ones(100,1);
        ver_Y = linspace(center(i)-hw_ch,center(i)+hw_ch);
        plot(ver_X,ver_Y,'r','LineWidth',1)
    end
    % ylim([900 1024])
    hold off
    %ICCD生画像をCH方向に足し合わせたもののフィッティング(目視で確認用)
    sum_data = sum(data,2);
    offset_sum = min(movmean(sum_data,10));
    sum_data = sum_data - offset_sum;
    f = fit(px,sum_data,'gauss5');
    coeff5=coeffvalues(f);
    % centerL = [coeff5(2),coeff5(5),coeff5(8),coeff5(11),coeff5(14)];
    % centerL = round(sort(centerL));
    subplot(2,1,2)
    plot(f,px,sum_data')
    title('ICCD Integrated Image')
    xlabel('X [px]')
    ylabel('Strength [a.u.]')
    xlim([1 1024])
    ylim([0 inf])
    hold off
end

%-------CHごとの線積分温度、線積分発光強度を計算--------
hw_plot_px = hw_lambda;%プロット範囲半幅
col_subp1 = 2;%サブプロット行数
raw_subp1 = 2;%サブプロット列数
n_subp1 = col_subp1*raw_subp1;
passive_Ti = zeros(numel(ch),1);%CH温度[eV]
passive_Timax = zeros(numel(ch),1);%CH温度95%信頼区間上限[eV]
passive_Timin = zeros(numel(ch),1);%CH温度95%信頼区間下限[eV]
passive_Em = zeros(numel(ch),1);%CH発光強度[a.u.]
spectra = zeros(2*hw_lambda+1,numel(ch));%CHスペクトル
for i=1:numel(ch)
    spectra(:,i) = sum(data(idx_l0(i)-hw_lambda:idx_l0(i)+hw_lambda,center(i)-hw_ch:center(i)+hw_ch),2)*relative(i);
    offset = min(movmean(spectra(:,i),20));
    spectra(:,i) = spectra(:,i) - offset;%オフセットを引く
    f = fit(lambda(hw_lambda+1-hw_fit:hw_lambda+1+hw_fit,i),spectra(hw_lambda+1-hw_fit:hw_lambda+1+hw_fit,i),'gauss1');
    coeff = coeffvalues(f);%フィッティング係数
    % coeff(3)/resolution/sqrt(2);%IDLのcoeff[2]に等しい(確認用)
    % Sigma(i) = coeff(3)/sqrt(2);
    confi = confint(f);%フィッティング係数の95%信頼区間(1行目下限, 2行目上限)
    passive_Ti(i) = 1.69e8*mass*(2*coeff(3)*sqrt(log(2))/lambda0)^2-Ti_instru_CH(i);
    passive_Timax(i) = 1.69e8*mass*(2*confi(2,3)*sqrt(log(2))/lambda0)^2-Ti_instru_CH(i);
    passive_Timin(i) = 1.69e8*mass*(2*confi(1,3)*sqrt(log(2))/lambda0)^2-Ti_instru_CH(i);
    passive_Em(i) = resolution*sum(spectra(:,i));
    if plot_CH_spectra
        idx_subp1 = mod(i-1,n_subp1)+1;%サブプロット位置番号
        if i == 1
            figure('Position',[0 500 400 300])
            sgtitle('Line Integrated Spectra')
        end
        subplot(col_subp1,raw_subp1,idx_subp1)
        if i <= n_subp1
            if i == 1
                fit_x = transpose(linspace(lambda0-hw_plot_px*resolution,lambda0+hw_plot_px*resolution,100));
                fit_y = feval(f,fit_x);
                pfit1 = plot(fit_x,fit_y,'r-',lambda(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i),spectra(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i),'b+');%フィティングあり
                pfit = repmat(pfit1,[1,1,n_subp1]);
            else
                fit_y = feval(f,fit_x);
                pfit(:,:,idx_subp1) = plot(fit_x,fit_y,'r-',lambda(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i),spectra(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i),'b+');%フィティングあり
            end
            xline(lambda0)
            title(['CH',num2str(ch(i))])
            xlabel('Wavelength [nm]')
            ylabel('Strength [a.u.]')
            xlim([lambda0-hw_plot_px*resolution lambda0+hw_plot_px*resolution])
            ylim([0 inf])
            legend('off')
        else
            fit_y = feval(f,fit_x);
            pfit(1,:,idx_subp1).YData = fit_y;
            pfit(2,:,idx_subp1).XData = lambda(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i);
            pfit(2,:,idx_subp1).YData = spectra(hw_lambda+1-hw_plot_px:hw_lambda+1+hw_plot_px,i);
            title(['CH',num2str(ch(i))])
            drawnow
        end
    end
end

%%----------線積分温度、線積分発光強度二次元分布を計算------------
if cal_LineInt
    Ti_LineInt = zeros(numel(p),numel(z));%二次元温度[eV]
    Em_LineInt = zeros(numel(p),numel(z));%二次元発光強度[a.u.]
    %死んだCHを補間
    for i = 1:numel(z)
        if i < numel(z)
            Ti_LineInt(:,i) = pchip(ch(separation(i)+1:separation(i+1))-ch(separation(i)+1)+1,passive_Ti(separation(i)+1:separation(i+1)),linspace(1,numel(p),numel(p)));
            Em_LineInt(:,i) = pchip(ch(separation(i)+1:separation(i+1))-ch(separation(i)+1)+1,passive_Em(separation(i)+1:separation(i+1)),linspace(1,numel(p),numel(p)));
        else
            Ti_LineInt(:,i) = pchip(ch(separation(i)+1:end)-ch(separation(i)+1)+1,passive_Ti(separation(i)+1:end),linspace(1,numel(p),numel(p)));
            Em_LineInt(:,i) = pchip(ch(separation(i)+1:end)-ch(separation(i)+1)+1,passive_Em(separation(i)+1:end),linspace(1,numel(p),numel(p)));
        end
    end
    %補間結果をプロット
    if plot_LineInt_interp
        figure('Position',[1000 1000 600 200])
        tiledlayout(1,2)
        nexttile;
        p1 = plot(linspace(1,numel(p),numel(p)),Ti_LineInt(:,1),'ro-');
        hold on
        p2 = plot(ch(separation(1)+1:separation(1+1))-ch(separation(1)+1)+1,passive_Ti(separation(1)+1:separation(1+1)),'bo-');
        xlabel('CH')
        ylabel('Ti [eV]')
        xlim([1 numel(p)])
        title('Interpoltaion of Ti')
        legend(p1,'Interpolated')
        nexttile;
        p3 = plot(linspace(1,numel(p),numel(p)),Em_LineInt(:,1),'ro-');
        hold on
        p4 = plot(ch(separation(1)+1:separation(1+1))-ch(separation(1)+1)+1,passive_Em(separation(1)+1:separation(1+1)),'bo-');
        xlabel('CH')
        ylabel('Emission [a.u.]')
        xlim([1 numel(p)])
        title('Interpoltaion of Em')
        legend(p3,'Interpolated')
        hold off
        for i = 1:numel(z)
            sgtitle([num2str(i),'/',num2str(numel(z))])
            p1.YData = Ti_LineInt(:,i);
            p3.YData = Em_LineInt(:,i);
            if i < numel(z)
                p2.XData = ch(separation(i)+1:separation(i+1))-ch(separation(i)+1)+1;
                p2.YData = passive_Ti(separation(i)+1:separation(i+1));
                p4.XData = p2.XData;
                p4.YData = passive_Em(separation(i)+1:separation(i+1));
                drawnow
            else
                p2.XData = ch(separation(i)+1:end)-ch(separation(i)+1)+1;
                p2.YData = passive_Ti(separation(i)+1:end);
                p4.XData = p2.XData;
                p4.YData = passive_Em(separation(i)+1:end);
                drawnow
            end
        end
    end
    %偶数列を上下反転(CHの順序が偶数列で反転するため)
    for i = 1:numel(z)
        if mod(i,2) == 0
            Ti_LineInt(:,i) = flipud(Ti_LineInt(:,i));
            Em_LineInt(:,i) = flipud(Em_LineInt(:,i));
        end
    end
    if plot_LineInt_result
        figure('Position',[0 0 700 350])
        tiledlayout(1,2)
        ax1 = nexttile;
        [~,h] = contourf(z,p,Ti_LineInt,100);
        h.LineStyle = 'none';
        daspect([1 1 1])
        colormap(ax1,jet)
        title('Line Integrated Ion Temperature')
        xlabel('Z [m]')
        ylabel('P [m]')
        c1 = colorbar;
        c1.Label.String = 'Ion Temperature [eV]';
        ax2 = nexttile;
        [~,h] = contourf(z,p,Em_LineInt,100);
        h.LineStyle = 'none';
        daspect([1 1 1])
        colormap(ax2,pink)
        title('Line Integrated Ion Emission')
        xlabel('Z [m]')
        ylabel('P [m]')
        c2 = colorbar;
        c2.Label.String = 'Ion Emission [a.u.]';
    end
end


%%----------アーベル変換線積分温度、発光二次元分布を計算------------
%三角グリッドの補間によりアーベル変換に用いる2次元(λ,P)スペクトルを用意
r = transpose(linspace(min(p),edge,num_r));
dr = r(2) - r(1);
lambdaA = linspace(-hw_lambdaA*resolution,hw_lambdaA*resolution,2*hw_lambdaA+1) + lambda0;
lambdaA = transpose(lambdaA);
spectra_interp = zeros(numel(lambdaA),numel(r),numel(z));
for i = 1:numel(z)
    if i < numel(z)
        tri_x = zeros(2*hw_lambdaA+1,separation(i+1)-separation(i)+1);
        if mod(i-1,2) == 0
            for k = 1:separation(i+1)-separation(i)
                idx_ch = separation(i)+k;
                tri_x(:,k) = lambda(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,idx_ch);
            end
            tri_x(:,end) = tri_x(:,end-1);
            tri_y = [p(ch(separation(i)+1:separation(i+1))-ch(separation(i)+1)+1);edge];
            tri_z = [spectra(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,separation(i)+1:separation(i+1)),zeros(numel(lambdaA),1)];
        else
            for k = 1:separation(i+1)-separation(i)
                idx_ch = separation(i)+k;
                tri_x(:,k+1) = lambda(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,idx_ch);
            end
            tri_x(:,1) = tri_x(:,2);
            tri_y = [edge;flipud(p(ch(separation(i)+1:separation(i+1))-ch(separation(i)+1)+1))];
            tri_z = [zeros(numel(lambdaA),1),spectra(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,separation(i)+1:separation(i+1))];
        end
    else
        tri_x = zeros(2*hw_lambdaA+1,numel(ch)-separation(i)+1);
        if mod(i-1,2) == 0
            for k = 1:numel(ch)-separation(i)
                idx_ch = separation(i)+k;
                tri_x(:,k) = lambda(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,idx_ch);
            end
            tri_x(:,end) = tri_x(:,end-1);
            tri_y = [p(ch(separation(i)+1:end)-ch(separation(i)+1)+1);edge];
            tri_z = [spectra(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,separation(i)+1:end),zeros(numel(lambdaA),1)];
        else
            for k = 1:numel(ch)-separation(i)
                idx_ch = separation(i)+k;
                tri_x(:,k+1) = lambda(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,idx_ch);
            end
            tri_x(:,1) = tri_x(:,2);
            tri_y = [edge;flipud(p(ch(separation(i)+1:end)-ch(separation(i)+1)+1))];
            tri_z = [zeros(numel(lambdaA),1),spectra(hw_lambda+1-hw_lambdaA:hw_lambda+1+hw_lambdaA,separation(i)+1:end)];
        end
    end
    buf_x = zeros(numel(tri_x),1);%CHごとの波長(CH数*ピクセル数)
    buf_y = zeros(numel(tri_x),1);%P
    buf_z = zeros(numel(tri_x),1);%スペクトル
    for j=1:numel(tri_x)
        index = mod(j-1,numel(tri_y));
        buf_x(j) = tri_x(ceil(j/numel(tri_y)),index+1);
        buf_y(j) = tri_y(index+1);
        buf_z(j) = tri_z(ceil(j/numel(tri_y)),index+1);
    end
    [grid_x, grid_y] = meshgrid(lambdaA,r);
    F = scatteredInterpolant(buf_x,buf_y,buf_z);
    grid_z = F(grid_x, grid_y);
    spectra_interp(:,:,i) = transpose(grid_z);
    if plot_2D_interp
        if i == 1
            figure('Position',[1000 0 500 300])
            [~,h] = contourf(grid_x,grid_y,grid_z);
            h.FaceAlpha = 0.7;
            h.LineStyle = 'none';
            colorbar
            colormap('jet')
            hold on
            T = delaunay(buf_x,buf_y);
            trp = triplot(T,buf_x,buf_y,'w');
            xlabel('Wavelength [nm]')
            ylabel('P [m]')
            xlim([min(lambdaA) max(lambdaA)])
            ylim([min(r) max(r)])
        else
            h.ZData = grid_z;
            T = delaunay(buf_x,buf_y);
            delete(trp)
            trp = triplot(T,buf_x,buf_y,'w');
            drawnow
        end
        title(['Interpolated Spectra at Z = ',num2str(z(i)), '[m]'])
    end
end
%---アーベル変換----
local_spectra = zeros(numel(lambdaA),numel(r),numel(z));
spectra_interp = movmean(spectra_interp,5,1);%二次元スペクトルのスムージング(波長方向)
spectra_interp = movmean(spectra_interp,round((numel(r)-1)/16),2);%二次元スペクトルのスムージング(P方向)
for k=1:numel(z)
    derivative = diff(spectra_interp(:,:,k),1,2)/dr;
    for l=1:numel(lambdaA)
        for i=1:numel(r)
            for j=i:numel(r)-1
                local_spectra(l,i,k) = local_spectra(l,i,k) - 1/pi*derivative(l,j)...
                    *log(r(j+1)*(1+sqrt(1-(r(i)/r(j+1))^2))/(r(j)*(1+sqrt(1-(r(i)/r(j))^2))));%Balandin's Abel inversion
            end
        end
    end
end
w_mov_l = 5;
local_spectra = movmean(local_spectra,w_mov_l,1);%局所スペクトルのスムージング(波長方向)
local_spectra = movmean(local_spectra,round((numel(r)-1)/16),2);%局所スペクトルのスムージング(R方向)
Em_local = zeros(numel(r),numel(z));
Ti_local = zeros(numel(r),numel(z));
Ti_local_max = zeros(numel(r),numel(z));
Ti_local_min = zeros(numel(r),numel(z));
Ti_instru_local = zeros(numel(z),1);
checker = ones(numel(r),numel(z));
RSQ = zeros(numel(r),numel(z));
for i=1:numel(z)
    if i <numel(z)
        Ti_instru_local(i) = sum(Ti_instru_CH(separation(i)+1:separation(i+1)+1)/(separation(i+1)-separation(i)));
    else
        Ti_instru_local(i) = sum(Ti_instru_CH(separation(i)+1:end)/(numel(ch)-separation(i)));
    end
end
%----フィッティング----
col_subp2 = 2;%サブプロット行数
raw_subp2 = 2;%サブプロット列数
n_subp2 = col_subp2*raw_subp2;
col_subp_f = 2;%サブプロット行数
raw_subp_f = 2;%サブプロット列数
n_subp_f = col_subp_f*raw_subp_f;
idx_f_fit = 1;
col_subp_t = 2;%サブプロット行数
raw_subp_t = 2;%サブプロット列数
n_subp_t = col_subp_t*raw_subp_t;
idx_t_fit = 1;
for j=1:numel(z)
    for i=1:numel(r)
        idx_fit = (j-1)*numel(r)+i;%スペクトル番号
        input = local_spectra(:,i,j);
        for l=1:numel(lambdaA)
            if input(l) < 0
                input(l) = -input(l) - min(abs(movmean(input,20)));
            end
        end
        try
            [f,gof] = fit(lambdaA(hw_lambdaA+1-hw_fit:hw_lambdaA+1+hw_fit),input(hw_lambdaA+1-hw_fit:hw_lambdaA+1+hw_fit),'gauss1');
            RSQ(i,j) = gof.rsquare;%決定係数(フィッティング精度指標)
            coeff = coeffvalues(f);%フィッティング係数
            % coeff(3)/resolution/sqrt(2);%IDLのcoeff[2]に等しい(確認用)
            confi = confint(f);%フィッティング係数の95%信頼区間(1行目下限, 2行目上限)
            Ti_local(i,j) = 1.69e8*mass*(2*coeff(3)*sqrt(log(2))/lambda0)^2-Ti_instru_local(j);
            Ti_local_max(i,j) = 1.69e8*mass*(2*confi(2,3)*sqrt(log(2))/lambda0)^2-Ti_instru_local(j);
            Ti_local_min(i,j) = 1.69e8*mass*(2*confi(1,3)*sqrt(log(2))/lambda0)^2-Ti_instru_local(j);
            Em_local(i,j) = resolution*sum(input);
            %外れ値(checker=0)の条件
            checker(i,j) = checker(i,j) * (abs(coeff(2)-lambda0) < 0.1);
            checker(i,j) = checker(i,j) * (coeff(1) > 0);
            checker(i,j) = checker(i,j) * (Em_local(i,j) > 0);
            checker(i,j) = checker(i,j) * (Ti_local(i,j) > 0);
            checker(i,j) = checker(i,j) * (RSQ(i,j) > 0.95);
            switch plot_2D_spectra
                case 'off'
                case 'all'
                    idx_subp2 = mod(idx_fit-1,n_subp2)+1;%サブプロット位置番号
                    if idx_fit == 1
                        figure('Position',[400 500 400 300]);
                        sgtitle('Local Spectra at (Z,R) [mm]')
                    end
                    subplot(col_subp2,raw_subp2,idx_subp2)
                    if idx_fit <= n_subp2
                        if idx_fit == 1
                            fit_x = lambdaA;
                            fit_y = feval(f,fit_x);
                            pfit1 = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                            pfit = repmat(pfit1,[1,1,n_subp2]);
                        else
                            fit_y = feval(f,fit_x);
                            pfit(:,:,idx_subp2) = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                        end
                        xline(lambda0)
                        title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                        xlabel('Wavelength [nm]')
                        ylabel('Strength [a.u.]')
                        xlim([min(lambdaA) max(lambdaA)])
                        ylim([0 inf])
                        legend('off')
                    else
                        fit_y = feval(f,fit_x);
                        pfit(1,:,idx_subp2).YData = fit_y;
                        pfit(2,:,idx_subp2).YData = input;
                        title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                        drawnow
                    end
                case 'good'
                    if checker(i,j) == 1
                        idx_subp_t = mod(idx_t_fit-1,n_subp_t)+1;%サブプロット位置番号
                        if idx_t_fit == 1
                            figure('Position',[400 500 400 300]);
                            sgtitle('Good Spectra at (Z,R) [mm]')
                        end
                        subplot(col_subp_t,raw_subp_t,idx_subp_t)
                        if idx_t_fit <= n_subp_t
                            if idx_t_fit == 1
                                fit_x = lambdaA;
                                fit_y = feval(f,fit_x);
                                pfit_t1 = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                                pfit_t = repmat(pfit_t1,[1,1,n_subp_t]);
                            else
                                fit_y = feval(f,fit_x);
                                pfit_t(:,:,idx_subp_t) = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                            end
                            title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                            xline(lambda0)
                            xlabel('Wavelength [nm]')
                            ylabel('Strength [a.u.]')
                            xlim([min(lambdaA) max(lambdaA)])
                            ylim([0 inf])
                            legend('off')
                        else
                            fit_y = feval(f,fit_x);
                            pfit_t(1,:,idx_subp_t).YData = fit_y;
                            pfit_t(2,:,idx_subp_t).YData = input;
                            title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                            drawnow
                        end
                        idx_t_fit = idx_t_fit + 1;
                    end
                case 'bad'
                    if checker(i,j) == 0
                        idx_subp_f = mod(idx_f_fit-1,n_subp_f)+1;%サブプロット位置番号
                        if idx_f_fit == 1
                            figure('Position',[400 500 400 300])
                            sgtitle('Bad Spectra at (Z,R) [mm]')
                        end
                        subplot(col_subp_f,raw_subp_f,idx_subp_f)
                        if idx_f_fit <= n_subp_f
                            if idx_f_fit == 1
                                fit_x = lambdaA;
                                fit_y = feval(f,fit_x);
                                pfit_f1 = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                                pfit_f = repmat(pfit_f1,[1,1,n_subp_f]);
                            else
                                fit_y = feval(f,fit_x);
                                pfit_f(:,:,idx_subp_f) = plot(fit_x,fit_y,'r-',lambdaA,input,'b+');
                            end
                            title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                            xline(lambda0)
                            xlabel('Wavelength [nm]')
                            ylabel('Strength [a.u.]')
                            xlim([min(lambdaA) max(lambdaA)])
                            ylim([0 inf])
                            legend('off')
                        else
                            fit_y = feval(f,fit_x);
                            pfit_f(1,:,idx_subp_f).YData = fit_y;
                            pfit_f(2,:,idx_subp_f).YData = input;
                            title(sprintf("(%.1f,%.1f)",z(j)*1e3,r(i)*1e3))
                            drawnow
                        end
                        idx_f_fit = idx_f_fit + 1;
                    end
                otherwise
                    warning('Input error in plot_2D_spectra.')%ICCD.lineの入力エラー
                    return;
            end
        catch ME
            checker(i,j) = 0;
            warning(sprintf("Fitting failed in (Z,R) = (%.1f,%.1f)",z(j)*1e3,r(i)*1e3));
        end
    end
end
%NaNを除去
for i=1:numel(r)
    for j=1:numel(z)
        if isnan(Ti_local(i,j))
            Ti_local(i,j) = 0;
            Ti_local_max(i,j) = 0;
            Ti_local_min(i,j) = 0;
            fprintf('Ti_2D(%d,%d) is NaN.',i,j);
        end
    end
end
%外れ値(checker=0)を補間
Ti_local_smooth = filloutliers(Ti_local,"linear");
Ti_local_smooth = movmean(Ti_local_smooth,5,1);
Ti_local_smooth = movmean(Ti_local_smooth,round((numel(r)-1)/16),2);
Ti_local_max_smooth = filloutliers(Ti_local_max,"linear");
Ti_local_max_smooth = movmean(Ti_local_max_smooth,5,1);
Ti_local_max_smooth = movmean(Ti_local_max_smooth,round((numel(r)-1)/16),2);
Ti_local_min_smooth = filloutliers(Ti_local_min,"linear");
Ti_local_min_smooth = movmean(Ti_local_min_smooth,5,1);
Ti_local_min_smooth = movmean(Ti_local_min_smooth,round((numel(r)-1)/16),2);
for i=1:numel(r)
    for j=1:numel(z)
        if checker(i,j) == 0
            Ti_local(i,j) = Ti_local_smooth(i,j);
            Ti_local_max(i,j) = Ti_local_max_smooth(i,j);
            Ti_local_min(i,j) = Ti_local_min_smooth(i,j);
            fprintf('Ti(%d,%d) was interpolated.\n',i,j);
        end
    end
end
Ti_local_offset = mean(Ti_local_min(end,:));
Ti_local = Ti_local - Ti_local_offset;
negative = find(Ti_local<0);
Ti_local(negative) = zeros(size(negative));

z_IDS = z;
r_IDS = r;
Ti_IDS = Ti_local;
Ti_max_IDS = Ti_local_max;
Ti_min_IDS = Ti_local_min;
Em_IDS = Em_local;

%温度分布データを保存
if save_Ti
    if not(exist([pathname.tempdata,'/',num2str(date)],'dir'))
        mkdir(sprintf("%s", pathname.tempdata), sprintf("%s", num2str(date)));
    end
    save([pathname.tempdata,'/',num2str(date),'/shot',num2str(IDS288ch.shot),'_', ...
        num2str(IDS288ch.trg),'us_w=',num2str(IDS288ch.exp_w),'.mat'], ...
        'z_IDS','r_IDS','Ti_IDS','Ti_max_IDS','Ti_min_IDS','Em_IDS')
end
end