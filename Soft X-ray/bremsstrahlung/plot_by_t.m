data = readmatrix('xpoint_data_240621_3kV.xlsx');
TF = 3;
transmission_n = 10;

% 各列の名前（適宜変更）
shot = data(:, 1);
time = data(:, 2);
energy = data(:, 3:transmission_n+2); % 10個のエネルギー
intensity = data(:, transmission_n+3:2*transmission_n+2); % 10個のintensity

% 時間の範囲設定
time_ranges = [460 470; 470 480;480 490;490 500];
%time_ranges = [450 470; 470 490; 490 510];
%time_ranges = [450 500];
n_ranges = size(time_ranges, 1);

% 平均・標準誤差の計算結果を保存する配列
mean_energy = zeros(n_ranges, transmission_n);
mean_intensity = zeros(n_ranges, transmission_n);
stderr_intensity = zeros(n_ranges, transmission_n);

% 外れ値のしきい値設定(この値以下の行は考慮しない)
threshold = 1e-4;
% 線形性のしきい値（R^2が0.99以上なら直線とみなす）
correlation_threshold = 0.50;

for i = 1:n_ranges
    % 時間範囲に該当する行を抽出
    range_indices = time >= time_ranges(i, 1) & time < time_ranges(i, 2);
    energy_subset = energy(range_indices, :);
    intensity_subset = intensity(range_indices, :);
    
    % 各エネルギーごとに処理
    for j = 1:10
        % 直線性のチェック（相関係数）
        corr_coef = corr(energy_subset(:, j), intensity_subset(:, j));
        if abs(corr_coef) < correlation_threshold
            % 条件を満たすデータのみで平均と標準誤差を計算
            valid_indices = intensity_subset(:, j) > threshold;
            filtered_intensity = intensity_subset(valid_indices, j);
            filtered_energy = energy_subset(valid_indices, j);
            
            mean_energy(i, j) = mean(filtered_energy);
            mean_intensity(i, j) = mean(filtered_intensity);
            stderr_intensity(i, j) = std(filtered_intensity) / sqrt(sum(valid_indices));
        else
            % 直線性が高い場合はNaNとしてスキップ
            mean_energy(i, j) = NaN;
            mean_intensity(i, j) = NaN;
            stderr_intensity(i, j) = NaN;
        end
    end
end
% 
% プロット
figure;
hold on;
for i = 1:n_ranges
    valid_points = ~isnan(mean_energy(i, :));
    errorbar(mean_energy(i, valid_points), mean_intensity(i, valid_points), stderr_intensity(i, valid_points), '-o', 'DisplayName', sprintf('%d-%d', time_ranges(i, 1), time_ranges(i, 2)));
end

x_min = 10;
x_max = 300;


%------------------------モデルプロット-------------------------
me = 9.11e-31; %kg
e = 1.6*10^(-19); %C
e0 = 8.85e-12; %F/m
mi = 6.64e-26; %kg
mr = me*mi/(me+mi);
TA_eV = 250; %[eV]
TB_eV = 1000;
% TA_K = TA_eV * 1.16045e4; %[K]
ne = 1*1e19; %m^-3
kB = 1.38e-23; %[J/K]
kBTA = TA_eV*e; %eV→J（kBはJ/K）
kBTB = TB_eV*e;
% lambdaD = sqrt(e0*kBT/(ne*e^2)); %m
% Lambda = 4*pi()*lambdaD^3*ne/3; %無次元数
Z = 1; %暫定的な価数
h = 6.63e-34; %プランク定数
% Ry = me*e^4/(2*h^2);
Ry = 13.6*e;
c = 3e8;
s = 2;
% vA_th = sqrt(8*kB*TA_K/me/pi()); %[m/s]
% vB_th = sqrt(8*kB*TB_K/me/pi());
% 
% EA_th_J = 1/2*me*vA_th^2; %[J]
% EB_th_J = 1/2*me*vB_th^2; %[J]

% e_brem = 10;
% f_brem = e_brem*e/h;
% w_brem = f_brem*2*pi();
% w_brem = 1e13:1e13:1e15;
% w_brem = 10.^(linspace(10,17,100));

f_brem = 10.^(linspace(10,17,100));
ve = 1e5:1e5:1e8;

dWA = log(me*ve.'.^2./(h*2*pi()*f_brem))./ve.';
dWB = log(me*ve.'.^2./(h*2*pi()*f_brem))./ve.';
dWB(ve < 1e6) = 0;  % Set dW to 0 where ve is smaller than 1e6

A = ve.^2.*exp(-me*ve.^2/(2*kBTA));
B = ve.^(2-2*s);
B = ve.^2.*exp(-me*ve.^2/(2*kBTB));

num_eA = sum(A);
num_eB = sum(B);

emissivityA = sum(dWA.*A.')/num_eA;
emissivityA(emissivityA<0) = 0;
emissivityB = sum(dWB.*B.')/num_eB;

energy = h*f_brem/e;

% emissivityB = energy.^(1-s)*10^(-5.5);

threshold = 88;
x = 1:1:100;
emissivitySA = emissivityA.*(1+exp((1*(x-threshold)))).^(-1);
emissivitySB = emissivityB.*(1+exp(-(1*(x-threshold)))).^(-1);
emissivityS = emissivitySA+emissivitySB;
plot(energy,emissivityS*1e8, 'DisplayName', sprintf('model'));
%---------モデルプロットした-------------

xlabel('Energy');
ylabel('Intensity');
title(strcat('Average intensity when TF is ', num2str(TF), 'kV')); % タイトル
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim([10 300]);
%ylim([1e-1 10]);
legend;
hold off;