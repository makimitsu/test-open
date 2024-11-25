function [] = transmission(n)
% Excelファイルのパスを指定
filename = 'transmission.xlsx';

% シート名リスト
sheets = {'al10','al25','mylar10', 'mylar25'};

% パラメータ設定
e_min = 10;  % 下限エネルギー (eV)
e_max = 300; % 上限エネルギー (eV)

% エネルギー範囲の分割境界を生成
edges = logspace(log10(e_min), log10(e_max), n+1);
% 各区間の平均を計算
means = zeros(1, n);
for i = 1:n
    means(i) = (edges(i) + edges(i + 1)) / 2;
end

% 各シートについて処理を行う
for k = 1:length(sheets)
    sheet = sheets{k};
    
    % データの読み込み
    data = readtable(filename, 'Sheet', sheet, 'VariableNamingRule', 'preserve');
    photon_energy = data{:, 1}; % 1列目: エネルギー (eV)
    transmission = data{:, 2};  % 2列目: 透過率
    
    % nx1の行列を初期化
    transmission_avg = zeros(n, 1);
    
    % 各エネルギー範囲ごとに透過率の平均を計算
    for i = 1:n
        % 現在のエネルギー範囲のインデックスを取得
        idx = photon_energy >= edges(i) & photon_energy < edges(i+1);
        
        % 該当する透過率があれば平均を計算し、なければ0を代入
        if any(idx)
            transmission_avg(i) = mean(transmission(idx));
        else
            transmission_avg(i) = 0; % データがない場合は0を代入
        end
    end
    
    % 結果を保存・表示
    % fprintf('Sheet: %s\n', sheet);
    % disp(transmission_avg);
    
    % 必要に応じて個別の変数名で保存
    eval([sheet, '_avg = transmission_avg;']); % mylar1.0u_avg, mylar2.5u_avgなどとして保存
end

transmissionfile = sprintf('transmission_data%d.mat',n);
transmission_matrix = [
    al10_avg';   % 1行目: al1.0u_avg
    al25_avg';   % 2行目: al2.5u_avg
    mylar10_avg'; % 3行目: mylar1.0u_avg
    mylar25_avg'  % 4行目: mylar2.5u_avg
];

%for Tikhonov
C = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1) - 2*diag(ones(n,1));
[U,S,V] = svd(transmission_matrix*(C^(-1)),'econ');
v = (C^(-1)*V);
[M,K] = size(transmission_matrix);
if K>M
    v = v(:,1:M);
end
s = (diag(S)).';
if M>K
    s = [s zeros(1,M-K)];
end

% 行列を.matファイルに保存
save(transmissionfile, 'transmission_matrix','means', 'U','s','v','M','K');

end
