addpath '/Users/shohgookazaki/Documents/GitHub/test-open/Soft X-ray/Four-view';
Newdata = true;

% 再構成条件の定義
newProjectionNumber = 50; %投影数＝視線数の平方根
newGridNumber = 90; %グリッド数（再構成結果の画素数の平方根）

[gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
          s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);

% ファントムテスト用の画像を準備（4視点分）
if Newdata
    data_set(N_projection,gm2d1,'1',1000);
    %data_set(N_projection,gm2d2,'2',1000);
    %data_set(N_projection,gm2d3,'3',1000);
    %data_set(N_projection,gm2d4,'4',1000);
end
