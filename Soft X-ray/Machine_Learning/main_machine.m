addpath '/Users/shohgookazaki/Documents/GitHub/test-open/Soft X-ray/Four-view';
Newdata = true;

% 再構成条件の定義
newProjectionNumber = 50; %投影数＝視線数の平方根
newGridNumber = 90; %グリッド数（再構成結果の画素数の平方根）

[gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
          s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);

datasetnum = 5000;

% ファントムテスト用の画像を準備（4視点分）
if Newdata
    data_set(N_projection,M,K,gm2d1,U1,s1,v1,'1',datasetnum);
    %data_set(N_projection,M,K,gm2d2,U2,s2,v2,'2',datasetnum);
    %data_set(N_projection,M,K,gm2d3,U3,s3,v3,,'3',datasetnum);
    %data_set(N_projection,M,K,gm2d4,U4,s4,v4,,'4',datasetnum);
    


    
end
