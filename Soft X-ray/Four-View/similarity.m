%Lf=gの類似度計算

newProjectionNumber = 50; %投影数＝視線数の平方根
newGridNumber = 90; %グリッド数（再構成結果の画素数の平方根）

[gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
          s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);
matrixPath = '/Users/shohgookazaki/Library/CloudStorage/GoogleDrive-shohgo-okazaki@g.ecc.u-tokyo.ac.jp/My Drive/OnoLab/data/result_matrix/cGAN/240621/shot4/3.mat';
sxrPath = '/Users/shohgookazaki/Library/CloudStorage/GoogleDrive-shohgo-okazaki@g.ecc.u-tokyo.ac.jp/My Drive/OnoLab/data/SXR_data/240621/shot4/3.mat';

load(matrixPath,'EE1','EE2','EE3','EE4');
load(sxrPath, 'sxr1', 'sxr2','sxr3','sxr4');
sim = simi(gm2d1,N_projection,EE1,sxr1);

matrixPath = '/Users/shohgookazaki/Library/CloudStorage/GoogleDrive-shohgo-okazaki@g.ecc.u-tokyo.ac.jp/My Drive/OnoLab/data/result_matrix/LF_LR/240621/shot4/3.mat';
load(matrixPath,'EE1','EE2','EE3','EE4');
sim2 = simi(gm2d1,N_projection,EE1,sxr1);

fprintf('cGANsim = %d\n', sim);
fprintf('LF_NLRsim = %d\n', sim2);

function sim = simi(gm2d,N_projection,EE, sxr)
    EE = EE(:);
    
    sxr_cal = gm2d*EE;
    
    n_p = N_projection;
    sxr_calcal = zeros(n_p);
    k=FindCircle(n_p/2);
    sxr_calcal(k) = sxr_cal;
    
    [ssimval, ssimmap] = ssim(sxr, sxr_calcal);
        
    sim = ssimval;
end