%Lf=gの類似度計算

date = '240621';
results = [];

newProjectionNumber = 50; %投影数＝視線数の平方根
newGridNumber = 90; %グリッド数（再構成結果の画素数の平方根）
if evalin('base', 'exist(''N_projection'', ''var'')')
    NP = evalin('base', 'N_projection');
    if NP ~= newProjectionNumber
        [gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
                  s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);
    end
else
    [gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
              s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);
end

matpath = getenv('SXR_MATRIX_DIR');
sxrpath = getenv('SXR_DATA_DIR');


for shotnum = 5%[1:2, 4:22, 24:55]
    for matnum = 5
        shot = num2str(shotnum);
        mat = num2str(matnum);

        matrixPath_cGAN = strcat(matpath, '/cGAN/', date, '/shot', shot, '/', mat, '.mat');
        matrixPath_LF_LR = strcat(matpath, '/LF_LR/', date, '/shot', shot, '/', mat, '.mat');

        sxrPath = strcat(sxrpath,'/', date, '/shot', shot, '/', mat, '.mat');
        load(sxrPath, 'sxr1');

        load(matrixPath_cGAN,'EE1');
        cGAN_sim = simi(gm2d1,N_projection,EE1,sxr1);
        
        cGAN_sxr = gm2d1*EE1(:);

        n_p = N_projection;
        cGAN_sxr_cal = zeros(n_p);
        k=FindCircle(n_p/2);
        cGAN_sxr_cal(k) = cGAN_sxr;
        cGAN_sxr_cal = cGAN_sxr_cal(:);

        load(matrixPath_LF_LR,'EE1');
        LF_LRsim = simi(gm2d1,N_projection,EE1,sxr1);

        LF_LR_sxr = gm2d1*EE1(:);
        LF_LR_sxr_cal = zeros(n_p);
        LF_LR_sxr_cal(k) = LF_LR_sxr;
        LF_LR_sxr_cal = LF_LR_sxr_cal(:);
        
        x = 1:2500;

        plot(x, cGAN_sxr_cal,'g', x, LF_LR_sxr_cal, 'b--o');
        
        
        results = [results; {shot, mat, cGAN_sim, LF_LRsim}];
    end
end
% Convert results cell array to table
resultsTable = cell2table(results, 'VariableNames', {'Shot', 'Mat','cGAN_normsim', 'LF_LRsim'});

% Write table to CSV file
writetable(resultsTable, 'SSIM_results.csv');

function sim = simi(gm2d,N_projection,EE, sxr)

    EE = EE(:);
    
    sxr_cal = gm2d*EE;
    
    sxr_calcal = getCircleData(sxr_cal,N_projection);
    
    sxr = double(sxr);

    [ssimval, ssimimage] = ssim(sxr, sxr_calcal);
        
    sim = ssimval;
end