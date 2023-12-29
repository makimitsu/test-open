function plot_sxr_multi(doSave,doFilter,NL,M,K,gm2d,U,s,v,EEorigin,range)

    % 生画像の取得
    [rawImage] = generateFiberImage(gm2d,flipud(EEorigin));

    % 非線形フィルターをかける（必要があれば）
    if doFilter
        [rawImage,~] = imnlmfilt(rawImage,'SearchWindowSize',11,'ComparisonWindowSize',5);
        imagesc(rawImage);
    end
    
    % ベクトル形式の画像データの読み込み
    VectorImage = rawImage';

    % 再構成計算
    EE = get_distribution(M,K,gm2d,U,s,v,VectorImage,NL);

    plot_save_sxr(range,EE,EEorigin);
end