% データ読み込み
load mri
D = squeeze(D);

% Plot
figure(1)
[~,handle_surf] = contour(D(:,:,1));
title('スカラー ボリューム データの可視化');

% スライダーバー作成
handle_slider = uicontrol('Style','slider','Position',[10 50 20 340],'Min',1,'Max',27,'Value',1);

% イベントリスナーの作成
addlistener(handle_slider,'Value','PostSet',@(event,obj) update(event,obj,D,handle_surf,handle_slider));

% 描画のアップデート処理用関数

function update(event,obj,D,handle_surf,handle_slider)
index = round(handle_slider.Value);
handle_surf.ZData = double(D(:,:,index));
end