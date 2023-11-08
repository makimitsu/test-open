max_i = 100;
p = zeros(2,max_i);
curv = zeros(2,max_i);%1s–ÚƒKƒ“ƒ}A2s–Ú‹È—¦
if draw_analisis
    figure('Position',[300 150 550 500])
end
for i = 1:max_i
    ganma = 1.1^(-i+20);
    curv(1,i) = ganma;
    ReF = inv(W.'*W + ganma * eye(numVx*numVy))*W.'*P;
    p(1,i) = log10(norm(W * ReF - P));
    p(2,i) = log10(norm(ReF));
    if draw_analisis
        plot(p(1,i),p(2,i), 'r+')
        hold on
    end
    if i>1 && i<max_i
        curv(2,i) = 2*det([p(:,i+1)-p(:,i), p(:,i-1)-p(:,i)])/...
            (norm(p(:,i+1)-p(:,i))*norm(p(:,i-1)-p(:,i))*norm(p(:,i+1)-p(:,i-1)));
    end
end
if draw_analisis
    title('L-curve','FontSize',20);
    xlabel('log10(WF-P)','FontSize',15);
    ylabel('log10(F)','FontSize',15);
    hold off
    figure('Position',[500 150 550 500])
    plot(p(1,:),curv(2,:))
end
[M,I] = max(curv.');
ganma = curv(1,I(2));%Å“K‚ÈƒKƒ“ƒ}‚ð—^‚¦‚é
ReF = inv(W.'*W + ganma * eye(numVx*numVy))*W.'*P;
