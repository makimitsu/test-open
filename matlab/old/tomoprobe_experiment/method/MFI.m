%‰Šú‰ð‚ð¶¬
run Tikhonov1.m
run goto0.m
j = 1;
res = 1;
if draw_analisis
    figure('Position',[300 150 550 500])
end

while res > 0.1 && j<50
    if draw_analisis
        plot(j,res, 'r+')
        hold on
    end
    buf = ReF;
    M = zeros(numVx*numVy);
    for k = 1:numVx*numVy
        if ReF(k) > 0
            M(k,k)=1/ReF(k);
        end
    end
    max_M = max(max(M));
    for k = 1:numVx*numVy
        if ReF(k) <= 0
            M(k,k)=max_M;
        end
    end
    penalty = gradVx.'*M*gradVx + gradVy.'*M*gradVy;
    ReF = inv(W.'*W + ganma * penalty)*W.'*P;
    res = norm(ReF - buf);
    j = j+1;
end
if draw_analisis
    xlabel('Iteration count','FontSize',15);
    ylabel('Residual','FontSize',15);
    hold off
end
