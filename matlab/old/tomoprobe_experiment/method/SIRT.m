rho = 10^(-3);%ŠÉ˜aƒpƒ‰ƒ[ƒ^
j = 1; 
res = 1;
ReF = zeros(numVx*numVy,1);
if draw_analisis
    figure('Position',[300 150 550 500])
end
while res > 1*10^-2 && j<100
    buf = ReF;
    ReF = ReF + rho*W.'*(P-W*ReF);
    run goto0.m
    res = norm(ReF - buf);
    if draw_analisis
        plot(j,res, 'r+')
        hold on
    end
    j = j+1;
end
if draw_analisis
    xlabel('Iteration count','FontSize',15);
    ylabel('Residual','FontSize',15);
    hold off
end
