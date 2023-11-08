max_i = 150;
p = zeros(2,max_i);
curv = zeros(2,max_i);%1�s�ڃK���}�A2�s�ڋȗ�
gradgradVx = zeros(numVx*numVy);
gradgradVy = zeros(numVx*numVy);
for i = 1:numVx*numVy
    x1 = 1 + idivide(i-1, int8(numVy), 'floor');
    y1 = int8(mod(i,numVy));
    if x1 ~= 1 && x1 ~= numVx
        for j = 1:numVx*numVy
            x2 = 1 + idivide(j-1, int8(numVy), 'floor');
            y2 = int8(mod(j,numVy));
            if y1 == y2
                if x2 == x1 - 1
                    gradgradVx(i,j) = 1/deltaVx^2;
                elseif x2 == x1
                    gradgradVx(i,j) = -2/deltaVx^2;
                elseif x2 == x1 + 1
                    gradgradVx(i,j) = 1/deltaVx^2;
                end
            end
        end
    end
    if y1 ~= 1 && y1 ~= numVy
        for j = 1:numVx*numVy
            x2 = 1 + idivide(j-1, int8(numVy), 'floor');
            y2 = int8(mod(j,numVy));
            if x1 == x2
                if y2 == y1 - 1
                    gradgradVy(i,j) = 1/deltaVy^2;
                elseif y2 == y1
                    gradgradVy(i,j) = -2/deltaVy^2;
                elseif y2 == y1 + 1
                    gradgradVy(i,j) = 1/deltaVy^2;
                end
            end
        end
    end
end
gradgrad = gradgradVx.'*gradgradVx + gradgradVy.'*gradgradVy;
if draw_analisis
    figure('Position',[300 150 550 500])
end
for i = 1:max_i
    ganma = 1.1^(-i+10);
    curv(1,i) = ganma;
    ReF = inv(W.'*W + ganma * gradgrad)*W.'*P;
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
    ylabel('log10(��^2F)','FontSize',15);
    hold off
    figure('Position',[500 150 550 500])
    plot(p(1,:),curv(2,:))
end
[M,I] = max(curv.');
ganma = curv(1,I(2));%�œK�ȃK���}��^����
ReF = inv(W.'*W + ganma * gradgrad)*W.'*P;
