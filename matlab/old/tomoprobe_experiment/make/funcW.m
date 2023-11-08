W = zeros(numTheta*numLambda,numVx*numVy);%F->P変換行列
plot_point = zeros(numTheta*numLambda,2);%プロット点
for i = 1:numTheta*numLambda
    %W(i,:)の対応シフト波長番号rを取得：対応シフト波長はShiftLambda(r,t)
    r = int8(mod(i,numLambda));
    if r == 0
        r = numLambda;
    end
    %W(i,:)の対応視線角度番号tを取得：対応角度はTheta(t)
    t = 1 + idivide(i-1, int8(numLambda), 'floor');
    %視線方向単位ベクトルUnitを計算
    Unit = [cos(Theta(t)) sin(Theta(t))];
    %速度空間上の観測者座標を計算
    plot_point(i,1) = ShiftLambda(r,t)/lambda0*Vc*Unit(1);
    plot_point(i,2) = ShiftLambda(r,t)/lambda0*Vc*Unit(2);
    for j = 1:numVx*numVy
        %W(:,j)の対応速度番号x,yを取得：対応速度はVx(x),Vy(y)
        x = 1 + idivide(j-1, int8(numVy), 'floor');
        y = int8(mod(j,numVy));
        if y == 0
            y = numVy;
        end
        V = [Vx(x), Vy(y)];
        D = dot(V, Unit);%視線方向速度
        %      W(i,j) = 1 - abs(-ShiftLambda(r)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc));%一次関数フィルタ
        W(i,j) = 1 - ((-ShiftLambda(r,t)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc)))^2;%二次関数フィルタ
        if W(i,j) < 0
            W(i,j) = 0;
        end
    end
end
