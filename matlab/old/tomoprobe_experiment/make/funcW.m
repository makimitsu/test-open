W = zeros(numTheta*numLambda,numVx*numVy);%F->P�ϊ��s��
plot_point = zeros(numTheta*numLambda,2);%�v���b�g�_
for i = 1:numTheta*numLambda
    %W(i,:)�̑Ή��V�t�g�g���ԍ�r���擾�F�Ή��V�t�g�g����ShiftLambda(r,t)
    r = int8(mod(i,numLambda));
    if r == 0
        r = numLambda;
    end
    %W(i,:)�̑Ή������p�x�ԍ�t���擾�F�Ή��p�x��Theta(t)
    t = 1 + idivide(i-1, int8(numLambda), 'floor');
    %���������P�ʃx�N�g��Unit���v�Z
    Unit = [cos(Theta(t)) sin(Theta(t))];
    %���x��ԏ�̊ϑ��ҍ��W���v�Z
    plot_point(i,1) = ShiftLambda(r,t)/lambda0*Vc*Unit(1);
    plot_point(i,2) = ShiftLambda(r,t)/lambda0*Vc*Unit(2);
    for j = 1:numVx*numVy
        %W(:,j)�̑Ή����x�ԍ�x,y���擾�F�Ή����x��Vx(x),Vy(y)
        x = 1 + idivide(j-1, int8(numVy), 'floor');
        y = int8(mod(j,numVy));
        if y == 0
            y = numVy;
        end
        V = [Vx(x), Vy(y)];
        D = dot(V, Unit);%�����������x
        %      W(i,j) = 1 - abs(-ShiftLambda(r)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc));%�ꎟ�֐��t�B���^
        W(i,j) = 1 - ((-ShiftLambda(r,t)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc)))^2;%�񎟊֐��t�B���^
        if W(i,j) < 0
            W(i,j) = 0;
        end
    end
end
