function xi = cholesky(A, b) % 変形コレスキー法
    [~, check] = chol(A); %正定値対称行列チェック
    if check ~= 0
        error('no cholesky');
    end
    [An, ~] = size(A);
    D = zeros(An);
    U = eye(An);
    D(1,1) = A(1,1);
    for j=2:An
        U(1,j) = A(1,j)/D(1,1);
    end
    for i=2:An
        sig = zeros(i-1,1);
        for m = 1:i-1
            sig(m) = U(m,i)*U(m,i)*D(m,m);
        end
        sig = sum(sig);
        D(i,i) = A(i,i) - sig;
        if i < An
            for j = i+1:An
                sig = zeros(i-1,1);
                for m=1:i-1
                    sig(m) = U(m,j)*D(m,m)*U(m,i);
                end
                sig = sum(sig);
                U(i,j) = (A(i,j)-sig)/D(i,i);
            end
        end
    end

    z = zeros(An,1);
    z(1,1) = b(1,1);
    for i=2:An
        sig = zeros(i-1,1);
        for m=1:i-1
            sig(m) = U(m,i)*z(m,1);
        end
        sig = sum(sig);
        z(i,1) = b(i,1)-sig;
    end

    xi = zeros(An,1);
    xi(An,1) = z(An,1)/D(An,An);
    for h=1:An-1
        i = An-h;
        sig = zeros(An-i,1);
        for m=i+1:An
            sig(m) = U(i,m)*xi(m,1);
        end
        sig = sum(sig);
        xi(i,1) = z(i,1)/D(i,i) - sig;
    end
end