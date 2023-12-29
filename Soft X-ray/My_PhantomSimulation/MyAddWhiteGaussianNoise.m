function I = MyAddWhiteGaussianNoise(I,sn)
    n = max(size(I));
    for i = 1:n
        I(i) = I(i) + I(i)/sn*rand();
    end
    I = awgn(I,sn*5,'measured');
end