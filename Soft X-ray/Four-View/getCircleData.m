function E_cal = getCircleData(E,n_p)
    E_cal = zeros(n_p);
    k = FindCircle(n_p/2);
    E_cal(k) = E;
    E_cal = double(E_cal);
end