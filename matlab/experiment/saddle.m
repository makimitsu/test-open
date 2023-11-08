for j = 1:size(psi,2)
    if j == 1
        [saddle_psi,saddle_r] = max(psi(:,j));
        saddle_z = j;
    else
        if max(psi(:,j)) < saddle_psi
            [saddle_psi,saddle_r] = max(psi(:,j));
            saddle_z = j;
        end
    end
end