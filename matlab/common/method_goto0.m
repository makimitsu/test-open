function [F] = method_goto0(F)
for i = 1:numel(F)
    if F(i) < 0
        F(i) = 0;
    end
end
end