function ret = maxabs(array)

max_val = max(array);
min_val = min(array);
if abs(max_val) >= abs(min_val)
    ret = max_val;
else
    ret = min_val;
end

end