function [hs,legStr] = make_legend(hs,h,legStr,legName)
if isempty(hs)
    hs = h;
else
    hs = [hs h];
end
if legStr == ""
    legStr = legName;
else
    legStr = [legStr legName];
end
end