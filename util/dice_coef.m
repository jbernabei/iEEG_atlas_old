function dsc = dice(s1,s2)
    nsect = length(intersect(s1,s2));
    dsc = 2*nsect./(length(s1)+length(s2));
end