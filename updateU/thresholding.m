function V = thresholding(V,thres)
%
Temp = max(abs(V)-thres,0);
V = sign(V).*Temp;