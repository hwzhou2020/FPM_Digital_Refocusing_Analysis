function Gini_C = GoG( im )
% edge sparsity-based autofocusing criterion
% 'Gini index of the gradient' (GoG)
% 'Edge sparsity criterion for robust holographic autofocusing'
% 2017_OL_YIBO ZHANG, AYDOGAN OZCAN.

[FX,FY] = gradient(im);
C       = sqrt(abs(FX).^2+abs(FY).^2);

N      = numel(C);
sum_C  = sum(C(:));
a      = sort(C(:));
k      = (1:N)';
Gini_C = 1-2*sum(a.*(N-k+0.5)/(sum_C*N));

end

