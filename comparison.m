%% Script for comparing TH*, TH-Δ, TH-FRAC
disp(['e    ϑ     ϑ-Δ     ϑf    ϑ*   α'])
for i = 1:10
    n =  20;     A = graph(n);
%     thstar = star_SDP(A);
%     thfrac = frac_SDP(A);
    thddel = ddel_SDP(A);
    disp(num2str([nnz(A)/2 theta(A) thddel alpha(A)]))
end