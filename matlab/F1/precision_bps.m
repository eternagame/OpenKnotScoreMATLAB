function precision = precision_bps( bps_pred, bps_ref )
%  precision = precision_bps( bps_pred, bps_ref )
%
% Fraction of predicted base pairs found in reference.
%
% Inputs:
%  bps_pred = Nx2 list of base pairs (i<j), predicted
%  bps_ref  = Nx2 list of base pairs (i<j), reference
%
% (C) R. Das, Stanford University & HHMI, 2024

count = 0;
N = size( bps_pred,1 );
if size(bps_pred,1) == 0 & size(bps_ref,1) == 0;; precision = nan; return; end;
if size(bps_pred,1) == 0; precision = 0; return; end;
if size(bps_ref,1) == 0; precision = 0; return; end;
for n = 1:N
    i = bps_pred(n,1);
    j = bps_pred(n,2);
    idx = find(bps_ref(:,1) == i ); 
    if length(idx) > 0 & bps_ref(idx,2) == j, count = count + 1; end;
end
precision = count/N;

