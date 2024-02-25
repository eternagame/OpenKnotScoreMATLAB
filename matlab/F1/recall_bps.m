function recall = recall_bps( bps_pred, bps_ref)
%  recall = recall_bps( bps_pred, bps_ref )
%
% Fraction of reference base pairs found in prediction.
%
% Inputs:
%  bps_pred = Nx2 list of base pairs (i<j), predicted
%  bps_ref  = Nx2 list of base pairs (i<j), reference
%
% (C) R. Das, Stanford University & HHMI, 2024

recall = precision_bps(bps_ref,bps_pred);

