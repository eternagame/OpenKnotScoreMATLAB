function [F1,prec,rec] = F1_bps(bps_pred, bps_ref, crossed_pair)
% [F1,prec,rec] = F1_bps(bps_pred, bps_ref, crossed_pair)
%
% F1 score: harmonic mean of precision and recall.
%
% Conventions for crossed_pair score & when to return nan vs. 0
%   developed by Shujun He
% 
% Inputs:
%  bps_pred = Nx2 list of base pairs (i<j), predicted
%  bps_ref  = Nx2 list of base pairs (i<j), reference
%  crossed_pair = calculate F1 over crossed pairs. (default 0)
%
% (C) R. Das, Stanford University & HHMI, 2024

if ~exist('crossed_pair','var'); crossed_pair = 0; end;

if crossed_pair
    bps_pred = detect_crossed_pairs(bps_pred);
    bps_ref = detect_crossed_pairs(bps_ref);
end

if length(bps_ref) == 0 & length(bps_pred)>0; F1 = 0; prec = 0; rec = 0; return; end;
if length(bps_ref) == 0 & length(bps_pred)==0; F1 = nan; prec = nan; rec = nan; return; end;

prec = precision_bps(bps_pred,bps_ref);
rec  = recall_bps(bps_pred,bps_ref);
if (prec+rec == 0); F1 = 0; return; end;

F1 = 2*prec*rec/(prec + rec);


