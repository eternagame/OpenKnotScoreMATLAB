function bps_out = detect_crossed_pairs(bps);
% bps_out = detect_crossed_pairs(bps);
%
% Just get the pseudoknots (pairs that cross some other pairs)
%
% Input
%  bps = [N x 2] list of base pairs (with i<j)
%
% Output
%  bps = [N x 2] list of just crossed base pairs (with i<j)
%
% (C) R. Das, Stanford University, 2024

bps_out = [];
crossed_res = figure_out_which_bps_are_crossed( bps );
for i = 1:size(bps,1)
    if any(bps(i,1) == crossed_res)
        bps_out = [bps_out; bps(i,1), bps(i,2)];
    end
end
