function bps_filtered = filter_bps(bps,goodpos)
% filter_bps(bps,goodpos)
%
% Inputs:
%   bps = [Nbpsx2] list of base pairs
%   goodpos = positions to keep
%
% Output
%   bps_filtered = [Nbps_filteredx2] list of base pairs where both
%      positions are in goodpos
%
% (C) R. Das, Stanford/HHMI 2024
bps_filtered = [];
for i = 1:size(bps,1)
    if ~any(bps(i,1)==goodpos); continue; end
    if ~any(bps(i,2)==goodpos); continue; end
    bps_filtered = [bps_filtered; bps(i,:)];
end