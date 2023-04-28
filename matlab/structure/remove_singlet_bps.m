function bps = remove_singlet_bps(bps);
% bps = remove_singlet_bps(bps);
%
% bps = [Nbp x 2] indicses into base pairs
%
% (C) R. Das, Stanford & HHMI, 2023.

stems = parse_stems_from_bps( bps );

stems_new = stems;
for i = 1:length(stems);
    stem = stems{i};
    if size(stem,1) == 1;  % singlet!
        idx = find( bps(:,1) == stem(1) & find(bps(:,2)==stem(2)));
        bps = bps( [1:(idx-1) (idx+1):end], :); % take out of bp list
    end
end

