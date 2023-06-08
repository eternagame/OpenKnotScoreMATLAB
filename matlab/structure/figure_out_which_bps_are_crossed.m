function crossed_res = figure_out_which_bps_are_crossed( bps );
% crossed_res = figure_out_which_bps_are_crossed( bps );
% figure out which pairs are crossed pairs.
% (C) R. Das, Stanford/HHMI 2023.

crossed_res = [];
Nbp = size(bps,1);
for k = 1:Nbp
    assert( bps(k,1) < bps(k,2) );
    for j = (k+1):Nbp        
        if ( bps(k,1) < bps(j,1) & bps(j,1) < bps(k,2) & bps(k,2) < bps(j,2) ) | ...
                ( bps(j,1) < bps(k,1) & bps(k,1) < bps(j,2) & bps(j,2) < bps(k,2) )
            crossed_res = [crossed_res, bps(k,1)];
            crossed_res = [crossed_res, bps(k,2)];
            crossed_res = [crossed_res, bps(j,1)];
            crossed_res = [crossed_res, bps(j,2)];
        end
    end
end
crossed_res = unique(crossed_res);
