function eterna_score = calc_eterna_score_classic( data, pred, BLANK_OUT5, BLANK_OUT3);
% calc_eterna_score_classic( data, pred, BLANK_OUT5, BLANK_OUT3);
%
% Adapted from calc_eterna_score_GUI.m and calc_eterna_score_RHIJU.m in
% eterna github repo's.
%
%  data = [Nres] data, normalized to go from 0 to 1 (~90th percentile) 
%  pred = [Nres] 0 and 1 to marked paired and unpaired.
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues 
%
% (C) R. Das, HHMI, Stanford University, 2023 
%

% In principle could float these and optimize by linear programming, as
% before, but that's pretty confusing actually.
threshold_SHAPE_fixed = 0.5;
min_SHAPE_fixed = 0.0;

goodbins = (1+BLANK_OUT5):(length(data)-BLANK_OUT3);
correct_hit = zeros(1,length(data));
for k = goodbins
    if ( pred(k) == 1 )
        if ( data(k) > (0.25*threshold_SHAPE_fixed + 0.75*min_SHAPE_fixed ) )
            correct_hit(k) = 1;
        end
    else
        if ( data(k) < threshold_SHAPE_fixed)
            correct_hit(k) = 1;
        end
    end
end
eterna_score = sum( correct_hit( goodbins ) )/ length( goodbins ) * 100;

