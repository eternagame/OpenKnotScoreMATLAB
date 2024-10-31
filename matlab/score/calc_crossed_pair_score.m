function [crossed_pair_score, crossed_pair_quality_score] = calc_crossed_pair_score(data, structure, BLANK_OUT5, BLANK_OUT3, REMOVE_SINGLETS );
% calc_crossed_pair_score( data, structure, BLANK_OUT5, BLANK_OUT3);
%
% crossed_pair_score = 
%  100 x
%   (number of residues in crossed pairs with data < 0.25) /
%     [0.7*(length of region with data - 20)]
%
%  crossed_pair_quality_score =
%  100 x (number of residues in crossed pairs with data < 0.25) /
%          ( number of residues modeled to be crossed pairs in structure )
%
%  data = [Nres] data, normalized to go from 0 to 1 (~90th percentile)
%  structure = [Nres] dot-parens notation for structure, with pseudoknots
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues
%  REMOVE_SINGLETS = remove singlet base pairs (Default 1) -- note time
%                      consuming!
%
% (C) R. Das, HHMI, Stanford University, 2023
%
if ~exist('REMOVE_SINGLETS','var'); REMOVE_SINGLETS = 1; end;

threshold_SHAPE_fixed_pair = 0.25;

crossed_pair_score = 0;
crossed_pair_quality_score = 0;
if all(structure=='x'); return; end;

% get base pairs
bps = convert_structure_to_bps_v3(structure);

% singlets are throwing things off.
if REMOVE_SINGLETS; bps = remove_singlet_bps( bps ); end;

crossed_res = figure_out_which_bps_are_crossed( bps );

bps_filtered = remove_bps_with_BLANK_OUT_regions( bps, length(data), BLANK_OUT5, BLANK_OUT3 );
crossed_res_filtered = figure_out_which_bps_are_crossed( bps_filtered );

num_crossed_pairs  = 0;

total_cross_res = 0;
max_count = 0;
for i = crossed_res
    if i <= BLANK_OUT5; continue; end;
    if i > length(data)-BLANK_OUT3; continue; end;
    max_count = max_count + 1;
    if data( i ) < threshold_SHAPE_fixed_pair
        if ~isempty(find(i == crossed_res_filtered))
            num_crossed_pairs = num_crossed_pairs + 1;
        else % the residue is part of a crossing pair only if base pairs to flanking regions are counted.
            num_crossed_pairs = num_crossed_pairs + 0.5;
        end
    end
end

data_region_length = length(data)- BLANK_OUT5 - BLANK_OUT3;
max_crossed_pairs = 0.7 * max(data_region_length-20, 20);
crossed_pair_score = 100 * min( num_crossed_pairs/max_crossed_pairs, 1.0);

crossed_pair_quality_score = 0;
if max_count > 0; crossed_pair_quality_score = 100 * (num_crossed_pairs/max_count); end;

%fprintf('%s %f\n',structure,crossed_pair_score)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bps_filtered = remove_bps_with_BLANK_OUT_regions( bps, Nres, BLANK_OUT5, BLANK_OUT3);
bps_filtered = [];
Nbp = size(bps,1);
for k = 1:Nbp
    assert( bps(k,1) < bps(k,2) );
    if bps(k,1) <= BLANK_OUT5; continue; end;
    if bps(k,2) <= BLANK_OUT5; continue; end;
    if bps(k,1) > Nres-BLANK_OUT3; continue; end;
    if bps(k,2) > Nres-BLANK_OUT3; continue; end;
    bps_filtered = [bps_filtered; bps(k,:)];
end



