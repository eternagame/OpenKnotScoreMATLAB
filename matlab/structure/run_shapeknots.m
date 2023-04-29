function structures_shapeknots = run_shapeknots( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3 );
% structures_shapeknots = run_shapeknots( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3 );
%
% Requires Biers to be installed.
%
% (C) R. Das, Stanford/HHMI 2023

structures_shapeknots = run_rnastructure( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3, 1)