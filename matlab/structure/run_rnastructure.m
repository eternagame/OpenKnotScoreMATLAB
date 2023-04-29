function structures_rnastructure = run_rnastructure( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3, shapeknots );
% structures_rnastructure = run_rnastructure( sequences, r_norm(:,:,2), good_idx, BLANK_OUT5, BLANK_OUT3, shapeknots );
%
% Requires Biers to be installed.
%
% Inputs
%  sequences = list of Nseq sequences
%  r_norm    = [Nseq x Nres] reactivity for all sequences
%  good_idx  = which indices in sequences to focus on
%  BLANK_OUT5 = ignore this many 5' end residues
%  BLANK_OUT3 = ignore this many 3' end residues
%  shapeknots = try to infer pseudoknots (default: 0).
%
% Outputs:
%  structures_rnastructure = cell of structures with length
%                            matching number of indices ingood_idx  
%
% (C) R. Das, Stanford/HHMI 2023


if ~exist( 'shapeknots','var'); shapeknots = 0; end;
if isempty(good_idx); good_idx = [1:size(r_norm,1)]; end;

structures_rnastructure = cell(length(good_idx),1);
parfor i = 1:length(good_idx)
    idx = good_idx(i);
    sequence = sequences{idx};
    N = length(sequence);
    which_pos = (BLANK_OUT5+1:N-BLANK_OUT3);
    area_shape = r_norm(idx,which_pos)/2;
    offset = 0;
    seqpos = which_pos;
    [structure, bpp, SHAPE_out ] = rna_structure( sequence, area_shape, offset, seqpos, [],0,shapeknots);
    structures_rnastructure{i} = structure;
end

