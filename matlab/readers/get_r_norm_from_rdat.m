function [r_norm,sequences,BLANK_OUT5,BLANK_OUT3,r_norm_err] = get_r_norm_from_rdat( r );
%
% Inputs
%  r = RDAT object with data for Ndesigns designs
%
% Outputs
%  r_norm    = [Ndesigns x Nres] data values. (Padded with NaN's for
%                missing residues at 5' and 3' ends.)
%  sequences = {cell of Nres strings} sequences 
%  BLANK_OUT5 = number of residues from 5' for which there is no data. 
%  BLANK_OUT3 = number of residues at 3' for which there is no data. 
%  r_norm_err = [Ndesigns x Nres] data error values. (Padded with NaN's for
%                missing residues at 5' and 3' ends.)
% 
% (C) R. Das, Stanford University


% need to make sure all sequences are same length! 
sequences = r.sequences;
if r.seqpos(2)<r.seqpos(1);
    r.seqpos = r.seqpos(end:-1:1);
    r.reactivity = r.reactivity(end:-1:1,:);
    r.reactivity_error = r.reactivity_error(end:-1:1,:);
end
BLANK_OUT5 = r.seqpos(1)-r.offset - 1;
BLANK_OUT3 = length(r.sequence) - (r.seqpos(end)-r.offset);
if BLANK_OUT5 < 0 | BLANK_OUT3 < 0
    warning(sprintf('BLANK_OUT5 = %d, BLANK_OUT3 = %d. not outputting anything',BLANK_OUT5,BLANK_OUT3));
    r_norm = []; r_norm_err = [];
    return;
end

% Eterna library with mixed lengths -- needs some special padding...
seqlen = cellfun(@length,r.sequences);
if length(unique(seqlen)) > 1 
    if max(seqlen) == length(sequences{end}) & length(r.seqpos)==size(r.reactivity,1) & all(r.reactivity(:,end)==0) 
        BLANK_OUT3_PREV = BLANK_OUT3;
        BLANK_OUT3 = max(seqlen(1:end-1)) - (r.seqpos(end)-r.offset);
        warning(sprintf('Only last sequence is long. resetting BLANK_OUT3 from %d to %d',BLANK_OUT3_PREV,BLANK_OUT3));
    end
    assert(length(r.sequence) == max(seqlen));
    assert(length(r.seqpos) == size(r.reactivity,1));
    for i = 1:size(r.reactivity,2)
        %badpos = [ (length(r.seqpos) + length(r.sequences{i})- (r.seqpos(end)-r.offset)) : size(r.reactivity,1)];
        badpos = [ length(r.seqpos) - BLANK_OUT3  -  (size(r.reactivity,1) - length(r.sequences{i})): size(r.reactivity,1)];
        r.reactivity(badpos,i) = NaN;
        r.reactivity_error(badpos,i) = NaN;
    end
end

r_norm = padarray( padarray(r.reactivity,BLANK_OUT5,NaN,'pre'), BLANK_OUT3,NaN,'post')';
r_norm_err = padarray( padarray(r.reactivity_error,BLANK_OUT5,NaN,'pre'), BLANK_OUT3,NaN,'post')';

