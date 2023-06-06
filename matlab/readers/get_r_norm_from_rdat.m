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

sequences = r.sequences;
BLANK_OUT5 = r.seqpos(1)-1;
BLANK_OUT3 = length(r.sequence) - r.seqpos(end);
r_norm = padarray( padarray(r.reactivity,BLANK_OUT5,NaN,'pre'), BLANK_OUT3,NaN,'post')';
r_norm_err = padarray( padarray(r.reactivity_error,BLANK_OUT5,NaN,'pre'), BLANK_OUT3,NaN,'post')';