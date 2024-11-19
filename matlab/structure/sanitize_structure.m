function structure = sanitize_structure( structure, REMOVE_SINGLETS )
% structure = sanitize_structure( structure, REMOVE_SINGLETS )
%
% Input
%  structure = structure in dot parens/brackets notation
%  REMOVE_SINGLETS = remove singlet base pairs (Default 0)
%
% Output
%  structure = structure in dot parens/brackets notation
%
% Requires Biers to be installed.
% (C) R. Das, Stanford, HHMI, 2023
if ~exist('REMOVE_SINGLETS','var') REMOVE_SINGLETS = 0; end;
if all(structure=='x')
    structure = '';
    return
end
structure = lower(structure);
structure = strrep(structure,',','');

bps = convert_structure_to_bps_v3( structure);
if REMOVE_SINGLETS;    bps = remove_singlet_bps( bps ); end

structure = convert_bps_to_structure_v3( bps, length(structure), 1 );
