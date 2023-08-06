function structure_map = get_structure_map( structure_sets, idx );
% structure_map = get_structure_map( structure_sets, idx );
%
% Inputs
%  structure_sets = cell of Nsets cells of Nstructures structure strings
%  idx = [Optional] which structures in a structure_set (cell of strings)
%               to look at. Rest will be zero-ed in output structure_map
%               Default: look at all Nstructuers structures.
%
% Output
%  structure_map = [Nstructures x Nres x Nsets] 0/1 for unpaired
%                        residues.
%
% (C) R. Das, HHMI, Stanford, 2023

structure_map = [];

%% NEED TO SET structure_map size based on maximum structure length in structure_sets.
Nstructures = 0;
Nres = 0;
for i = 1:length( structure_sets );
    structure_set = structure_sets{i};
    Nstructures = max(Nstructures,length(structure_set));
    design_indices = 1:length(structure_set);
    if exist('idx','var'); design_indices = idx; end;
    for j = design_indices
        Nres = max(Nres,length(structure_sets{i}{j}));
    end
end

structure_map = zeros(Nstructures,Nres,length(structure_sets));
for count = 1:length( structure_sets );
    structure_set = structure_sets{count};
    design_indices = 1:length(structure_set);
    if exist('idx','var'); design_indices = idx; end;

    % Quickly identify '.' across all structures
    structure_append = strjoin(pad(structure_set(design_indices),Nres),'');
    unpaired_pos = strfind(structure_append,'.');
    % note transpose needed to get from 1D unpaired_pos into the right
    % places in 2D
    structure_map_for_set = zeros(Nres,length(design_indices)); 
    structure_map_for_set(unpaired_pos) = 1;

    % move 'em in. 
    structure_map(design_indices,:,count) = structure_map_for_set';

    % old code -- takes longer for lots of structures
    %     for i = design_indices % loop over designs
    %          structure_map(i,strfind(structure_sets{count}{i},'.'),count) = 1;
    %     end
end
