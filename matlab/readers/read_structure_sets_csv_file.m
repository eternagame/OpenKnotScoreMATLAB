function [x,structure_tags, structure_sets,structure_map] = read_structure_sets_csv_file( structure_sets_csv_file, ordered_sequences )
% [x,structure_tags, structure_sets,structure_map] = read_structure_sets_csv_file( structure_sets_csv_file, ordered_sequences );
%
% Inputs
%  structure_sets_csv_file = csv file with columns like "_mfe" holding
%                        structure predictions from different packages in dot bracket notation.
%  ordered_sequences = [Optional] list of sequences. If provided, structures read
%          in from .csv file will be reordered based on sequence column to match ordering in
%          ordered_sequence
%
% Outputs
% x = csv file in MATLAB Table object
% structure_tags = tags for each structure/mfe column
% structure_sets = cell of cell of strings of predicted structures
% structure_map = [Ndesign x Nres x Npackage] matrix of 0,1 for
%               paired/unpaired in each package structure prediction
%
% (C) R. Das, HHMI/Stanford University 2023.

x = readtable(structure_sets_csv_file);
structure_tags = {}; structure_sets = {}; count = 0;
for n = 1:length(x.Properties.VariableNames);
    tag = x.Properties.VariableNames{n};
    structures = table2cell(x(:,n));
    if contains(tag,'eternafold_threshknot_'); continue; end;
    if contains(tag,'sequence'); continue; end;
    is_structure_col = 1;
    for i = 1:size(x,1) 
        if ~ischar(structures{i}) | (~contains(structures{i},'.') & ~contains(structures{i},'(') & ~contains(structures{i},'x'))
            is_structure_col = 0;
        end
    end
    if ~is_structure_col; continue; end;
    count = count + 1;
    structure_tags{count} = strrep(strrep(tag,'__mfe',''),'_mfe','');
    structure_sets{count} = structures;
    fprintf( 'Sanitizing %d structures for %s...\n',size(x,1),structure_tags{count});
    for i = 1:size(x,1) % loop over designs
        structure_sets{count}{i} = sanitize_structure( structure_sets{count}{i} );
    end
end


% may need to do a reordering if sequence order in structure file does not
% match sequences 
if exist('ordered_sequences','var') & length( ordered_sequences ) > 0
    assert(length(unique(ordered_sequences))==length(ordered_sequences)); % check uniqueness
    sequence_col = find(strcmp(x.Properties.VariableNames,'sequence'));
    assert(~isempty(sequence_col));
    sequences = table2cell(x(:,sequence_col));
    assert(length(unique(sequences))==length(sequences)); % check uniqueness
    reorder = [];
    for i = 1:length(ordered_sequences)
        idx = find( strcmp(sequences, ordered_sequences{i} ));
        assert( ~isempty(idx));
        reorder(i) = idx;
    end
    for n = 1:length(structure_sets)
        structure_sets{n} = structure_sets{n}(reorder);
    end
end

structure_map = get_structure_map( structure_sets );
