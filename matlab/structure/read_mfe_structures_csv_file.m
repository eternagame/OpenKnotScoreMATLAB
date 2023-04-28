function [x,mfe_tags, mfe_structures,mfe_structure_map] = read_mfe_structures_csv_file( mfe_structures_csv_file, ordered_sequences );
% [x,mfe_tags, mfe_structures,mfe_structure_map] = read_mfe_structures_csv_file( mfe_structures_csv_file, ordered_sequences );
%
% Inputs
%  mfe_structures_csv_file = csv file with columns like "_mfe" holding
%                        structure predictions from different packages in dot bracket notation.
%  ordered_sequences = [Optional] list of sequences. If provided, structures read
%          in from .csv file will be reordered based on sequence column to match ordering in
%          ordered_sequence
%
% Outputs
% x = csv file in MATLAB Table object
% mfe_tags = tags for each structure/mfe column
% mfe_structures = cell of cell of strings of predicted structures
% mfe_structure_map = [Ndesign x Nres x Npackage] matrix of 0,1 for
%               paired/unpaired in each package structure prediction
%
% (C) R. Das, HHMI/Stanford University 2023.

x = readtable(mfe_structures_csv_file);
mfe_tags = {}; count = 0;
for n = 1:length(x.Properties.VariableNames);
    tag = x.Properties.VariableNames{n};
    structures = table2cell(x(:,n));
    if contains(tag,'_mfe') | ...
            contains(tag,'structure') | ...
            ( contains(tag,'hfold') & ...
            ~contains(tag,'time')  | ...
            strcmp(tag,'eterna_nupack') ) | ...
            (ischar(structures{1}) & contains(structures{1},'('))
        if ~ischar( structures{1}); continue; end;
        if contains(tag,'eternafold_threshknot_'); continue; end;
        count = count + 1;
        mfe_tags{count} = strrep(strrep(tag,'__mfe',''),'_mfe','');
        mfe_structures{count} = structures;
        fprintf( 'Sanitizing %d structures for %s...\n',size(x,1),mfe_tags{count});
        for i = 1:size(x,1) % loop over designs
            mfe_structures{count}{i} = sanitize_structure( mfe_structures{count}{i} );
        end
    end
end

% may need to do a reordering if sequence order in structure file does not
% match sequences 
if exist('ordered_sequences','var')
    assert(length(unique(ordered_sequences))==length(ordered_sequences)); % check uniqueness
    sequence_col = find(strcmp(x.Properties.VariableNames,'sequence'));
    assert(~isempty(sequence_col));
    sequences = table2cell(x(:,sequence_col));
    assert(length(unique(sequences))==length(sequences)); % check uniqueness
    reorder = [];
    for i = 1:length(sequences)
        idx = find( strcmp(ordered_sequences, sequences{i} ));
        assert( ~isempty(idx));
        reorder(i) = idx;
    end
    for n = 1:length(mfe_structures)
        mfe_structures{n} = mfe_structures{n}(reorder);
    end
end

mfe_structure_map = get_mfe_structure_map( mfe_structures );
