function [x,structure_tags, structure_sets,structure_map,found_structure_idx] = read_structure_sets_csv_file( structure_sets_csv_file, ordered_sequences, sanitize_structures )
% [x,structure_tags, structure_sets,structure_map,found_structure_idx] = read_structure_sets_csv_file( structure_sets_csv_file, ordered_sequences, sanitize_structures );
%
% Inputs
%  structure_sets_csv_file = csv file with columns like "_mfe" holding
%                        structure predictions from different packages in dot bracket notation.
%  ordered_sequences = [Optional] list of sequences. If provided, structures read
%          in from .csv file will be reordered based on sequence column to match ordering in
%          ordered_sequence
%  sanitize_structures = sanitize dot-bracket structure (Default 1)
%
% Outputs
% x = csv file in MATLAB Table object
% structure_tags = tags for each structure/mfe column
% structure_sets = cell of cell of strings of predicted structures
% structure_map = [Ndesign x Nres x Npackage] matrix of 0,1 for
%               paired/unpaired in each package structure prediction
% found_structure_idx = which ordered_sequences found a match
%
% (C) R. Das, HHMI/Stanford University 2023.
warning('off');
if ~exist('sanitize_structures','var') sanitize_structures = 1; end;
tic
x = readtable(structure_sets_csv_file);
structure_tags = {}; structure_sets = {}; count = 0;
fprintf('\n');
for n = 1:length(x.Properties.VariableNames);
    tag = x.Properties.VariableNames{n};
    structures = table2cell(x(:,n));
    if contains(tag,'eternafold_threshknot_'); continue; end;
    if contains(tag,'sequence'); continue; end;
    is_structure_col = 1;
    for i = 1:size(x,1) 
        if strcmp(structures{i},'ERR'); structures{i}=repmat('.',1,length(structures{1})); end;
        if ~ischar(structures{i}) | (~contains(structures{i},'.') & ~contains(structures{i},'(') & ~contains(structures{i},'x') )
            is_structure_col = 0;
            break;
        end
    end
    if ~is_structure_col; continue; end;
    count = count + 1;
    structure_tags{count} = strip(strrep(strrep(tag,'__mfe',''),'_mfe',''),'_');
    structure_sets{count} = structures;
    if sanitize_structures
        fprintf( 'Sanitizing %d structures for %s...\n',size(x,1),structure_tags{count});
        for i = 1:size(x,1) % loop over designs
            structure_sets{count}{i} = sanitize_structure( structure_sets{count}{i} );
        end
    end
end
fprintf('\n');
toc

% may need to do a reordering if sequence order in structure file does not
% match sequences 
found_structure_idx = [];
if exist('ordered_sequences','var') & length( ordered_sequences ) > 0
    tic
    fprintf( 'Matching into ordered sequences...\n')
    assert(length(unique(ordered_sequences))==length(ordered_sequences)); % check uniqueness
    sequence_col = find(strcmp(x.Properties.VariableNames,'sequence'));
    assert(~isempty(sequence_col));
    sequences = table2cell(x(:,sequence_col));
    assert(length(unique(sequences))==length(sequences)); % check uniqueness
    count = 0;
    structure_sets_from_csv = structure_sets;
    blank_sequences = {};
    for i = 1:length(ordered_sequences)
        blank_sequences{i} = repmat('.',1,length(ordered_sequences{i}));
    end
    structure_sets = {}; 
    for n = 1:length(structure_tags); structure_sets{n} = blank_sequences; end
    reorder = [];
    % need to use map/dictionary for speed.
    d = containers.Map(sequences,[1:length(sequences)]);
    for i = 1:length(ordered_sequences)
        if d.isKey( ordered_sequences{i} )
            reorder( d(ordered_sequences{i}) ) = i;
        end
    end
    for n = 1:length(structure_sets)
        structure_sets{n}(reorder) = structure_sets_from_csv{n};
    end
    fprintf( 'Matched %d csv sequences into %d out of %d ordered sequences\n',...
        length(structure_sets_from_csv{1}),length(reorder),length(ordered_sequences));
    found_structure_idx = sort(reorder);
    toc
end
tic
fprintf( 'Getting structure map...\n')
structure_map = get_structure_map( structure_sets );
toc