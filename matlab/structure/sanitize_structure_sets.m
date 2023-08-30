function  structure_sets = sanitize_structure_sets(structure_sets, structure_tags, REMOVE_SINGLETS, PARALLELIZE);
% structure_sets = sanitize_structure_sets(structure_sets, structure_tags);
%
% Inputs
%  structure_sets = cell of cell of strings of predicted structures
%  structure_tags = tags for each structure/mfe column
%  REMOVE_SINGLETS = remove singlet base pairs (Default 0)
%  PARALLELIZE = use parfor to parallelize (Default 0).
%                Note: this doesn't seem worth it unless
%                you have like 100,000 sequences due to long
%                set up time.
%
% Outputs
%  structure_sets = cell of cell of strings of predicted structures,
%                      sanitized
%
% (C) R. Das, Stanford University, HHMI
if ~exist('REMOVE_SINGLETS','var') REMOVE_SINGLETS = 0; end;
if ~exist('PARALLELIZE','var') PARALLELIZE = 0; end;
tic
structure_tags_defined = exist( 'structure_tags','var');

if PARALLELIZE
    parfor count = 1:length(structure_sets)
        if structure_tags_defined; fprintf( 'Sanitizing %d structures for %s...\n',length(structure_sets{count}),structure_tags{count}); end
        for i = 1:length(structure_sets{count}) % loop over designs
            structure_sets{count}{i} = sanitize_structure( structure_sets{count}{i}, REMOVE_SINGLETS);
        end
    end
else
    for count = 1:length(structure_sets)
        if structure_tags_defined; fprintf( 'Sanitizing %d structures for %s...\n',length(structure_sets{count}),structure_tags{count}); end
        for i = 1:length(structure_sets{count}) % loop over designs
            structure_sets{count}{i} = sanitize_structure( structure_sets{count}{i}, REMOVE_SINGLETS);
        end
    end
end
toc