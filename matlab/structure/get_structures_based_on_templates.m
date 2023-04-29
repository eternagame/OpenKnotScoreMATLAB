function [structures,idx_with_template] = get_structures_based_on_templates( sequences,templates_file,baseline_structure,template_sequence_tag,template_structure_tag,template_design_id_tag, ids);
%
% [structures,idx_with_template] = get_structures_based_on_templates( sequences,templates_file,baseline_structure,template_sequence_tag,template_structure_tag,template_design_id_tag, ids);
%
% Manually process a .csv with template sequences and (dot-bracket)
% structures and output structures for all design sequences
%
% Optional: if .csv file has manually annotated, e.g., Eterna ID's
%  associated with templates, check match with automatic annotation.
%
% Inputs
%  sequences = {cell of Ndesign strings} sequences for which we want
%                 template-based structures
%  templates_file  = .csv file with columns that include sequence of
%           templates, structure of templates.
%
% Optional inputs
%  baseline_structure = string of length Nres (dot-parens) with baseline
%     structure, which might include, e.g., flanking hairpins.
%  template_sequence_tag = column name in .csv file with template sequences, e.g.
%                   'Sequence'.
%  template_structure_tag = column name in .csv file with template sequences, e.g.
%                   'DotBracketFromUTEP'.
%  template_design_id_tag = column name in .csv file with manually curated design 
%                         IDs. This and later columns are assumed to
%                         potentially hold 0, 1, or more design IDs
%  ids = [Ndesign integers] design IDs for input sequences (allows for cross-check)
%
% Outputs
%  structures = {cell of Ndesign strings} dot-parens structures with
%           template structures matched in.
%  idx = list of integers: index of sequences for which templates were
%           found.
%
% (C) R. Das, Stanford, HHMI 2023
warning('off');
x = readtable( templates_file);
fprintf('Read in file %s\n',templates_file);
fprintf('Total number of templates to match: %d\n',size(x,1));

%% check template
template_sequences  = getfield_from_description(x, template_sequence_tag);
template_structures = getfield_from_description(x, template_structure_tag);
design_ids = []; 

N = length(sequences{1});
if ~exist( 'baseline_structure','var') | isempty(baseline_structure)
    baseline_structure = repmat('.',1,N);
end

templates_in_sequence = {};
template_with_sequence = [];
for i = 1:length(sequences)
    sequence = sequences{i};
    structure = baseline_structure;
    template_in_sequence = [];
    for j = 1:length(template_sequences)
        idx = strfind(sequence, upper(template_sequences{j}));
        if isempty(idx); continue; end;
        template_in_sequence = [template_in_sequence,j];
        template_with_sequence(j) = i;
        for m = idx;
            template_structure = template_structures{j};
            structure(m + [0:length(template_structure)-1]) = sanitize_structure(strrep(template_structure,'-','.'));
        end
    end
    if length(template_in_sequence) == 0; % no templates found in here
        % structure = '';  % this is leading to issues in CSV readin
        structure = repmat('.',1,N);
    end
    structures{i} = structure;
    templates_in_sequence{i} = template_in_sequence;
end

idx_with_template = find(cellfun(@length,templates_in_sequence));

fprintf('\nTotal number of sequences with templates identified: %d\n',length(idx_with_template));

if ~exist('ids','var') | isempty(ids); return; end;
if ~exist('template_design_id_tag','var') | isempty(template_design_id_tag); return; end;


% Check by template -- did we find matches for all the ones that were also
% found manually?
fprintf('\n')
manual_template_with_sequence_idx = [];
manual_design_ids = convert_to_numbers( getfield_from_description(x,template_design_id_tag));
manual_template_with_sequence_idx = find(~isnan(manual_design_ids));
template_with_sequence_idx = find(template_with_sequence>0);

i = setdiff(template_with_sequence_idx,manual_template_with_sequence_idx);
fprintf('Template found by automatic scan but not manually: %d\n',length(i));
if length(i)>0; fprintf(' %d',i);fprintf('\n'); x(i,:)
end;

i = setdiff(manual_template_with_sequence_idx,template_with_sequence_idx);
fprintf('Template found manually but not by automatic scan: %d\n',length(i))
if length(i)>0; fprintf(' %d',i);fprintf('\n');x(i,:);end;

% Check by sequence -- did we find template matches for all the ones that were also
% found manually?
for i = 1:length(sequences); manual_templates_in_sequence{i} = []; end;

manual_design_col = find(strcmp(x.Properties.VariableDescriptions, template_design_id_tag));
for k = (manual_design_col):length(x.Properties.VariableDescriptions)
    col = convert_to_numbers(getfield_from_description(x, x.Properties.VariableDescriptions{k}));
    for i = 1:length(template_sequences)
        id = col(i);
        if ~isnan(id);  
            idx = find(ids==id);
            if isempty( idx ); fprintf( 'Was not able to find manually assigned design id %d in ids!\n',id); continue;end;
            manual_templates_in_sequence{idx} = [manual_templates_in_sequence{idx},i]; 
        end
    end
end

manual_templates_in_sequence_idx = find(cellfun(@length,manual_templates_in_sequence));
templates_in_sequence_idx = find(cellfun(@length,templates_in_sequence));


i = setdiff(templates_in_sequence_idx,manual_templates_in_sequence_idx);
fprintf('Sequence with template found by automatic scan but not manually: %d\n%d\n', length(i));
if length(i)>0; fprintf(' %d',i);fprintf('\n'); end;
i = setdiff(manual_templates_in_sequence_idx,templates_in_sequence_idx);
fprintf('Sequence with template found manually but not in automatic scan: %d\n%d\n',length(i));
if length(i)>0; fprintf(' %d',i);fprintf('\n');end;

%%%%%%%%%%%%%%%%%
function f = getfield_from_description( x, desc)
idx = find(strcmp(x.Properties.VariableDescriptions,desc));
f = getfield(x,x.Properties.VariableNames{idx});

%%%%%%%%%%%%%%%%%
function manual_design_ids_new = convert_to_numbers( manual_design_ids );
if isnumeric( manual_design_ids); 
    manual_design_ids_new = manual_design_ids;
    return; 
end;

manual_design_ids_new = [];
for i = 1:length(manual_design_ids)
    id = str2num( manual_design_ids{i} );
    if isempty(id); manual_design_ids_new(i) = NaN; else manual_design_ids_new(i) = id; end;
end

