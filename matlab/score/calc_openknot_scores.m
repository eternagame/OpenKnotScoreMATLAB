function [openknot_info_structs,structures,openknot_scores,eterna_scores,crossed_pair_scores,crossed_pair_quality_scores,best_cc,all_eterna_classic_scores] = calc_openknot_scores( r_norm, structure_sets, good_idx, mfe_tags, BLANK_OUT3, BLANK_OUT5, headers, make_plot, REMOVE_SINGLETS);
% [openknot_info_structs,structures,openknot_scores,eterna_scores,crossed_pair_scores,crossed_pair_quality_scores,structures,best_cc,all_eterna_classic_scores] = calc_openknot_scores( r_norm, structure_sets, good_idx, mfe_tags, BLANK_OUT3, BLANK_OUT5, headers, make_plot, REMOVE_SINGLETS);
%
% Inputs:
%  r_norm = [Ndesign x Nres] Reactivity matrix, normalized.
%               normalized to go from 0 to 1 (~90th percentile)
%  structure_sets = [Npackages x Ndesign] cell of cell of strings of predicted structures
%  good_idx = [Nidx] index of designs for which to show heatmap
%  mfe_tags = cell of string, name of each package
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues
%  headers = cell of Ndesign strings describing each design (titles for
%  plot)
%  REMOVE_SINGLETS = remove singlet base pairs (Default 1) -- note time
%                      consuming!
%
% Outputs
%  openknot_info_structs = cell of 'structs' with more detailed information,
%       including alternative structures and their scores. Structs include
%       all of following
%  openknot_scores = [Nidx] OpenKnot score = 1/2( Eterna + Crossed_pair)
%                         (averaged over all structures with EternaScore
%                         within 5 of best.) Range: 0-100.
%  eterna_scores   = ([Nidx] classic Eterna score
%                         (averaged over all structures with Eterna Score
%                         within 5 of best.) Range: 0-100.
%  crossed_pair_scores = [Nidx] num of crossed_pair residues with
%                         SHAPE<0.5, normalized to 0.7 x region length.
%                         (averaged over all structures with EternaScore
%                         within 5 of best.) Range: 0-100.
%  crossed_pair_quality_scores = [Nidx] num of crossed pair residues with
%                         SHAPE < 0.5, normalized to total number of
%                         predicted crossed pairs Range: 0-100.
%  best_structures = (cell with Nidx strings) very best fit structure
%  best_cc = [Nidx] correlation coefficients for best fit structure
%
% (C) Rhiju Das, Stanford, HHMI, 2023

if length(good_idx) == 0; good_idx = [1:size(r_norm,1)]';end;
if ~exist('make_plot','var'); make_plot = 1; end;
if ~exist('REMOVE_SINGLETS','var'); REMOVE_SINGLETS = 1; end;

openknot_info_structs = repmat({struct()},1,length(good_idx));
structures = repmat({''},1,length(good_idx));
openknot_scores = zeros(length(good_idx),1);
eterna_scores = zeros(length(good_idx),1);
crossed_pair_scores = zeros(length(good_idx),1);
best_cc = zeros(length(good_idx),1);
crossed_pair_quality_scores = zeros(length(good_idx),1);
all_eterna_classic_scores = repmat({[]},1,length(good_idx));

% parfor incurs long set up time. better to pre-remove singlets with 
%  sanitize_structure_sets and then run calc_openknot_scores with
%  REMOVE_SINGLETS = 0.
% parfor i = 1:length(good_idx) 
for i = 1:length(good_idx)
    idx = good_idx(i);
    structure_sets_for_idx = {};
    for n = 1:length(structure_sets); structure_sets_for_idx{n}=structure_sets{n}(idx); end;
    [openknot_info_struct, eterna_classic_scores] = calc_openknot_score( r_norm(idx,:,:), structure_sets_for_idx, 1, mfe_tags, BLANK_OUT3, BLANK_OUT5, headers(idx), make_plot, REMOVE_SINGLETS );

    openknot_info_structs{i} = openknot_info_struct;
    structures{i} = openknot_info_struct.best_fit.structures{1};
    openknot_scores(i) = openknot_info_struct.score.openknot_score;
    eterna_scores(i) = openknot_info_struct.score.eterna_classic_score;
    crossed_pair_scores(i) = openknot_info_struct.score.crossed_pair_score;
    crossed_pair_quality_scores(i) = openknot_info_struct.score.crossed_pair_quality_score;
    best_cc(i) = openknot_info_struct.score.best_cc;
    all_eterna_classic_scores{i} = eterna_classic_scores;
end

function [openknot_info_struct, eterna_classic_scores] = calc_openknot_score( r_norm, structure_sets, idx, mfe_tags, BLANK_OUT3, BLANK_OUT5, headers, make_plot, REMOVE_SINGLETS )

data = r_norm(idx,:,1);

for n = 1:length(mfe_tags)
    structure = structure_sets{n}{idx};
    structure = strrep(structure,'x','.');
    if REMOVE_SINGLETS; structure_sets{n}{idx} = sanitize_structure( structure, 1); end;
    structure_sets{n}{idx} = structure;
end
structure_map_for_idx = get_structure_map( structure_sets, idx );

[cc,cc_sort_idx] = get_corr_coeff( r_norm(idx,:,:), structure_map_for_idx, 1, mfe_tags, BLANK_OUT3, BLANK_OUT5, 'Pearson', (make_plot-1) );

% eterna score classic:
for n = 1:length(mfe_tags)
    pred = squeeze(structure_map_for_idx(:,:,n));
    eterna_classic_scores(n) = calc_eterna_score_classic( data, pred, BLANK_OUT5, BLANK_OUT3);
    [crossed_pair_scores(n),crossed_pair_quality_scores(n)] = calc_crossed_pair_score(data, structure_sets{n}{idx}, BLANK_OUT5, BLANK_OUT3, REMOVE_SINGLETS);
end
if make_plot
    plot( cc, eterna_classic_scores,'.')
    xlabel('Correlation coefficient');
    ylabel( 'Eterna score (classic)');
    title(headers{idx},'interp','none');
end

% use eterna score instead of correlation coefficient...
eterna_classic_score = max(eterna_classic_scores);
best_model_idx = find( eterna_classic_score == eterna_classic_scores);

crossed_pair_score = max(crossed_pair_scores(best_model_idx));
best_model_idx = find( eterna_classic_score == eterna_classic_scores & crossed_pair_scores == crossed_pair_score);

structure = structure_sets{best_model_idx(1)}{idx};
best_cc = cc(best_model_idx(1));

% try ensembling any models that are close, based on Eterna Score
eterna_score_band = 5;
best_model_idx = find( eterna_classic_scores >= (eterna_classic_score - eterna_score_band) );
[~,sortidx] = sort(eterna_classic_scores(best_model_idx));
best_model_idx = best_model_idx( sortidx(end:-1:1) );
eterna_classic_score = mean( eterna_classic_scores( best_model_idx) );
crossed_pair_score   = mean( crossed_pair_scores( best_model_idx) );
crossed_pair_quality_score   = mean( crossed_pair_quality_scores( best_model_idx) );
% openknot_score = 0.5 * ( eterna_classic_score + crossed_pair_score );
openknot_score = 0.5 * ( eterna_classic_score + crossed_pair_quality_score );
best_struct_eterna_scores = eterna_classic_scores( best_model_idx );

best_structs = {};
best_struct_tags = {};
for n = best_model_idx;
    if make_plot; fprintf( '%30s %s\n',mfe_tags{n},structure_sets{n}{idx}); end;
    best_structs = [best_structs,{structure_sets{n}{idx}}];
    best_struct_tags = [best_struct_tags,{mfe_tags{n}}];
end

fprintf( 'OpenKnot: %5.2f Eterna classic: %5.2f  Crossed pair score: %5.2f  Crossed quality: %5.2f CC: %4.2f. Best model for %s. \n',...
    openknot_score, ...
    eterna_classic_score,...
    crossed_pair_score, ...
    crossed_pair_quality_score,...
    best_cc, headers{idx});

openknot_info_struct = struct();
openknot_info_struct.score = struct();
openknot_info_struct.score.openknot_score = openknot_score;
openknot_info_struct.score.eterna_classic_score = eterna_classic_score;
openknot_info_struct.score.crossed_pair_score = crossed_pair_score;
openknot_info_struct.score.crossed_pair_quality_score = crossed_pair_quality_score;
openknot_info_struct.score.best_cc = best_cc;
openknot_info_struct.best_fit = struct();
openknot_info_struct.best_fit.tags = best_struct_tags;
openknot_info_struct.best_fit.structures = best_structs;
openknot_info_struct.best_fit.eterna_classic_scores = best_struct_eterna_scores;




