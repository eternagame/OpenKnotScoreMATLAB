function make_heatmap_with_structures_for_design( idx, r_norm, structure_map, structure_sets,structure_tags, pkg_sort_idx, headers, sequences, BLANK_OUT5, BLANK_OUT3, conditions, create_figure_window );
% make_heatmap_with_structures_for_design( idx, r_norm, structure_map, structure_sets,structure_tags, pkg_sort_idx, headers, sequences, BLANK_OUT5, BLANK_OUT3, conditions );
%
% Inputs
%  idx = [integer] index of design for which to show heatmap
%  r_norm = [Ndesign x Nres x Nconditions] Reactivity matrix, normalized.
%  structure_map = [Ndesign x Nres x Npackages] 0/1 map
%            of paired/unpaired for each predicted structure
%  structure_sets = [Npackages x Ndesign] cell of cell of strings of predicted structures
%  structure_tags = cell of string, name of each package
%  pkg_sort_idx = permutation of packages (e.g., [2, 4, 1, 3]) to order
%             packages (Give empty set [] to order by SHAPE)
%  headers = cell of Ndesign strings describing each design (titles for
%  plot)
%  sequences = cell of sequences for Ndesigns
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues 
%  conditions = [Nconditions] cell of strings with labels for each of the
%       experimental conditions in r_norm
%  create_figure_window = Create new figure window (Default 1). If 0, 
%                           suppresses the re-focusing of MATLAB to the
%                           window and allows background image rendering.
%
% (C) R. Das, HHMI/Stanford University 2023.

if ~exist('conditions','var') conditions = {'SHAPE, no Mg2+','SHAPE, +Mg2+'}; end
if ~exist('create_figure_window','var') create_figure_window = 1; end

labels = {'sequence'};
reverse_profile_order = [size(r_norm,3):-1:1];
d = shiftdim(r_norm(idx,:,reverse_profile_order),1)';
conditions{1} = ['*',conditions{1}];
labels = [labels,conditions(reverse_profile_order)];
Ndata = size(r_norm,3);

if isempty(structure_map)
    structure_map_for_idx = get_structure_map( structure_sets, idx);
else
    structure_map_for_idx = structure_map(idx,:,:);
end

if isempty(pkg_sort_idx)
    r_norm_for_scoring = r_norm(idx,:,1); % just use first data set
    % r_norm_for_scoring = mean(r_norm(idx,:,:),3)
    for n = 1:size(structure_map_for_idx,3);
        eterna_scores(n) = calc_eterna_score_classic( r_norm_for_scoring, structure_map_for_idx(:,:,n), BLANK_OUT5, BLANK_OUT3);
    end
    [~,pkg_sort_idx] = sort(eterna_scores);
    for n = 1:size(structure_map_for_idx,3);
        structure_tags{n} = sprintf( '%s EternaClassic %5.1f', structure_tags{n}, eterna_scores(n) );
    end
end

best_structure = '';
best_fit_idx = find( contains( structure_tags, 'best_fit' ) );
[~,sortidx] = sort(-pkg_sort_idx(best_fit_idx));
if ~isempty(sortidx); best_structure = structure_sets{ best_fit_idx(sortidx(1)) }{idx}; end;

show_structures=repmat({''},1,Ndata);
for n = pkg_sort_idx(end:-1:1)
    d = [d; 0.4*structure_map_for_idx(:,:,n)];
    show_structures = [show_structures,structure_sets{n}{idx}];

    mfe_tag = structure_tags{n};
    if length(best_structure) > 0 & strcmp(best_structure,structure_sets{n}{idx}); mfe_tag = ['*** ',mfe_tag]; 
    elseif exist('eterna_scores','var') & eterna_scores(n) > eterna_scores(pkg_sort_idx(end)) - 5;
        mfe_tag = ['*',mfe_tag]; 
    end
    labels = [labels,mfe_tag];
end
imagesc(d,[-2 2]);
colormap(redwhiteblue(-1,1));

for i = 1:length(show_structures)
    show_structure = show_structures{i};
    for n = 1:length(show_structure)
        if ~strcmp(show_structure(n),'.')
            text(n,i,show_structure(n),'interpreter','none','HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
end

sequence = sequences{idx};
for n = 1:length(sequence)
    text(n,0,sequence(n),'interpreter','none','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',8);
end

N = size(d,2);
% gray out regions where there shouldn't be data.
% gray out flanking regions where primers bind
rectangle('Position',[0 0.5 BLANK_OUT5+0.5 Ndata],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
rectangle('Position',[N-BLANK_OUT3+0.5 0.5 BLANK_OUT3+0.5 Ndata],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
% gray out profiles for which we are missing data (signaled by nan)
nan_profile_idx = find(isnan(d(:,BLANK_OUT5+1)));
for i = nan_profile_idx'
    rectangle('Position',[0 i-0.5 N+0.5 1],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
end
% gray out G and U for DMS data, since DMS mainly hits A and C?
for i = find(contains(conditions,'DMS') & ~contains(conditions,'DMS_N7'))
    gu_idx = union( strfind(sequence,'G'), strfind(sequence,'U') );
    for k = gu_idx
        rectangle('Position',[k-0.5 reverse_profile_order(i)-0.5 1 1],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
    end
end
% gray out A, C, and U for DMS_N7
for i = find(contains(conditions,'DMS_N7'))
    acu_idx = find(sequence=='A' | sequence=='C' | sequence == 'U');
    for k = acu_idx
        rectangle('Position',[k-0.5 reverse_profile_order(i)-0.5 1 1],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
    end
end

% Axis labels
xlabel( 'Position');
h=title(strsplit(strrep(headers{idx},sprintf('\t'),' '),';'));
set(h,'interpreter','none')
set(gca,'ytick',[0:size(d,1)],'yticklabels',labels,'xtick',[0:10:N],'xticklabel',[0:10:N]);
set(gca,'fontweight','bold','TickLabelInterpreter','none','tickdir','out');
box off
ylim([-1 size(d,1)+0.5])
xlim([0.5 size(d,2)+0.5]);
%pause;

h=colorbar();
colorTitleHandle = get(h,'Title');
set(colorTitleHandle,'String','SHAPE')