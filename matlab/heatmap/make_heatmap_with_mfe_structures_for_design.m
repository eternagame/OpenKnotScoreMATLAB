function make_heatmap_with_mfe_structures_for_design( idx, r_norm, mfe_structure_map, mfe_structures,mfe_tags, pkg_sort_idx, headers, sequences, BLANK_OUT5, BLANK_OUT3, tags, create_figure_window );
% make_heatmap_with_mfe_structures_for_design( idx, r_norm, mfe_structure_map, mfe_structures,mfe_tags, pkg_sort_idx, headers, sequences, BLANK_OUT5, BLANK_OUT3, tags );
%
% Inputs
%  idx = [integer] index of design for which to show heatmap
%  r_norm = [Ndesign x Nres x Nconditions] Reactivity matrix, normalized.
%  mfe_structure_map = [Ndesign x Nres x Npackages] 0/1 map
%            of paired/unpaired for each predicted structure
%  mfe_structures = [Npackages x Ndesign] cell of cell of strings of predicted structures
%  mfe_tags = cell of string, name of each package
%  pkg_sort_idx = permutation of packages (e.g., [2, 4, 1, 3]) to order
%             packages
%  headers = cell of Ndesign strings describing each design (titles for
%  plot)
%  sequences = cell of sequences for Ndesigns
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues 
%  tags = [Nconditions] cell of strings with labels for each of the
%       experimental conditions in r_norm
%  create_figure_window = Create new figure window (Default 1). If 0, 
%                           suppresses the re-focusing of MATLAB to the
%                           window and allows background image rendering.
%
% (C) R. Das, HHMI/Stanford University 2023.

if ~exist('tags','var') tags = {'SHAPE, no Mg2+','SHAPE, +Mg2+'}; end
if ~exist('create_figure_window','var') create_figure_window = 1; end

labels = {'sequence'};
d = shiftdim(r_norm(idx,:,:),1)';
labels = [labels,tags];
Ndata = size(r_norm,3);

if isempty(pkg_sort_idx)
    %[all_corr_coef, pkg_sort_idx] = get_corr_coeff( mean(r_norm,3), mfe_structure_map, idx, mfe_tags, BLANK_OUT3, BLANK_OUT5);
    for n = 1:size(mfe_structure_map,3);
        eterna_scores(n) = calc_eterna_score_classic( mean(r_norm(idx,:,:),3), squeeze(mfe_structure_map(idx,:,n)), BLANK_OUT5, BLANK_OUT3);
    end
    [~,pkg_sort_idx] = sort(eterna_scores);
end

best_structure = '';
best_fit_idx = find( contains( mfe_tags, 'best_fit' ) );
[~,sortidx] = sort(-pkg_sort_idx(best_fit_idx));
if ~isempty(sortidx); best_structure = mfe_structures{ best_fit_idx(sortidx(1)) }{idx}; end;

show_structures=repmat({''},1,Ndata);
for n = pkg_sort_idx(end:-1:1)
    d = [d; 0.4*mfe_structure_map(idx,:,n)];
    show_structures = [show_structures,mfe_structures{n}{idx}];

    mfe_tag = mfe_tags{n};
    if length(best_structure) > 0 & strcmp(best_structure,mfe_structures{n}{idx}); mfe_tag = ['*** ',mfe_tag]; end;
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
    text(n,0,sequence(n),'interpreter','none','HorizontalAlignment','center','VerticalAlignment','middle');
end

N = size(d,2);
% gray out regions where there shouldn't be data.
rectangle('Position',[0 0.5 BLANK_OUT5+0.5 Ndata],'EdgeColor','none','FaceColor',[0.5 0.5 0.5]);
rectangle('Position',[N-BLANK_OUT3+0.5 0.5 BLANK_OUT3+0.5 Ndata],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
xlabel( 'Position');
h=title(headers(idx));
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