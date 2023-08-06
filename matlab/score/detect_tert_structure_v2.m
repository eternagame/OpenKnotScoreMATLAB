function protected_loop2_score = detect_tert_structure_v2( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3, openknot_info_structs, headers, r_norm_err, SHOW_PLOT   );
% detect_tert_structure_v2( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3, structure_sets_extended, structure_tags_extended, headers, r_norm_err, SHOW_PLOT   );
%
% (C) R. Das, HHMI/Stanford 2023

if ~exist( 'SHOW_PLOT'); SHOW_PLOT = 1; end;
if isempty( good_idx ); good_idx = [1:length(sequences)]; end;
if SHOW_PLOT; set(figure(15),'color','white','position',[100 500 1200 250]); clf; end;
N = length(sequences{1});
which_pos = [BLANK_OUT5+1:N-BLANK_OUT3];
rprof_err = [];
protected_loop2_scores = [];
for q = 1:length(good_idx)
    i = good_idx(q);
    sequence  = sequences{i};
    structures = openknot_info_structs{i}.best_fit.structures;
    % find non-empty structure:
    for m = 1:length(structures)
        structure_full = openknot_info_structs{i}.best_fit.structures{m};
        structure_tag = openknot_info_structs{i}.best_fit.tags{m};
        if contains(structure_full, '('); break; end;
    end
    structure = remove_pk( structure_full );
    rprof = r_norm(i,:);
    if exist('r_norm_err','var'); rprof_err = r_norm_err(i,:); end;
    protected_loop = intersect(intersect( find(rprof<0.25),strfind(structure,'.')),which_pos);
    protected_loop2 = intersect(protected_loop,protected_loop+1);
    protected_loop2 = union(protected_loop2,protected_loop2-1);
    exposed_pair = intersect(setdiff( find(rprof>0.25),strfind(structure,'.')),which_pos);
    protected_loop2_score(q) = 100*length(protected_loop2)/length(which_pos);
    if SHOW_PLOT
        plot( rprof, 'k','linew',2 ); hold on;
        plot([0 length(sequence)],[0 0],'k')
        plot([0 length(sequence)],0.25*[1 1],'k:')
        for m = 1:length(structure_full)
            text(m,-0.2,structure_full(m),'horizontalalign','center');
            text(m,-0.4,sequence(m),'horizontalalign','center');
        end
        text(length(structure)+1,-0.2,structure_tag,'horizontalalign','left','interpreter','none')
        plot_res = which_pos; clr = [0 0 0]; plot_symbols(rprof,rprof_err,plot_res, clr);
        plot_res = protected_loop; clr = [0.2 0.3 0.2]; plot_symbols(rprof,rprof_err,plot_res, clr);
        plot_res = protected_loop2;clr = [0 0.7 0]; plot_symbols(rprof,rprof_err,plot_res, clr);
        plot_res = exposed_pair;   clr = [1 0 0];  plot_symbols(rprof,rprof_err,plot_res, clr);
        ylim([-1 3])
        hold off
        box off
        title( headers{i},'interpreter','none')
        if length(good_idx)>1 & i ~= good_idx(end); pause;end;
    end
end


%%%%%%%
function structure = remove_pk( structure_full );
structure=structure_full;
pk_pos = setdiff( [1:length(structure)], [strfind(structure,'.'),strfind(structure,'('),strfind(structure,')')] );
structure(pk_pos) = '.';

%%%%%%%
function plot_symbols(rprof,rprof_err,plot_res,clr); 

plot(plot_res,rprof(plot_res),'o','color',clr,'markerfacecolor',clr);

if length(rprof_err) == 0; return; end;
for i = 1:length(plot_res)
    plot( plot_res(i)*[1 1], rprof(plot_res(i)) + rprof_err(plot_res(i))*[-1 1],'color',clr,'linew',1.5 );
end

