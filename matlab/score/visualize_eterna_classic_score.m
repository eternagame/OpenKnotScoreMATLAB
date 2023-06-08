function visualize_eterna_classic_score( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3, structures, headers, crossed_pair_only );
% detect_tert_structure( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3, structures, headers, crossed_pair_only  );
%
% (C) R. Das, HHMI/Stanford 2023

set(figure(15),'color','white','position',[100 500 1200 250]); clf;
N = length(sequences{1});
which_pos_general = [BLANK_OUT5+1:N-BLANK_OUT3];
rprof_err = [];
if ~exist( 'crossed_pair_only', 'var') crossed_pair_only = 0; end;


for i = good_idx'
    sequence =sequences{i};
    structure = structures{i};
    rprof = r_norm(i,:);
    which_pos = which_pos_general;
    if crossed_pair_only
        bps = convert_structure_to_bps_v2(structure);
        % singlets are throwing things off.
        bps = remove_singlet_bps( bps );
        crossed_res = figure_out_which_bps_are_crossed( bps );
        which_pos = intersect(crossed_res, which_pos);
    end

    plot( rprof, 'k','linew',2 ); hold on;
    plot([0 length(sequence)],[0 0],'k')
    plot([0 length(sequence)],0.125*[1 1],'k:')
    for m = 1:length(structure)
        text(m,-0.2,structure(m),'horizontalalign','center');
    end
    text(length(structure)+1,-0.2,'structure','horizontalalign','left')
    
    protected_loop = intersect(intersect( find(rprof<0.125),strfind(structure,'.')),which_pos);
    exposed_pair   = intersect(intersect( find(rprof>0.5),find(structure~='.')),which_pos);
    plot_res = which_pos; clr = [0 0 0]; plot_symbols(rprof,rprof_err,plot_res, clr); 
    plot_res = protected_loop; clr = [1 0.5 0]; plot_symbols(rprof,rprof_err,plot_res, clr); 
    plot_res = exposed_pair;   clr = [1 0 0];  plot_symbols(rprof,rprof_err,plot_res, clr); 
    ylim([-1 3])
    hold off
    box off
    title( headers{i},'interpreter','none')
    xlabel('Position');

    goodpos = setdiff(setdiff(which_pos,exposed_pair),protected_loop);
    fprintf('Length of region probed: %d, Number of good positions: %d. Eterna classic score: %6.1f\n', length(which_pos),length(goodpos),100*length(goodpos)/length(which_pos));
    if length(good_idx)>1 & i ~= good_idx(end); pause;end;
end

%%%%%%%
function plot_symbols(rprof,rprof_err,plot_res,clr); 

plot(plot_res,rprof(plot_res),'o','color',clr,'markerfacecolor',clr);

if length(rprof_err) == 0; return; end;
for i = 1:length(plot_res)
    plot( plot_res(i)*[1 1], rprof(plot_res(i)) + rprof_err(plot_res(i))*[-1 1],'color',clr,'linew',1.5 );
end

