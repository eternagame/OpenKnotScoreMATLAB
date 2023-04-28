function detect_tert_structure( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3, structures_rnastructure, headers, structures,  structures_shapeknots, r_norm_err  );
% detect_tert_structure( sequences, r_norm, good_idx, BLANK_OUT5, BLANK_OUT3, structures_rnastructure, headers, structures,  structures_shapeknots, r_norm_err  );
%
% (C) R. Das, HHMI/Stanford 2023

set(figure(15),'color','white','position',[100 500 1200 250]); clf;
N = length(sequences{1});
which_pos = [BLANK_OUT5+1:N-BLANK_OUT3];
rprof_err = [];
for i = good_idx'
    sequence =sequences{i};
    structure = structures{i};
    structure_rnastructure = structures_rnastructure{i};
    structure_pk = '';
    if i <= length(structures_shapeknots); structure_pk = structures_shapeknots{i}; end
    if length(structure_pk) == 0; structure_pk = repmat('.',length(sequence),1); end;
    rprof = r_norm(i,:);
    if exist('r_norm_err','var'); rprof_err = r_norm_err(i,:); end;
    plot( rprof, 'k','linew',2 ); hold on;
    plot([0 length(sequence)],[0 0],'k')
    plot([0 length(sequence)],0.25*[1 1],'k:')
    for m = 1:length(structure)
        text(m,-0.2,structure_rnastructure(m),'horizontalalign','center');
        text(m,-0.4,structure_pk(m),'horizontalalign','center');
        text(m,-0.6,structure(m),'horizontalalign','center');
        text(m,-0.8,sequence(m),'horizontalalign','center');
    end
    text(length(structure)+1,-0.2,'SHAPE-RNAstructure','horizontalalign','left')
    text(length(structure)+1,-0.4,'SHAPE-knots','horizontalalign','left');
    text(length(structure)+1,-0.6,'Eterna','horizontalalign','left');
    
    protected_loop = intersect(intersect( find(rprof<0.25),strfind(structure_rnastructure,'.')),which_pos);
    protected_loop2 = intersect(protected_loop,protected_loop+1);
    protected_loop2 = union(protected_loop2,protected_loop2-1);
    exposed_pair = intersect(setdiff( find(rprof>0.25),strfind(structure_rnastructure,'.')),which_pos);
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

%%%%%%%
function plot_symbols(rprof,rprof_err,plot_res,clr); 

plot(plot_res,rprof(plot_res),'o','color',clr,'markerfacecolor',clr);

if length(rprof_err) == 0; return; end;
for i = 1:length(plot_res)
    plot( plot_res(i)*[1 1], rprof(plot_res(i)) + rprof_err(plot_res(i))*[-1 1],'color',clr,'linew',1.5 );
end

