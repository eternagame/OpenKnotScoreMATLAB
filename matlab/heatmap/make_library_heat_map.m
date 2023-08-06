function make_library_heat_map( r_norm, good_idx, structure_map, headers, BLANK_OUT5, BLANK_OUT3, tags_conditions);
% make_library_heat_map( r_norm, good_idx, structure_map, headers, BLANK_OUT5, BLANK_OUT3, tags_conditions);
%
% Inputs
%  r_norm = [Ndesign x Nres x Nconditions] Reactivity matrix, normalized.
%             If Nconditions = 1, will make heatmap side by side with
%                   structure map
%             If Nconditions > 1, will tile the heatmaps
%  good_idx = [list of integers] index of designs for which to show heatmap
%  structure_map = [Ndesign x Nres] 0/1 map
%            of paired/unpaired for each predicted structure. Give [] to
%            not show it.
%  headers = cell of Ndesign strings describing each design (titles for
%  plot)
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues 
%
% (C) R. Das, HHMI/Stanford, 2023

Nconditions = size(r_norm,3);
N = size(r_norm,2);
Nplots = Nconditions;
if ~isempty(structure_map) Nplots = Nconditions+1; end;

set(figure(2),'Position',[100 100 max(50+Nconditions*150,800) 700]); clf;
set(gcf,'color','white')

plot_width = 0.65/Nplots;
for i = 1:Nconditions
    %subplot(1,Nplots,i)
    %set(gca,'Position',[0.5 0.1 0.2 0.8]);
    axes('Position',[0.3+(i-1)*plot_width, 0.1, plot_width, 0.8]);

    imagesc(r_norm(good_idx,:,i),[-2 2])
    colormap([0.7 0.7 0.7; redwhiteblue(-1,1)])
    set(gca, 'TickLabelInterpreter','none','fontweight','bold','ticklength',[0.01,0.25],'tickdir','out' );
    if i == 1
        for n = 1:length(good_idx); header = headers{good_idx(n)}; labels{n} = header(1:min(50,length(header)));end;
        set(gca,'ytick',[1:length(good_idx)],'yticklabel', labels);
    else
        set(gca,'ytick',[]);
    end
    if length(good_idx)< 100; make_lines_horizontal; end
    xlim([BLANK_OUT5 N-BLANK_OUT3])
    title( 'SHAPE data')
    if exist( 'tags_conditions', 'var');
        h = title( strsplit(tags_conditions{i},'_'),'interp','none' );
    end;
    xlabel('Position');
    box off
end

if ~isempty(structure_map)
    axes('Position',[0.3+Nconditions*plot_width 0.1 plot_width 0.8])
    structure_map(:,1:BLANK_OUT5) = 0;
    structure_map(:,end-BLANK_OUT3+1:end) = 0;
    imagesc(0.02*structure_map(good_idx,:),[-0.05 0.05])
    colormap([0.7 0.7 0.7; redwhiteblue(-0.1,0.1)])
    set(gca,'ytick',[1:length(good_idx)],'yticklabel', [],'fontweight','bold','ticklength',[0.01,0.25],'tickdir','out')
    set(gca,'yticklabel',[])
    xlim([BLANK_OUT5 N-BLANK_OUT3])
    xlabel('Position');
    %colorbar()
    if length(good_idx)< 100; make_lines_horizontal; end
    title( 'Predicted (unpaired)')
    box off
end