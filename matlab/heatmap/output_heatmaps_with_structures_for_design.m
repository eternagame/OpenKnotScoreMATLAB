function output_heatmaps_with_structures_for_design( image_dir, good_idx, r_norm, structure_map, structure_sets, structure_tags, pkg_sort_idx, headers, sequences, ids, titles, authors, BLANK_OUT5, BLANK_OUT3, tags)
% output_heatmaps_with_structures_for_design( image_dir, good_idx, r_norm, structure_map, structure_sets, structure_tags, pkg_sort_idx, headers, sequences, ids, titles, authors, BLANK_OUT5, BLANK_OUT3, tags)
%
% Input files
%  image_dir = Directory to output images.
%  idx = [integer] index of design for which to show heatmap
%  r_norm = [Ndesign x Nres x Nconditions] Reactivity matrix, normalized.
%  structure_map = [Ndesign x Nres x Npackages] 0/1 map
%            of paired/unpaired for each predicted structure
%  structure_sets = [Npackages x Ndesign] cell of cell of strings of predicted structures
%  structure_tags = cell of string, name of each package
%  pkg_sort_idx = permutation of packages (e.g., [2, 4, 1, 3]) to order
%             packages. Give [] to sort by Eterna score to first r_norm dataset.
%              Give -1 to instead output Eterna score of each r_norm dataset.
%  headers = cell of Ndesign strings describing each design (titles for
%  plot). Provide [] to auto-fill from ids,titles,authors.
%  sequences = cell of sequences for Ndesigns
%  ids = [Ndesign integers] ids 
%  titles  = {cell of Ndesign strings} design titles.
%  authors = {cell of Ndesign strings} authors
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues 
%  tags = strings that provide experimental names for each of the
%        Nconditions in r_norm
%
% (C) R. Das, Stanford University, 2023 
% 
if ~exist(image_dir,'dir') mkdir( image_dir ); end;
if ~exist('tags','var') tags = {'SHAPE, no Mg2+','SHAPE, +Mg2+'}; end
set(figure(5),'position',[100 500 1200 150+12*length(structure_sets)]);
set(gcf,'color','white')
clf

if isempty(headers)
    for i = 1:length(ids)
        headers{i} = sprintf('%d\t%s\t%s',ids(i),titles{i},authors{i});
    end
end

for n = 1:length(good_idx)
    idx = good_idx(n);
    make_heatmap_with_structures_for_design( idx, r_norm, structure_map, structure_sets, structure_tags, pkg_sort_idx, headers, sequences, BLANK_OUT5, BLANK_OUT3, tags, 0 );
    if ~isempty(ids) & ~isnan(ids(idx))
        image_file = sprintf('%s/%d_%s_%s.png',image_dir,ids(idx),strrep(titles{idx},'/','-'),authors{idx});
    else
        image_file = sprintf('%s/%s.png',image_dir,cleanup(headers{idx}));
    end
    fprintf( 'Making %d of %d: %s ... \n ', n,length(good_idx),image_file );    
    print( image_file,'-dpng','-r300' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name = cleanup(name);
name = strip(name,'_');
name = strrep(name,',','');
name = strrep(name,'/','_');
name = strrep(name,' ','_');
name = strrep(name,'(','');
name = strrep(name,')','');
name = strrep(name,sprintf('\t'),' ');
name = strip(name);

