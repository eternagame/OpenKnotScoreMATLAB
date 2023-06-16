function [all_corr_coef, pkg_sort_idx] = get_corr_coeff( r, structure_map, good_idx, structure_tags, BLANK_OUT3, BLANK_OUT5, corr_type, num_show, clip);
% [all_corr_coef, pkg_sort_idx] = get_corr_coeff( r, structure_map, good_idx, structure_tags, BLANK_OUT3, BLANK_OUT5, corr_type, num_show, clip);
%
% Note that values are clipped so that negative values go to zero; 
%  high values can be clipped to maximum specified by the 'clip' variable.
%
% Inputs
%  r = [Ndesign x Nres x Ncond] Reactivity matrix. Assume last of Ncond has
%  SHAPE data.
%  structure_map = [Ndesign x Nres x Npackages] 0/1 map
%            of paired/unpaired for each predicted structure
%  good_idx = [list of integers] index of designs to use for correlation
%            coefficients
%  structure_tags = cell of string, name of each package
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues 
%  corr_type  = correlation coefficient type (default:'Pearson')
%  num_show   = number of top algorithms to show (default: all)
%  clip       = maximum value to clip at [Default no clipping]
%
% Outputs
%   all_corr_coef = [Npackages] correlation coefficients for each package
%   pkg_sort_idx = [Npackages] permutation of package indices that sort from best to
%                    worst by correlation coefficient
%
% (C) R. Das, HHMI/Stanford University 2023.

if ~exist( 'corr_type','var') | length(corr_type)==0; corr_type = 'Pearson'; end;
if ~exist( 'num_show','var') | num_show == 0; num_show = length(structure_tags); end;
if ~exist( 'clip','var') | clip == 0; clip=Inf; end;
idx_shape_condition = size(r, 3);

which_pos = (BLANK_OUT5+1):(size(r,2)-BLANK_OUT3);
r_subset = min(max(r(good_idx, which_pos, idx_shape_condition),0),clip);
for n = 1:size( structure_map, 3)    
    structure_map_subset = structure_map(good_idx,which_pos,n);
    gp = find (~isnan(r_subset(:)));
%    R = corrcoef( r_subset(gp), structure_map_subset(gp));
    all_corr_coef(n) = corr(r_subset(gp), structure_map_subset(gp),'Type',corr_type);
end

all_corr_coef( isnan( all_corr_coef) ) = -1;
[~,pkg_sort_idx] = sort(all_corr_coef,'descend');
fprintf('Using correlation type %s',corr_type);
if clip < Inf; fprintf('and maximum value (clip) of %8.3f',clip); end;
fprintf('\n');
for n = pkg_sort_idx(1:num_show)
    fprintf( '%7.4f %s\n',all_corr_coef(n),structure_tags{n});
end
