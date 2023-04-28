function [all_corr_coef, pkg_sort_idx] = get_corr_coeff( r, mfe_structure_map, good_idx, mfe_tags, BLANK_OUT3, BLANK_OUT5);
% [all_corr_coef, pkg_sort_idx] = get_corr_coeff( r, mfe_structure_map, good_idx, mfe_tags, BLANK_OUT3, BLANK_OUT5);
%
% Inputs
%  r = [Ndesign x Nres x Ncond] Reactivity matrix. Assume last of Ncond has
%  SHAPE data.
%  mfe_structure_map = [Ndesign x Nres x Npackages] 0/1 map
%            of paired/unpaired for each predicted structure
%  good_idx = [list of integers] index of designs to use for correlation
%            coefficients
%  mfe_tags = cell of string, name of each package
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues 
%
% Outputs
%   all_corr_coef = [Npackages] correlation coefficients for each package
%   pkg_sort_idx = [Npackages] permutation of package indices that sort from best to
%                    worst by correlation coefficient
%
% (C) R. Das, HHMI/Stanford University 2023.

idx_shape_condition = size(r, 3);

which_pos = (BLANK_OUT5+1):(size(r,2)-BLANK_OUT3);
r_subset = max(r(good_idx, which_pos, idx_shape_condition),0);
for n = 1:size( mfe_structure_map, 3)    
    mfe_structure_map_subset = mfe_structure_map(good_idx,which_pos,n);
    gp = find (~isnan(r_subset(:)));
    R = corrcoef( r_subset(gp), mfe_structure_map_subset(gp));
    all_corr_coef(n) = R(1,2);
end

all_corr_coef( isnan( all_corr_coef) ) = -1;
[~,pkg_sort_idx] = sort(all_corr_coef);
for n = pkg_sort_idx(end:-1:1)
    fprintf( '%7.4f %s\n',all_corr_coef(n),mfe_tags{n});
end
