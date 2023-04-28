function mfe_structure_map = get_mfe_structure_map( mfe_structures, idx );
% mfe_structure_map = get_mfe_structure_map( mfe_structures );
% (C) R. Das, HHMI, Stanford, 2023

mfe_structure_map = [];
for count = 1:length( mfe_structures );
    design_indices = 1:length(mfe_structures{count});
    if exist('idx','var'); design_indices = idx; end;
    for i = design_indices % loop over designs
        mfe_structure_map(i,:,count) = zeros(1,length(mfe_structures{1}{1}));
        if length(mfe_structures{count}{i})> 1000
            count
            i
        end
        mfe_structure_map(i,strfind(mfe_structures{count}{i},'.'),count) = 1;
    end
end
