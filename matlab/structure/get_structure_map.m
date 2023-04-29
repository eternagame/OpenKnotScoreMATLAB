function structure_map = get_structure_map( structure_sets, idx );
% structure_map = get_structure_map( structure_sets );
% (C) R. Das, HHMI, Stanford, 2023

structure_map = [];

%% NEED TO SET structure_map sixze based on maximum structure length in structure_sets.
N = 0;
for i = 1:length( structure_sets );
    for j = 1:length(structure_sets{i})
        N = max(N,length(structure_sets{i}{j}));
    end
end

for count = 1:length( structure_sets );
    design_indices = 1:length(structure_sets{count});
    if exist('idx','var'); design_indices = idx; end;
    for i = design_indices % loop over designs
        structure_map(i,:,count) = zeros(1,N);
        structure_map(i,strfind(structure_sets{count}{i},'.'),count) = 1;
    end
end
