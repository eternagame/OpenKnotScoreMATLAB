function [F1,structures_F1] = get_F1_scores(sequences, structure_sets, structure_tags, structure_tags_subset)

N = length(sequences);
idx = 1:N;
structures_F1 = cell(N,1);
F1 = zeros(N,1);

if ~exist( 'structure_tags_subset', 'var'); structure_tags_subset = structure_tags; end;
which_tags = [];
for n = 1:length(structure_tags_subset)
    m = find(strcmp(structure_tags,structure_tags_subset{n}));
    assert( ~isempty(m) );
    which_tags = [which_tags, m];
end

parfor k = idx
    structure_set = {};
    for n = 1:length(which_tags); structure_set{n} = structure_sets{which_tags(n)}{k}; end
    Nstruct = length(structure_set);
    F1_avg = zeros(1,Nstruct); 
    for m = 1:Nstruct
        ref_str = structure_set{m};
        F1_test = [];
        for n = setdiff(1:length(structure_set),m)
            test_str = structure_set{n};
            [sens, ppv] = calc_sens_ppv_bp(test_str, ref_str);;
            %fprintf('ref %s\ntarget %s\nSens %5.1f PPV %5.1f\n\n',ref_str,test_str,100*sens,100*ppv);
            F1_test(n) = 2/(1/sens + 1/ppv);
        end
        F1_avg(m) = mean(F1_test);
        %fprintf('%5.1f %s\n',100*F1_avg(m),ref_str);
    end
    [F1_max,max_idx] = max(F1_avg);
    F1(k) = F1_max;
    structures_F1{k} = structure_set{max_idx};
end
