function output_structures_csv(csv_file_name,structure_tags,structure_sets,sequences);
%
% output_structures_csv(csv_file_name,structure_tags,structure_sets,sequences);
%
%
% (C) R. Das, Stanford University and HHMI, 2023

if ischar( structure_tags ); structure_tags = {structure_tags}; end;
if iscell( structure_sets ) & length(structure_sets) > 0 & ischar( structure_sets{1} ); structure_sets = {structure_sets}; end;

assert( length( structure_tags) == length(structure_sets));
N_structure_sets = length( structure_sets );
assert( N_structure_sets > 0 );

fid = fopen( csv_file_name, 'w');

Nstructures = length( structure_sets{1} );
for i = 1:N_structure_sets;
    assert( length(structure_sets{i}) == Nstructures);
end

Nres = length(sequences{1});
if length(sequences) > Nstructures
    for j = 1:N_structure_sets;
        for i = (Nstructures+1):length(sequences)
            structure_sets{j}{i}='';
        end
    end
    Nstructures = length(sequences);
end

% replace blank structures with '.....' or downstream structure
% processing can go awry.
for j = 1:N_structure_sets;
    for i = 1:length(sequences)
        if length(structure_sets{j}{i}) == 0
            structure_sets{j}{i}=repmat('.',1,Nres) ;
        end
    end
end


% header
fprintf(fid,'sequence,%s\n',strjoin(structure_tags,','));

% rows
for j = 1:Nstructures
    fprintf(fid,'%s',sequences{j});
    for i = 1:N_structure_sets;
        fprintf(fid,',%s',structure_sets{i}{j});
    end
    fprintf(fid,'\n',sequences{j});
end
fclose(fid);
fprintf('Created %s with %d sets of %d structures.\n',csv_file_name,N_structure_sets, Nstructures);





