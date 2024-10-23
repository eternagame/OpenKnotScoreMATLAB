function [structure, helix_map] = convert_bps_to_structure3( bps, nres, ORDER_LAYERS_BY_HELIX_LENGTH)
% structure = convert_bps_to_structure3( bps, nres )
%
%  Takes in list of base pairs & structure length and 
%    outputs structure in dot-bracket format.
%
% INPUTS
%  bps  = [REQUIRED] list of base pairs (i,j) (Nbp X 2 matrix )
%  nres = [REQUIRED] total number of residues
%  ORDER_LAYERS_BY_HELIX_LENGTH = When there are pseudoknots prioritize use of (), then [], then {} 
%            for the set of helices with most pairs (Default 1). If set to
%            0, use (), [], based on first encounter.
%
%  (C) R. Das, Stanford University, 2009-2015, 2017
%  (C) R. Das, HHMI 2023-24
%
if ( nargin < 2 ) help( mfilename ); return; end;
if ~exist( "ORDER_LAYERS_BY_HELIX_LENGTH", 'var'); ORDER_LAYERS_BY_HELIX_LENGTH = 1; end;

if length( bps ) > 0 & size( bps, 2 ) ~= 2
    fprintf( 'bps must be Nbps x 2 array\n' );
    return;
end

% initialize structure with '.'
structure = '';
for i = 1:nres; structure = [structure,'.']; end;
left_delim  = '([{<abcdefghijklmnopqrstuvwxyz';
right_delim = ')]}>ABCDEFGHIJKLMNOPQRSTUVWXYZ';

if length( bps ) == 0; return; end;

% break up into helices (stems)
stems = parse_stems_from_bps( bps );

% convert into helix_map format (this is historical):
% Looks like partner i of first base pair, partner j of first base pair,
% length of helix
helix_map = [];
for i = 1:length(stems)
    helix_map = [helix_map; stems{i}(1,1), stems{i}(1,2), size( stems{i}, 1 ) ];
end


% order helices by length, longest to shortest
if ORDER_LAYERS_BY_HELIX_LENGTH
    [~,idx] = sort( helix_map(:,3) );
    idx = fliplr( idx' );
else
    % order helices by order of appearance in structure
    [helix_map_sort,idx] = sort(helix_map(:,1));
end

% "layers" are helices that can be connected by (), then by [], then by {},
% ...
helix_layers = {};
Ndelim = length(left_delim);
for n = 1:length( idx )
    i = idx( n );
    found_layer = 0;
    test_helix = helix_map(i,:);
    for j = 1:min(Ndelim,length(helix_layers))
        if ( not_pseudoknotted( helix_layers{j}, test_helix ) )
            helix_layers{j} = [ helix_layers{j}; test_helix ];
            found_layer = 1; break;
        end
    end
    if ( ~found_layer ) % new layer
        helix_layers = [helix_layers, helix_map(i,:) ];
    end
end

if ORDER_LAYERS_BY_HELIX_LENGTH
    % reorder so that helix layers with the most base pairs are first; 
    %  for ties, prioritize order of appearance in the structure.
    for j = 1:length( helix_layers )
        num_bps(j,:) = [sum( helix_layers{j}(:,3) ),  -helix_layers{j}(1,1) ];
    end
    [~,sorted_idx] = sortrows( num_bps, 'descend' );
    helix_layers = helix_layers( sorted_idx );
end

% ready to create output string in dot-parens notation!
for j = 1:length( helix_layers )
    helix_layer = helix_layers{j};
    for i = 1:size( helix_layer, 1 )
        structure = pk_bracket_substitute(structure, helix_layer(i,:), left_delim(j), right_delim(j) );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = not_pseudoknotted( helix_layer, test_helix )

ok = 1;
for i = 1:size( helix_layer, 1 )
    % "crossing" base pairs
    if ( helix_layer(i,1) < test_helix(1) & ...
         helix_layer(i,2) < test_helix(2) & ...
         helix_layer(i,2) > test_helix(1) ) ...
         | (helix_layer(i,1) > test_helix(1) & ...
            helix_layer(i,2) > test_helix(2) & ...
            test_helix(2) > helix_layer(i,1) )
        %fprintf( 'Crossing: existing helix %d-%d with test %d-%d\n',helix_layer(i,1),helix_layer(i,2), test_helix(1),test_helix(2) );  
        ok = 0;
        return;
    end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str_return = pk_bracket_substitute(str_input, helix_layer_sub, left_delim, right_delim)
str_return = [...
    str_input(1:(helix_layer_sub(1)-1)), ...
    repmat(left_delim,1,helix_layer_sub(3)), ...
    str_input((helix_layer_sub(1)+helix_layer_sub(3)):(helix_layer_sub(2)-helix_layer_sub(3))),...
    repmat(right_delim,1,helix_layer_sub(3)), ...
    str_input((helix_layer_sub(2)+1):end)];



    