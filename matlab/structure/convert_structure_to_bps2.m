function bps = convert_structure_to_bps2( structure );
%  bps = convert_structure_to_bps2( structure );
%
%  Updated to handle weird structures somewhat gracefully...
%
%  INPUTS
%  structure = structure in dot-parens notation, using characters .()[]{}
%              can also specify helices by other characters, like ...aaa...aaa....
%  
%  OUTPUT
%  bps  = matrix of base pairs.
% (C) Rhiju Das, Stanford University, 2011-2012, 2017. 
% (C) Rhiju Das, HHMI, 2023.
%
%

bps = [];
bps = get_bps( structure, bps, '(', ')' );
bps = get_bps( structure, bps, '[', ']' );
bps = get_bps( structure, bps, '{', '}' );
bps = get_bps( structure, bps, '<', '>' );

% allow any other chars to represent additional base pairs
other_chars = setdiff( structure, '()[]{}<>.' );
for n = 1:length( other_chars )
    pos = find( structure == other_chars(n) );
    if mod(length(pos),2) ~= 0
        fprintf( 'Non-even number of character %s in structure!? %s\n',other_chars(n),structure )
        continue;
    end
    assert( mod(length(pos),2) == 0 );
    N = length(pos)/2;
    i = pos(1 : N); % first N
    j = pos(2*N: -1: (N+1) ); % last N are assumed to be pairs
    if ( N > 1 ) % each of the two sets of characters must be contiguous
        if ~all( i(1:end-1)+1 == i(2:end) ) | ~all( j(1:end-1)-1 == j(2:end) )
            fprintf( 'Non-contiguous stem from character %s in structure! %s\n',other_chars(n),structure )
        end
    end
    bps = [bps; i' j' ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bps = get_bps( structure, bps, left_delim, right_delim );

LEFT_BRACKET = [];
for i = 1:length(structure )
  switch structure(i)
   case left_delim
    LEFT_BRACKET = [LEFT_BRACKET, i];
   case right_delim
    if length(LEFT_BRACKET)==0
        fprintf('Error in matching %s in structure %s\n',left_delim, structure);
        return;
    end
    bps = [bps; LEFT_BRACKET(end), i];
    LEFT_BRACKET = LEFT_BRACKET( 1 : end-1 );
  end
end

