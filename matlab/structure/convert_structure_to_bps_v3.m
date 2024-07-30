function bps = convert_structure_to_bps_v3( structure );
%  bps = convert_structure_to_bps3( structure );
%
%  Updated to handle weird structures somewhat gracefully...
%  Updated to handle a/A, b/B as delimiters.
%
%  INPUTS
%  structure = structure in dot-parens notation, using characters .()[]{}
%              can also specify helices by other characters, like ...aaa...aaa....
%  
%  OUTPUT
%  bps  = matrix of base pairs.
% (C) Rhiju Das, Stanford University, 2011-2012, 2017. 
% (C) Rhiju Das, HHMI, 2023-24.
%
%

bps = [];
bps = get_bps( structure, bps, '(', ')' );
bps = get_bps( structure, bps, '[', ']' );
bps = get_bps( structure, bps, '{', '}' );
bps = get_bps( structure, bps, '<', '>' );
letters = 'abcdefghijklmnopqrstuvwxyz';
for i = 1:length(letters)
    bps = get_bps( structure, bps, letters(i), upper(letters(i)) );
end

% allow any other chars to represent additional base pairs
other_chars = setdiff( structure, ['()[]{}<>.-:_,~',letters,upper(letters)] );
for n = 1:length( other_chars )
    pos = find( structure == other_chars(n) );
    if mod(length(pos),2) ~= 0
        warning( sprintf('Non-even number of character %s in structure!? %s\n',other_chars(n),structure ))
        continue;
    end
    assert( mod(length(pos),2) == 0 );
    N = length(pos)/2;
    i = pos(1 : N); % first N
    j = pos(2*N: -1: (N+1) ); % last N are assumed to be pairs
    bps = [bps; i' j' ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bps,errorcode] = get_bps( structure, bps, left_delim, right_delim );

LEFT_BRACKET = [];
errorcode = 0;
bps_input = bps; % save in case of error
for i = 1:length(structure )
  switch structure(i)
   case left_delim
    LEFT_BRACKET = [LEFT_BRACKET, i];
   case right_delim
    if length(LEFT_BRACKET)==0
        errorcode = 1;
        if contains('abcdefghijklmnopqrstuvwxyz',left_delim) & upper(left_delim)==right_delim
            warning(sprintf('Trying switched order of %s and %s.',left_delim,right_delim));
            [bps,errorcode] =  get_bps( structure, bps_input, right_delim, left_delim );
        end
        if errorcode 
            fprintf('Error in matching %s at position %d in structure %s \n',left_delim, i, structure );
        end
        return;
    end
    bps = [bps; LEFT_BRACKET(end), i];
    LEFT_BRACKET = LEFT_BRACKET( 1 : end-1 );
  end
end

