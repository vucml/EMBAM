function arr = catcell(c)
% CATCELL(C) Concatenates the row vectors in a cell array into a single
%            vector
%
% INPUTS:
%     c:  a cell array containing only row vectors or scalars
% 
% OUTPUTS:
%   arr:  a row vector containing the values found in c

  % for now, I'm not going to do type or sanity checking; this is a
  % one-liner designed to help re-obtain the old-style output of
  % conditional_transitions.
  arr = [c{:}];
%endfunction
