function mask = make_mask_only_final2d(recalls_matrix)
%  mask = make_mask_only_final2d(recalls_matrix)
%  Makes a mask of the same shape as recalls_matrix, which is true at 
%  the final recalled items
%

% sanity:
if ~exist('recalls_matrix', 'var')
  error('You must pass an recalls matrix.')
elseif ndims(recalls_matrix) ~= 2
  error('recalls_matrix must be two-dimensional.')
end

rows = size(recalls_matrix, 1);
mask = ~make_blank_mask(recalls_matrix);

for i = 1:rows
  mask(i, :) = make_mask_only_final1d(recalls_matrix(i, :));
end
%endfunction
