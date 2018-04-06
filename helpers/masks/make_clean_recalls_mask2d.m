function mask = make_clean_recalls_mask2d(recalls)
%MAKE_CLEAN_RECALLS_MASK2D  Exclude repeats and intrusions.
%
%  Makes a mask of the same shape as recalls which is false at 
%  positions (i,j) if recalls(i,j) is an intrusion, 
%  repeat, or empty cell.
%
%  mask = make_clean_recalls_mask2d(recalls)
  
% sanity:
if ndims(recalls) ~= 2
  error('recalls must be two-dimensional.')
end

mask = false(size(recalls));

list_length = max(recalls(:));
for i = 1:list_length
  [row, col] = find(recalls == i);
  [first_row, first_ind] = unique(row, 'first');
  mask(sub2ind(size(recalls), first_row, col(first_ind))) = true;
end
