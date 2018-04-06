function mask_rec = update_recalls_mask(recalls, mask_rec, mask_pres)
% UPDATE_RECALLS_MASK
%
% function mask_rec = update_recalls_mask(recalls, mask_rec, mask_pres)



if ~exist('mask_rec','var')
  error('You must provide a recalls mask.');
elseif ~exist('mask_pres','var')
  error('You must provide a presentation mask.');
elseif ~exist('recalls','var')
  error('You must provide a recalls matrix.');
elseif size(mask_rec,1) ~= size(mask_pres,1)
  error('The two masks must have the same number of rows.')
elseif size(recalls,1) ~= size(mask_pres,1)
  error(['The recalls matrix and the mask must have the same number' ...
	 ' of rows.']);
end

for i=1:size(mask_rec,1)
  mask_these = find(mask_pres(i,:) == 0);
  for j=1:length(mask_these)
    mask_rec(i,recalls(i,:)==mask_these(j)) = 0;
  end
end





