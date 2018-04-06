function cl = cat_lag(from_pt, to_pt, params)
%CAT_LAG   Lag between two items, conditional on their category.
%
%  cl = cat_lag(from_pt, to_pt, params)
%
%  INPUTS:
%  from_pt:  serial position of the item being transitioned from.
%
%    to_pt:  serial position of the item being transitioned to.
%
%   params:  struct with fields:
%             pres_catlabels - vector giving numeric labels for category
%                              by serial position.
%             cat_type       - 1 for within-category transitions;
%                              2 for between-category transitions.
%
%  OUTPUTS:
%       cl:  category lag. When calculating lag, only the included
%            transitions (e.g. within-category transitions) are counted;
%            the other items are clipped out before calculating lag.

% find included items (the from_pt is always included)
if params.cat_type == 1
  cat_mask = params.pres_catlabels == params.pres_catlabels(from_pt);
else
  cat_mask = params.pres_catlabels ~= params.pres_catlabels(from_pt);
  cat_mask(from_pt) = true;
end

if cat_mask(to_pt)
  % get category position (serial position with the invalid items
  % clipped out), get lag based on that
  cat_pos_map = cumsum(cat_mask);
  cl = cat_pos_map(to_pt) - cat_pos_map(from_pt);
else
  % the to_pt isn't included, so lag is undefined
  cl = [];
end

