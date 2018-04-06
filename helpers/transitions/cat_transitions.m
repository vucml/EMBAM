function [transits_array] = cat_transitions(serial_position, ...
				            prior_recalls, ~, ...
				            params)
%CAT_TRANSITIONS   Get possible transitions, conditional on category.
%
%  Returns the possible lags from a given serial position, excluding
%  lags to serial positions which have already been recalled,
%  conditional on those lags being of the same category (or a different
%  category) from the just-recalled item. Lag only includes items
%  that meet the criteria (within- or between-category), rather
%  than usual lag that counts all serial positions.
%
%  NOTE: this function is meant to be passed as a condition function to
%  conditional_transitions(); its arguments are dictated by the
%  requirements of that function.
%
%  [transits_array] = cat_transitions(serial_position, prior_recalls, ...
%                                          transition, params)
%
%  INPUTS:
%  serial_position:  the serial position from which possible transitions
%                    should be calculated.  If this value is less than 1
%                    (i.e., the 'recall' was an intrusion or empty
%                    cell), an empty array is returned.
%
%    prior_recalls:  a row vector of serial positions which have already
%                    been recalled; transitions to these serial
%                    positions are excluded from the output.
%
%       transition:  the current transition value (accepted here to meet
%                    the requirements of conditional_transitions(), but
%                    not used)
%
%           params:  structure with the following fields:
%                     cat_type       - 1 for within-category
%                                      transitions, 2 for between-
%                                      category transitions.
%                     pres_catlabels - numeric category label for each
%                                      item.
%                     to_mask_pres   - array that is true for valid
%                                      serial positions to transition
%                                      to.
%
%  OUTPUTS:
%   transits_array:  a row vector of possible transitions from the
%                    current serial position (excluding transitions to
%                    previously-recalled serial positions).
%                    transitions are returned as positions relative
%                    to the originating item, ignoring invalid
%                    intervening items.
%
%  EXAMPLES:
%  >> sp = 4;
%  >> prior_recalls = [3 2];
%  >> params.cat_type = 1;
%  >> params.pres_catlabels = [1 0 1 1 0 1];
%  >> params.to_mask_pres = true(size(params.pres_catlabels));
%  
%  % with no prior recalls 
%  >> cat_transitions(sp, [], -1, params)
%  ans =
%     -2    -1    1
%     
%  % with two prior recalls, transitions to those positions are excluded:
%  >> cat_transitions(sp, prior_recalls, -1, params)
%  ans = 
%     -2     1

pres_catlabels = params.pres_catlabels;
cat_type = params.cat_type;
from_mask_pres = params.from_mask_pres;
to_mask_pres = params.to_mask_pres;

% transitions from an intrusion, empty cell, or repeat are never
% allowed
if serial_position < 1 || any(serial_position == prior_recalls) || ...
       ~from_mask_pres(serial_position)
  transits_array = [];
  return
end

% get a mask for included positions, based on category. The "from"
% item is always included
if cat_type == 1
  cat_mask = pres_catlabels == pres_catlabels(serial_position);
else
  cat_mask = pres_catlabels ~= pres_catlabels(serial_position);
  cat_mask(serial_position) = true;
end

% category position, for each item
cat_pos_map = cumsum(cat_mask);

% include items fitting the category transition
total_mask = cat_mask & to_mask_pres;

% exclude the "from" item
total_mask(serial_position) = false;

% exclude prior recalls
total_mask(prior_recalls(prior_recalls > 0)) = false;

% calculate all possible lags
transits_array = cat_pos_map(total_mask) - cat_pos_map(serial_position);

