function [transits_array] = possible_transitions(serial_position, ...
						 prior_recalls, ...
                                                 ~, params)
%POSSIBLE_TRANSITIONS   Possible transitions, given previous recalls.
%
%  Returns the possible lags from a given serial position, excluding
%  lags to serial positions which have already been recalled.
%
%  NOTE: this function is meant to be passed as a condition function to
%  conditional_transitions(); its arguments are dictated by the
%  requirements of that function.
%
%  [transits_array] = possible_transitions(serial_position, ...
%      prior_recalls, transition, params)
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
%                    not used).
%
%           params:  a structure with fields:
%                     'from_mask_pres' - true for positions that are
%                                        valid to transition from.
%                     'to_mask_pres'   - true for positions that are
%                                        valid to transition to.
%
%  OUTPUTS:
%   transits_array:  a row vector of possible transitions from the
%                    current serial position (excluding transitions to
%                    previously-recalled serial positions)
%
%  EXAMPLES:
%  >> sp = 4;
%  >> prior_recalls = [3 2];
%  >> params.list_length = 6;
%  
%  % with no prior recalls, all transitions for list_length 6 are possible:
%  >> possible_transitions(sp, [], -1, params)
%  ans =
%     -3   -2   -1   1   2
%     
%  % with two prior recalls, transitions to those positions are excluded:
%  >> possible_transitions(sp, prior_recalls, -1, params)
%  ans = 
%     -3    1    2

to_mask_pres = params.to_mask_pres;
from_mask_pres = params.from_mask_pres;

% transitions from an intrusion, a repeat, or a masked out position
% are invalid, so there are no possible transitions from them
if serial_position < 1 || ...
   any(serial_position == prior_recalls) || ...
   ~from_mask_pres(serial_position)
  transits_array = [];
  return
end

% exclude the "from" item
to_mask_pres(serial_position) = false;

% exclude prior recalls
to_mask_pres(prior_recalls(prior_recalls > 0)) = false;

% calculate all possible lags
transits_array = find(to_mask_pres) - serial_position;

