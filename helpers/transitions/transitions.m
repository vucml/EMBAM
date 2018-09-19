function trans_row = transitions(data_row, from_mask, to_mask, transit_func, ...
				 step, params)
%TRANSITIONS   Calculate transitions between elements of a vector.
%
%  trans_row = transitions(data_row, from_mask, to_mask, transit_func,
%                          step, params)  
%
%  INPUTS:
%      data_row:  numeric row vector.
%
%     from_mask:  logical array the same size as data_row. false at
%                 positions i where the transition from data_row(i) to
%                 data_row(i + step) should be excluded (NaNs in the
%                 output).
%
%       to_mask:  same as from_mask, but is false at positions i where
%                 transitions from data_row(i - step) to data_row(i)
%                 should be excluded.
%
%  transit_func:  handle to a function that computes the transition
%                 between the values x (data_row(i)) and y
%                 (data_row(i + step)) in data_row. Should be of the
%                 form:
%                  x_to_y_trans = transit_func(x, y, params)
%
%          step:  integer indicating the number of elements ahead to
%                 consider a transition.
%
%        params:  structure containing additional inputs to
%                 transit_func.
%
%  OUTPUTS:
%     trans_row:  transition values. trans_row(i) gives the transition
%                 between data_row(i) and data_row(i + step).
%  EXAMPLES:
%  Using masks:
%  >> data_row = [24 23 19 2 5];
%  >> all = true(size(data_row));
%  >> not_19 = [true true false true true];
%  >> transitions(data_row, not_19, all, @distance, 1, struct)
%      -1 -4 NaN 3
%  >> transitions(data_row, not_19, not_19, @distance, 1, struct)
%      -1 NaN NaN 3
%
%  Defining a transit_func:
%  >> is_positive_lag = @(x, y, params) (y > x)
%  >> transitions(data_row, all, all, is_positive_lag, 1, struct)
%      0 0 0 1
%
%  Step size:
%  >> transitions(data_row, all, all, @distance, 2, struct())
%      -5 -21 -14

% all the work is done by conditional_transitions...just call
% it with a dummy condition function
trans_row = conditional_transitions(data_row, from_mask, to_mask, ...
                                    transit_func, @allow_all_transits, ...
                                    step, params);
%endfunction

function possible_transit = allow_all_transits(from_pt, priors, transit, params)
% Helper to allow all transitions (i.e., the 'transparent condition')
  possible_transit = transit;
%endfunction

