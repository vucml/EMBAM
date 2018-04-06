function [trans_row, cond_cell] = conditional_transitions(data_row, ...
          from_mask, to_mask, transit_func, condition, step, params)
%CONDITIONAL_TRANSITIONS   Calculate transitions that meet conditions.
%
%  This function behaves similarly to transitions, with two important
%  differences: first, that the inclusion in the output of a particular
%  transition between two elements is conditional upon that transition
%  being in a vector returned by a condition function; and second, that
%  two values (one of actual transitions, one of conditions) are
%  returned.  Effectively, this function allows you to dynamically
%  determine which transitions to count, as opposed to the static
%  approach encouraged by transitions() (i.e., using a mask to exclude
%  certain positions).
%
%  A particular transition between two items, as computed by
%  transit_func, will be included in the output array of transitions
%  IF AND ONLY IF its value is found in the array returned by the
%  condition function.  Otherwise, a NaN will be placed transitions
%  array.
%
%  [trans_row, cond_cell] = conditional_transitions(data_row, from_mask,
%                           to_mask, transit_func, condition, step, params)
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
%     condition:  handle to a function which returns an array of
%                 possible transitions. Must have the form:
%
%                 poss_trans = condition(current, prev, trans, params)
%
%                 where current is the current element of data_row,
%                 prev is an array of the previous elements, trans is
%                 the current transition returned by transit_func, and
%                 params is a structure specifying options. If no
%                 transitions are possible, condition should output an
%                 empty array.
%
%          step:  integer indicating the number of elements ahead to
%                 consider a transition.
%
%        params:  structure specifying options for transit_func and/or
%                 condition.
%
%  OUTPUTS:
%     trans_row:  transition values. trans_row(i) gives the transition
%                 between data_row(i) and data_row(i + step).
%
%     cond_cell:  cell array of the possible transitions returned by the
%                 condition function. The corresponding possible
%                 transitions for trans_row(i) are stored in
%                 cond_cell{i}.
%
%  EXAMPLES:
%  Use a condition function that includes only small lags:
%  >> data_row = [22 19 23 4 5 6];
%  >> exclude_big_lags = @(current, prev, trans, params) ([-3:3]);
%  >> mask = true(size(data_row));
%  >> [trans, conds] = conditional_transitions(data_row, mask, mask, ...
%                      @distance, exclude_big_lags, 1, struct);
%  >> trans
%      -3 NaN NaN 1 1
%  >> celldisp(conds)
%  conds{1} =
%      -3 -2 -1 0 1 2 3
%  ...

% sanity checks
if ~islogical(from_mask) || ~islogical(to_mask)
  error('Masks must be logical variables');
end
if step < 1
  error('Non-positive steps are not supported');
end
if ~isa(condition, 'function_handle')
  error('condition must be a function handle');
end

% initialization: trans_row should have as many elements as
% data_row - step, as should cond_cell, since that is the maximum
% number of transitions and conditions we can calculate
row_length = size(data_row, 2);  
trans_row = NaN(1, row_length - step);
cond_cell = cell(1, row_length - step);

for i = 1:row_length - step
  if from_mask(i) && to_mask(i + step)
    % this transition is not masked out
    from_pt = data_row(i);
    to_pt = data_row(i + step);

    % calculate the current transition and condition and append
    % them to trans_row and cond_cell 
    transition = transit_func(from_pt, to_pt, params);
    poss_trans = condition(from_pt, data_row(1:i-1), transition, params);

    if ~isempty(transition) && any(transition == poss_trans)
      % the transition meets the condition; include it in trans_row
      trans_row(i) = transition;
    end
    
    % whether or not the transition meets the condition, we update
    % cond_cell and priors
    cond_cell{i} = poss_trans;
  end
end

%endfunction
  