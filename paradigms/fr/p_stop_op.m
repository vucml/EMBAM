function [p_stops,denoms,r_index] = p_stop_op(recalls,time_mat,...
                              rec_length,exit_time_thresh,...
                              subjects,mask)
%P_STOP_OP   Probability of stopping recall by output position.
%
% [p_stops,denoms] = p_stop_op(recalls,time_mat,rec_length,exit_time_thresh,...
%                              subjects,mask)
%
% INPUTS:
%           recalls: a matrix whose elements are serial positions of recalled
%                    items.  The rows of this matrix should represent recalls
%                    made by a single subject on a single trial.
%
%          time_mat: a matrix whose elements are the millisecond times of the 
%                    recalled items.  The rows of this matrix should represent 
%                    one recall sequence made by a subject on a single trial.
%
%        rec_length: a column vector of the length, in milliseconds, of the 
%                    recall period for each trial.
%
%  exit_time_thresh: a scalar, in ms, representing time required between 
%                    the final recall in a trial and the end of the recall 
%                    period.  If the final recall occurred less than 
%                    exit_time_thresh away from the end of the recall period, 
%                    the trial will not be used in the analysis (the idea 
%                    being that perhaps the subject simply ran out of time 
%                    but was not done recalling).  A trial will only be 
%                    included if the final recall is greater than 
%                    exit_time_thresh away from the end of the recall period
%                    AND the time between the final recall and the end of
%                    of the recall period is greater than all of the
%                    inter-response times on that trial.  NOTE: If you do not
%                    want to exclude any trials based on these criteria, do not
%                    pass values for time_mat, rec_length, and exit_time_thresh.
%
%          subjects: a column vector which indexes the rows of recalls
%                    with a subject number (or other identifier).  That is, 
%                    the recall trials of subject S should be located in
%                    recalls(find(subjects==S), :)
%
%              mask: a logical matrix the same shape as recalls.  The mask
%                    should be true for any item in the condition of interest.
%                    If NOT given, a clean recalls mask is used (i.e., only
%                    correct recalls will be analyzed.)
%
% OUTPUTS:
%           p_stops: a matrix of stopping probabilities, with columns
%                    representing output positions, and rows representing
%                    subjects.
%
%            denoms: a matrix a denominator values that went into the 
%                    probability calculations.
%
% Notes about the mask: You can use the mask input to analyze different types
%                       of recalls (e.g., correct recalls, repetitions,
%                       intrusions).  There are masking functions that will
%                       create these.
%
%                       To only look at correct responses:
%                       mask = make_clean_recalls_mask2d(recalls)
%
%                       To only look at repetitions:
%                       mask = make_mask_only_reps2d(recalls);
%
%                       To look at intrusions, you can create a mask using an
%                       intrusions matrix (which must be the same size as 
%                       recalls, where a positive integer indicates a prior
%                       list intrusions, and a -1 indicates an extra list
%                       intrusion)
%
%                       To only look at PLIs:
%                       mask = make_mask_only_pli2d(intrusions);
%
%                       To only look at XLIs:
%                       mask = make_mask_only_xli2d(intrusions);
%
%                       You can use the repetition mask, PLI mask, and XLI mask
%                       to create a mask for all incorrect recalls.
%
% EXAMPLE:
% [p_stops,denom] = p_stop_op(recalls,time_mat,rec_length,12000,subjects,mask)

% sanity checks
if ~exist('recalls','var')
    error('You must pass a recalls matrix.')
elseif any([isempty(rec_length) isempty(exit_time_thresh) isempty(time_mat)]) ...
    & ~all([isempty(rec_length) isempty(exit_time_thresh) isempty(time_mat)])
    error(['You must pass a time_mat, recall_length scalar, and an ' ...
           'exit_time_thresh scalar, or all must be empty.'])
elseif ~exist('subjects','var')
    error('You must pass a subjects vector.')
elseif ~exist('mask','var')
    disp('No mask input. Using make_clean_recalls_mask2d as default mask')
    mask = make_clean_recalls_mask2d(recalls);
end

% make a mask for the final recalled item in each trial
mask_final_rec = make_mask_only_final2d(recalls);

% if desired, only include trials where the final recall occurred greater than
% exit_time_thresh away from the end of the recall period AND was longer than
% the longest IRT
if all([~isempty(rec_length) ~isempty(exit_time_thresh) ~isempty(time_mat)])

   % get the irts
   irts = diff(time_mat,1,2);
   irts(irts<=0) = NaN;

   % calculate the exit times (the time between the last recall and the end of
   % the recall period)
   time_mat_trans = time_mat';
   mask_final_rec_trans = mask_final_rec';
   exit_times = rec_length - time_mat_trans(mask_final_rec_trans);

   % find where the exit times are greater than exit_time_thresh AND longer
   % than the longest IRT
   good_trials_exit = exit_times > exit_time_thresh;
   good_trials_irt = exit_times > max(irts,[],2);

   % if there was only one recall in a trial, the above lines will not catch it as,
   % there are no IRTs.  First, find trials with only one recorded time.  Then
   % perform the exit time check.
   one_rec_trials = sum(logical(time_mat),2)==1;
   good_trials_one_rec = good_trials_exit & one_rec_trials;
   
   % combine good_trials_one_rec with good_trials_exit and good_trials_irt
   good_trials = (good_trials_exit & good_trials_irt) | good_trials_one_rec;

   % update masks to only include good trials
   mask_final_rec(~good_trials,:) = 0;
   mask(~good_trials,:) = 0;
end

[p_stops_denom,r_index] = apply_by_index(@p_stop_for_subj,...
                         subjects,...
                         1,...
                         {mask,mask_final_rec});
p_stops = p_stops_denom(:,:,1);
denoms = p_stops_denom(:,:,2);
%endfunction

function subj_p_stop = p_stop_for_subj(mask,mask_final_rec)
%
% does that actual probability of stopping calculation, which is the 
% number of times recall stopped over the number of times recall could have
% stopped, calculated for each output position
%

% To calculate the probability of stopping, we need to know the number of
% responses for each output position.  This is simply the sum of the mask over
% the rows.  For example, if you have a clean recalls mask, the sum will be
% the number of correct responses for each output position.
denom = sum(mask,1);

% We also need to know where recall stopped and whether the final recall was a
% member of the masked recalls.  This is the sum over rows of where both mask
% and mask_final_rec are true.
numerator = sum(mask & mask_final_rec,1);

% and divide to get probability of stopping for each output position
subj_p_stop = numerator./denom;
denom(denom==0) = NaN;
subj_p_stop(:,:,2) = denom;
%endfunction