function [p_recalls] = p_rec_core(recalls_matrix, subjects, list_length, ...
			     rec_mask, pres_mask)
%P_REC_CORE   Recall probability.
%
% Unlike P_REC, this function does not have any error-checking. It
% assumes all of the inputs listed below are given and formatted correctly.
% 
% p_recalls = p_rec_core(recalls_matrix, subjects, list_length,
%                   rec_mask, pres_mask)
%
% INPUTS:
%  recalls_matrix:  a matrix whose elements are serial positions of recalled
%                   items.  The rows of this matrix should represent recalls
%                   made by a single subject on a single trial.
%
%        subjects:  a column vector which indexes the rows of recalls_matrix
%                   with a subject number (or other identifier).  That is, 
%                   the recall trials of subject S should be located in
%                   recalls_matrix(find(subjects==S), :)
%
%     list_length:  a scalar indicating the number of serial positions in the
%                   presented lists.  serial positions are assumed to run 
%                   from 1:list_length.
%
%        rec_mask:  a logical matrix of the same shape as 
%                   recalls_matrix, which is false at positions (i, j) where
%                   the value at recalls_matrix(i, j) should be excluded from
%                   the calculation of the probability of recall.
% 
%       pres_mask:  a logical matrix of the same shape as the
%                   item presentation matrix (i.e., num_trials x
%                   list_length).  This mask should be true at
%                   positions (t, sp) where an item in the condition
%                   of interest was presented at serial position sp on
%                   trial t; and false everywhere else. note that this does
%                   NOT take into account the excluded positions from which
%                   items are recalled.
%
%  OUTPUTS:
%        p_recall:  a vector of probablities.  Its rows are indexed by
%        subject.

p_recalls = apply_by_index(@prec_for_subj, ...
			   subjects, ...
			   1, ...
			   {recalls_matrix, rec_mask, pres_mask}, ...
			   list_length);

function subj_p_recall = prec_for_subj(recalls, rec_mask, pres_mask, ...
				       list_length)
  % helper for p_rec:
  % calculates the probability of recall for one subject's recall
  % trials.  Probability of recall is defined as:
  %  (# unmasked items recalled) / (# of unmasked presentations)
  
  unmasked_recalls = recalls(rec_mask);
  
  rec_counts = sum(collect(unmasked_recalls, 1:list_length));
  pres_counts = collect(pres_mask, true); 
  
  subj_p_recall = rec_counts / pres_counts;