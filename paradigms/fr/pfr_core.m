function [p_recalls] = pfr_core(recalls_matrix, subjects, list_length, ...
			   rec_mask, pres_mask)
%PFR_CORE   Probability of first recall.
%
%  Computes probability of recall by serial position for the
%  first output position.  Wraps spc.m. Unlike pfr.m, this function does 
%  not have any error-checking. It assumes all of the inputs listed below 
%  are given and formatted correctly.
%  
%  p_recall = pfr_core(recalls_matrix, subjects, list_length, rec_mask, pres_mask)
%
%  INPUTS: 
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
%                   trial t; and false everywhere else.
% 
%  OUTPUTS:
%        p_recall:  a matrix of probablities.  Its columns are indexed by
%                   serial position and its rows are indexed by subject.
%


% mask out everything but the first output position.
rec_mask(:,2:end) = 0;

% spc does the work here
p_recalls = spc_core(recalls_matrix, subjects, list_length, rec_mask, pres_mask);