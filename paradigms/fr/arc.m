function [arc_scores] = arc(catrecalls_matrix, subjects, mask)
%ARC  Adjusted Ratio of Clustering.
%
% ARC quantifies the degree of clustering in a recall sequence.
% ARC = (R - E(R)) / (maxR - E(R)), where R is the total number of
% observed category repetitions, maxR is the maximum number of possible
% category repetitions, and E(R) is the expected number of chance
% repetitions. Based on Roenker et al., 1971.
%
% arc_scores = arc(catrecall_matrix, subjects, mask);
%
% INPUTS:
% catrecalls_matrix:  a matrix whose elements are the category
%                     identities of recalled items.  The rows of
%                     this matrix should represent recalls made by
%                     a single subject on a single trial.  After
%                     all valid recalls, the terminal positions
%                     should be filled by negative values, which
%                     will be masked.
% 
%          subjects:  a column vector which indexes the rows of
%                     catrecalls_matrix with a subject number (or
%                     other identifier).  That is, the recall
%                     trials of subject S should be located in
%                     catrecalls_matrix(find(subjects==S), :) 
%
% OUTPUTS:
%        arc_scores:  a vector of scores.  Each score is the mean
%                     arc score for a particular subject.
%
% EXAMPLES:
% >> catrecalls_matrix = [0 1 1 1 2 1 0 0 NaN NaN NaN; ...
%                         2 2 2 1 1 0 1 NaN NaN NaN NaN];
% >> subjects = [1; 2]; cats = [0 1 2];
% >> arc_scores = arc(catrecall_matrix, subjects, cats);

% sanity checks
if ~exist('catrecalls_matrix', 'var')
  error('You must pass a recalls matrix.')
elseif ~exist('subjects', 'var')
  error('You must pass a subjects vector.')
elseif size(catrecalls_matrix, 1) ~= length(subjects)
  error('catrecalls matrix must have the same number of rows as subjects.')
elseif ~exist('mask', 'var')
  % create not nan mask if none was given
  mask = ~isnan(catrecalls_matrix);
end
if size(mask) ~= size(catrecalls_matrix)
  error('recalls_matrix and mask must have the same shape.')
end

arc_scores = apply_by_index(@arc_for_subj, ...
			    subjects, ...
			    1, ...
			    {catrecalls_matrix, mask});
%endfunction

function subj_arc_score = arc_for_subj(cat_recalls, mask)
  % Helper for arc:
  % calculates the adjusted ratio of clustering for one subject's
  % recall trials; returns a score.
  
  num_trials = size(cat_recalls,1);
  trial_arc_scores = NaN(1,num_trials);
  
  for i=1:num_trials
    
    cat_trial = cat_recalls(i,mask(i,:));
    cats = unique(cat_trial);
    num_cat = length(cats);
    
    % R is observed cat repetitions
    R = sum(diff(cat_trial)==0);
    % N is total # of items recalled
    N = length(cat_trial);
    % maxR is total possible cat repetitions
    maxR = N - num_cat;
    % E_R is expected or chance number of category repetitions
    n = zeros(1,num_cat);
    for j=1:num_cat
      n(j) = sum(cat_trial==cats(j));
    end    
    E_R = (sum(n.^2)/N) - 1;
    trial_arc_scores(i) = (R - E_R) / (maxR - E_R);
    %keyboard
  end

  subj_arc_score = nanmean(trial_arc_scores);
