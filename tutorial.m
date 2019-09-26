% basic script to run all

%% recursively find the EMBAM directory
cwd = pwd;
[~,name,~] = fileparts(cwd);

while ~strcmp(name, 'EMBAM')
    cd ..
    cwd = pwd;
    [~,name,~] = fileparts(cwd);
end
addpath(genpath(cwd));

% or simply:
% init_embam();


%% load asymFR3 data in which three conditions: toronto noun pool lists, uncategorized lists, categorized lists
load data_asymFR3.mat

% get sublist of just toronto word pool
toronto = trial_subset(data.pres.listtype(:,1)==2,data);

% load semantic mat (870x870) where length(unique(data.pres_itemnos(:))) == 870
glove = load('asym_glove_raw.mat');

%% dist fact
dist_facts = dist_fact(toronto.rec_itemnos, toronto.pres_itemnos, toronto.subject, glove.sem_mat);


%% run lbc score % lbc won't apply for AsymFR3 dataset
lbc_score = lbc(data.pres.listtype, data.rec.listtype, data.subject); 


%% run temporal factor score
tf_toronto = temp_fact(toronto.recalls, toronto.subject, toronto.listLength);


%% serial position curve
% create a matrix of [subjects X probability recall by serial pos]
spc_toronto = spc(toronto.recalls, toronto.subject, toronto.listLength);
plot_spc(spc_toronto);
title('serial position curve');


%% probability of first recall
% createa matrix of [subject X probability recall by serial pos]
% pfr calls pfr_core, which depends on spc_core
pfr1 = pfr(toronto.recalls, toronto.subject, toronto.listLength);
plot_spc(pfr1);
title('probability of first recall');


%% lag-CRP
lag_crps = crp(toronto.recalls, toronto.subject, toronto.listLength);
plot_crp(lag_crps);
title('basic crp');


%% semantic crp
[bin_mean, sem_freq_binned_all, slope_subj, slope] = sem_crp(toronto, glove.sem_mat,5);


%% arc  NEED TO CHECK
arc_scores = arc(data.rec.listtype, data.subject);
