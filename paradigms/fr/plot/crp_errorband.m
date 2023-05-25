function [nph, pph] = crp_errorband(crp_ci, listlength)
% CRP_ERRORBAND
%
% Adds an alpha transparent error-band to a lag-CRP plot by drawing two
% polygons, one for the negative lags and one for the positive lags.
% If you add this to the figure before plotting the lag-CRP itself, it will
% be arranged behind the lag-CRP markers (graphically important if they
% are, e.g., white-filled shapes)
%
% crp_ci: I've been using the output of bootstrap_ci, this is expecting 2
% rows by (LL*2)-1 columns.
% listlength: code uses this to grab the errorbars surrounding lag 0.
%
% Note: Should be easy to customize to allow a different range of lags. 
% Returns the handles for the negative polygon and positive polygon so 
% some customization can be done in calling script.
%

% vertices of polygon errorband
neg_lags_upper = [-5 -4 -3 -2 -1];
neg_lags_lower = [-1 -2 -3 -4 -5];

for i=1:length(neg_lags_upper)
    neg_vertices_upper(i,1) = neg_lags_upper(i);
    neg_vertices_upper(i,2) = crp_ci(2,listlength+neg_lags_upper(i));
    neg_vertices_lower(i,1) = neg_lags_lower(i);
    neg_vertices_lower(i,2) = crp_ci(1,listlength+neg_lags_lower(i));
end
neg_vertices = [neg_vertices_upper; neg_vertices_lower];
neg_poly = polyshape(neg_vertices);
nph = plot(neg_poly);
set(nph,'EdgeAlpha',0,'FaceColor',[0.5 0.5 0.5]);
hold on

pos_lags_upper = [1 2 3 4 5];
pos_lags_lower = [5 4 3 2 1];

for i=1:length(pos_lags_upper)
    pos_vertices_upper(i,1) = pos_lags_upper(i);
    pos_vertices_upper(i,2) = crp_ci(2,listlength+pos_lags_upper(i));
    pos_vertices_lower(i,1) = pos_lags_lower(i);
    pos_vertices_lower(i,2) = crp_ci(1,listlength+pos_lags_lower(i));
end
pos_vertices = [pos_vertices_upper; pos_vertices_lower];
pos_poly = polyshape(pos_vertices);
pph = plot(pos_poly);
set(pph,'EdgeAlpha',0,'FaceColor',[0.5 0.5 0.5]);
hold on



