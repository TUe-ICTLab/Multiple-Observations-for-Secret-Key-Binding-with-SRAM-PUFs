%% example call for sim_MOhds
clear all
close all

rng('shuffle')             

% settings
n_sim = 10^3;
R = 3; % 1/R
t_max = 5; % max enrollments

for i = 1:5 
    % Divide in parts to store inbetween results or to divide over 
    % calculations over multiple workers. !Don't forget to initialize the
    % random generator for each worker!
    fprintf('iteration %d, started: %s\n',i,datetime)
    [n_errors,sim_settings] = sim_MOhds(n_sim,1/R,t_max);
    save(sprintf('sim_results/MO_R%d_E%d_t%d_w1part%d.mat',R,log10(n_sim),t_max,i));
end



%% example call for sim_SMOhds
clear all
close all

rng('shuffle')             

% settings
n_sim = 10^3;
R = 3; % 1/R
t_min = 2; % t0 enrollments
t_max = 7; % t_max = t0+ max reconstructions

for i = 1:5
    % Divide in parts to store inbetween results or to divide over 
    % calculations over multiple workers. !Don't forget to initialize the
    % random generator for each worker!
    fprintf('iteration %d, started: %s\n',i,datetime)
    [n_errors,hist_n_fails,sim_settings] = sim_SMOhds(n_sim,1/R,t_min,t_max);
    save(sprintf('sim_results/SMO_R%d_E%d_t%d_t%d_w1part%d.mat',R,log10(n_sim),t_min,t_max,i));
end

