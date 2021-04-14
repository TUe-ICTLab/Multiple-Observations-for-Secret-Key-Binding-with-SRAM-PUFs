%% example plot results for sim_MOhds
clear all
close all 


% settings
n_sim = 10^6;
R = 6; % 1/R
t_max = 7; % max enrollments
storetable = true; % store results in .txt file for plotting with pgfplots

n_errors_tot = 0; % total nr. errors
n_sim_tot = 0;
for i = 1:5
    load(sprintf('sim_results/MO_R%d_E%d_t%d_w1part%d.mat',R,log10(n_sim),t_max,i));
    n_errors_tot = n_errors_tot+n_errors;
    n_sim_tot = n_sim_tot+sim_settings.n_sim;
end
FER = n_errors_tot./n_sim_tot;

% plot results
figure;
plot(1:sim_settings.t_max,FER,'*-');
xlabel('t enrollment observations');
ylabel('FER')
set(gca, 'YScale', 'log')
grid on
title('MO scheme')

% save as .txt
if storetable
    t_enroll = (1:t_max)';FER=FER(t_enroll)';
    T = table(t_enroll,FER); % let op het moeten columns zijn
    writetable(T,sprintf('sim_results/MOR%d.txt',R),'Delimiter','\t');
end


%% example plot results for sim_SMOhds
clear all
%close all        

% settings
n_sim = 10^6;
R = 6; % 1/R
t_min = 4; % t0 enrollments
t_max = 7; % t_max = t0+ max reconstructions
storetable = true; % store results in .txt file for plotting with pgfplots

n_errors_tot = 0; % total nr. errors
n_sim_tot = 0;
hist_n_fails_tot = 0; % hist. of simulated SRAMPUFs that failed 0, 1, 2, ..
% reconstructions
for i = 1:5
    load(sprintf('sim_results/SMO_R%d_E%d_t%d_t%d_w1part%d.mat',R,log10(n_sim),t_min,t_max,i));
    n_errors_tot = n_errors_tot+n_errors;
    hist_n_fails_tot = hist_n_fails_tot + hist_n_fails;
    n_sim_tot = n_sim_tot+sim_settings.n_sim;
end
FER = n_errors_tot./n_sim_tot;
% plot results
figure;
plot(1:sim_settings.t_max,FER,'*-');
xlabel('t0 (enrollments) + alpha (reconstructions)');
ylabel('FER')
set(gca, 'YScale', 'log')
grid on
title('SMO scheme')

%hist_n_fails_tot % verify whether there are SRAMPUFs that are 'bad', i.e.
% SRAMPUFs that fail reconstruction all the time and thus never update HD
    
if storetable
    t_enroll = (t_min:t_max)';FER=FER(t_enroll)';
    t_reconstruct = (1:length(t_enroll))';
    T = table(t_enroll,t_reconstruct,FER); % let op het moeten columns zijn
    writetable(T,sprintf('sim_results/SMOR%d.txt',R),'Delimiter','\t');
end




