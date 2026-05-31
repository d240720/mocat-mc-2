clc; clear;

addpath(genpath('/Users/dchen/Documents/GitHub/mocat-mc-2'))
addpath(genpath('/Users/dchen/Documents/GitHub/MOCAT-SSEM'));

ICfile = '2020.mat';
seed = 1;

%% Run single MC simulation
cfg = setup_MCconfig(seed, ICfile);
[nS, nD, nN, nB, mat_sats] = main_mc(cfg, seed);

%% Load output data
matfile = sprintf('TLEIC_year%i_rand%i.mat', cfg.time0.Year, seed);
data = load(matfile);

n_time = size(data.numObjects, 1);
t_years = (0:n_time-1)' * cfg.dt_days / 365.25;

S  = sum(data.S_MC, 2);   % active satellites
D  = sum(data.D_MC, 2);   % derelicts
N  = sum(data.N_MC, 2);   % debris
total = S + D + N;

%% Plot
figure('Position', [100 100 1000 700]);

subplot(2,2,1);
plot(t_years, total, 'k', 'LineWidth', 1.5);
xlabel('Time (years)'); ylabel('Count');
title('Total Objects'); grid on;

subplot(2,2,2);
plot(t_years, S, 'b', 'LineWidth', 1.5); hold on;
plot(t_years, D, 'r', 'LineWidth', 1.5);
xlabel('Time (years)'); ylabel('Count');
title('Active Satellites vs Derelicts');
legend('Active','Derelict'); grid on;

subplot(2,2,3);
plot(t_years, N, 'm', 'LineWidth', 1.5);
xlabel('Time (years)'); ylabel('Count');
title('Debris'); grid on;

subplot(2,2,4);
plot(t_years, double(data.count_coll), 'r', 'LineWidth', 1.5); hold on;
plot(t_years, double(data.count_expl), 'b', 'LineWidth', 1.5);
xlabel('Time (years)'); ylabel('Cumulative Count');
title('Collisions vs Explosions');
legend('Collisions','Explosions'); grid on;

sgtitle(sprintf('MOCAT-MC Single Run (Seed %d, IC: %s)', seed, ICfile));