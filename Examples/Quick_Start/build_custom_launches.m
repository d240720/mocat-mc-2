function repeatLaunches = build_custom_launches(cfg, varargin)
% BUILD_CUSTOM_LAUNCHES  Construct a synthetic repeatLaunches mat_sats matrix
%
% Usage:
%   repeatLaunches = build_custom_launches(cfg)
%   repeatLaunches = build_custom_launches(cfg, 'launches_per_year', 200, ...)
%
% Required input:
%   cfg   - MCconfig struct (needs cfg.missionlifetime, cfg.YEAR2DAY, cfg.time0)
%
% Optional name-value pairs:
%   'pool_size'           Total objects in repeat pool    [default: 750]
%                         Effective rate = pool_size / (launchRepeatYrs span in years)
%                         e.g. pool_size=1000 with [2018,2022] => 200 launches/yr
%   'alt_bins'            Altitude bin centers [km]       [default: [400,550,700,1000,1200]]
%   'alt_weights'         Fraction of launches per bin    [default: equal]
%   'mass_mean'           Mean satellite mass [kg]        [default: 150]
%   'mass_std'            Std dev of mass [kg]            [default: 50]
%   'radius_mean'         Mean satellite radius [m]       [default: 0.5]
%   'pmd_compliance'      Fraction of sats that are controlled [default: 1.0]
%   'missionlife'         Mission lifetime [years]        [default: cfg.missionlifetime]
%   'seed'                RNG seed for reproducibility    [default: 42]
%
% Output:
%   repeatLaunches  - mat_sats matrix [N x 24] representing one year of launches
%                     to be repeated cyclically by main_mc matsat mode
%
% Column layout (standard MOCAT mat_sats):
%   1:  idx_a           semi-major axis [Earth radii]
%   2:  idx_ecco        eccentricity
%   3:  idx_inclo       inclination [rad]
%   4:  idx_nodeo       RAAN [rad]
%   5:  idx_argpo       argument of perigee [rad]
%   6:  idx_mo          mean anomaly [rad]
%   7:  idx_bstar       bstar drag term
%   8:  idx_mass        mass [kg]
%   9:  idx_radius      radius [m]
%   10: idx_error       error flag
%   11: idx_controlled  controlled flag [0/1]
%   12: idx_a_desired   desired semi-major axis [Earth radii]
%   13: idx_missionlife mission lifetime [years]
%   14: idx_constel     constellation flag [0/1]
%   15: idx_date_created date created [Julian date]
%   16: idx_launch_date  launch date [Julian date]
%   17-19: idx_r        position vector [km]
%   20-22: idx_v        velocity vector [km/s]
%   23: idx_objectclass object class [1=payload]
%   24: idx_ID          object ID
%
% Example — 200 launches/yr, mostly at 550km and 1200km:
%   rl = build_custom_launches(cfg, ...
%       'launches_per_year', 200, ...
%       'alt_bins',    [400, 550, 700, 1200], ...
%       'alt_weights', [0.1, 0.4, 0.2, 0.3]);
%
% Daniel Chen / research use — MOCAT-MC Quick_Start compatible
% -------------------------------------------------------------

    %% Parse inputs
    p = inputParser;
    addRequired(p, 'cfg');
    addParameter(p, 'pool_size', 750);
    addParameter(p, 'alt_bins',    [400, 550, 700, 1000, 1200]);
    addParameter(p, 'alt_weights', []);   % empty = uniform
    addParameter(p, 'mass_mean',   150);
    addParameter(p, 'mass_std',    50);
    addParameter(p, 'radius_mean', 0.5);
    addParameter(p, 'pmd_compliance', 1.0);
    addParameter(p, 'missionlife', cfg.missionlifetime);
    addParameter(p, 'seed', 42);
    parse(p, cfg, varargin{:});
    opt = p.Results;

    rng(opt.seed);

    %% Derived quantities
    Re      = 6378.135;          % Earth radius [km]
    mu      = 3.986004418e5;     % GM [km^3/s^2]
    N       = opt.pool_size;
    n_bins  = numel(opt.alt_bins);
    n_repeat_years = diff(cfg.launchRepeatYrs) + 1;
    effective_rate = N / n_repeat_years;

    % Normalise altitude weights
    if isempty(opt.alt_weights)
        weights = ones(1, n_bins) / n_bins;
    else
        assert(numel(opt.alt_weights) == n_bins, ...
            'alt_weights must have same length as alt_bins');
        weights = opt.alt_weights / sum(opt.alt_weights);
    end

    % How many sats per altitude bin
    counts = round(weights * N);
    counts(end) = N - sum(counts(1:end-1));  % fix rounding residual

    %% Base Julian date — spread launches over full launchRepeatYrs window
    jd0 = juliandate(datetime(cfg.launchRepeatYrs(1), 1, 1));
    jd_end = juliandate(datetime(cfg.launchRepeatYrs(2), 12, 31));
    span_days = jd_end - jd0;
    YEAR2DAY = cfg.YEAR2DAY;
    fprintf('launch date window: %i-Jan-%i to %i-Dec-%i (%.0f days)\n', ...
        1, cfg.launchRepeatYrs(1), 31, cfg.launchRepeatYrs(2), span_days);

    %% Build mat_sats rows
    repeatLaunches = zeros(N, 24);
    row = 1;

    for b = 1:n_bins
        nb   = counts(b);
        if nb <= 0; continue; end
        alt  = opt.alt_bins(b);           % km
        a_km = Re + alt;                  % semi-major axis [km]
        a_er = a_km / Re;                 % semi-major axis [Earth radii]

        % Orbital velocity magnitude [km/s]
        v_circ = sqrt(mu / a_km);

        % Random orbital elements
        inclo  = deg2rad(unifrnd(0,   90,  nb, 1));   % inclination 0-90 deg
        nodeo  = deg2rad(unifrnd(0,   360, nb, 1));   % RAAN
        argpo  = deg2rad(unifrnd(0,   360, nb, 1));   % arg of perigee
        mo     = deg2rad(unifrnd(0,   360, nb, 1));   % mean anomaly
        ecco   = zeros(nb, 1);                         % circular orbit

        % Physical properties
        mass   = max(1, opt.mass_mean + opt.mass_std * randn(nb,1));
        radius = max(0.01, opt.radius_mean * ones(nb,1));

        % bstar: approximate from mass and radius
        % bstar = 0.5 * Cd * A / m * rho_ref (MOCAT convention)
        bstar  = 0.5 * 2.2 * radius.^2 ./ mass * 0.157;

        % Controlled flag — Bernoulli draw
        controlled = double(rand(nb,1) < opt.pmd_compliance);

        % Spread launch dates uniformly over full launchRepeatYrs window
        launch_jd = jd0 + rand(nb,1) * span_days;

        % Position/velocity: approximate circular orbit in ECI
        % (main_mc recomputes these via SGP4/propagator anyway)
        r_mag = a_km;
        r_vec = [r_mag * ones(nb,1), zeros(nb,1), zeros(nb,1)];
        v_vec = [zeros(nb,1), v_circ * ones(nb,1), zeros(nb,1)];

        % Assemble rows
        for i = 1:nb
            repeatLaunches(row, :) = [ ...
                a_er,           ...  % 1  idx_a
                ecco(i),        ...  % 2  idx_ecco
                inclo(i),       ...  % 3  idx_inclo
                nodeo(i),       ...  % 4  idx_nodeo
                argpo(i),       ...  % 5  idx_argpo
                mo(i),          ...  % 6  idx_mo
                bstar(i),       ...  % 7  idx_bstar
                mass(i),        ...  % 8  idx_mass
                radius(i),      ...  % 9  idx_radius
                0,              ...  % 10 idx_error
                controlled(i),  ...  % 11 idx_controlled
                a_er,           ...  % 12 idx_a_desired (same as a for circular)
                opt.missionlife,...  % 13 idx_missionlife
                0,              ...  % 14 idx_constel
                jd0,            ...  % 15 idx_date_created
                launch_jd(i),   ...  % 16 idx_launch_date
                r_vec(i,1), r_vec(i,2), r_vec(i,3), ...  % 17-19 idx_r
                v_vec(i,1), v_vec(i,2), v_vec(i,3), ...  % 20-22 idx_v
                1,              ...  % 23 idx_objectclass (1 = payload)
                row + 1000000   ...  % 24 idx_ID (offset to avoid collisions with IC IDs)
            ];
            row = row + 1;
        end
    end

    fprintf('build_custom_launches: %i objects (%.0f/yr) across %i altitude bins\n', N, effective_rate, n_bins);
    for b = 1:n_bins
        fprintf('  alt=%4ikm  n=%i  (%.0f%%)\n', ...
            opt.alt_bins(b), counts(b), weights(b)*100);
    end
end