function void = main_birthDeathSimulation_selection(void)

clear all; close all; clc;

format long g

epi_params.n_replicates = 1; %10;           % number of replicate simulations under this parameterization
epi_params.timeStart = 0;               % start time is 0; end time depends on when the n^th individual is sampled
epi_params.time_strobe = 1;             % log variables every time_strobe units of time (in days)
epi_params.timeEnd = 55;                % start time is 0; end time depends on when the n^th individual is sampled

epi_params.lambda = 0.3;                % per capita birth rate (per day)
epi_params.delta = 0.1;                 % per capita death rate (i.e.. becoming uninfectious rate) (per day); delta = (psi + mu), where mu is the recovery through clearance and death and psi is the sampling rate
epi_params.I_init = 1;                  % initial infected population size

mlist = 0.3; %[0.3 3.0];
sigmalist = 0; %[0 0.02 0.04 0.08 0.16 0.32 0.64];

for m = mlist
    epi_params.m = m;   % per genome, per transmission mutation rate, from Park, Martin, + Koelle (2023) Nature Communications

    for s = sigmalist
        epi_params.sigma = s;       % deleterious fitness effect
        
        epi_params.filename_pre = strcat('simData_tEnd55_m', int2str(epi_params.m*100), '_sel', int2str(epi_params.sigma*100), '_');

        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % derived parameters:

        epi_params.R0 = epi_params.lambda/epi_params.delta;  % basic reproduction number (= 3; see Koelle et al. (2022) Science)
        epi_params.r = epi_params.lambda - epi_params.delta; % intrinsic growth rate (= 0.20 per day; see Koelle et al. (2022) Science)

        epi_params

        for i = 1:epi_params.n_replicates
            epi_params.filename = strcat(epi_params.filename_pre, int2str(i));
            SimulateBirthDeathExact_selection(epi_params);
        end
    end
end
