function [tlist, Ilist, I_all_matrix, t_end, I_end] = SimulateStochasticBirthDeathExact_selection(epi_params, I_all_matrix, I_matrix_curr)

% initializing current population states:
t_curr = epi_params.timeStart; I_curr = epi_params.I_init;

tlist = []; Ilist = [];

ID = I_curr + 1; t_strobe = Inf;

while 1

    if t_strobe > epi_params.time_strobe
        tlist = [tlist; t_curr]; Ilist = [Ilist; I_curr]; 
        t_strobe = 0;
        [t_curr I_curr]
    end
    
    % specify rates:
    w(1) = epi_params.lambda*I_curr;         % birth
    w(2) = epi_params.delta*I_curr;         % death 
    
    if sum(w) == 0 % no more infecteds in population
        break;
    end
    
    deltaT = (1/sum(w))*log(1/rand);
    t_curr = t_curr + deltaT;
    t_strobe = t_strobe + deltaT;
    
    locs = find(cumsum(w)/sum(w) >= rand); event = locs(1);
    switch event
        case 1  % infection
            k_mutations = I_matrix_curr(:,5);
            fitness_I_curr = (1-epi_params.sigma).^k_mutations; % Haigh
            chosen_loc = randsample(length(fitness_I_curr'),1,true,fitness_I_curr');
            parent_ID = I_matrix_curr(chosen_loc,1);
            n_muts_at_transmission = poissrnd(epi_params.m);
            n_muts_cumulative = n_muts_at_transmission + I_matrix_curr(chosen_loc,5);
            insert = [ID t_curr parent_ID n_muts_at_transmission n_muts_cumulative];
            I_matrix_curr = [I_matrix_curr; insert];
            ID = ID + 1;
            I_curr = I_curr + 1;
        case 2   % death
            chosen_loc = randi(I_curr);         % location in I_matrix_curr which is the parent
            insert = [I_matrix_curr(chosen_loc,:) t_curr 0];
            I_all_matrix = [I_all_matrix; insert];
            I_matrix_curr(chosen_loc,:) = [];
            I_curr = I_curr - 1;
        otherwise
            error('not a valid event')
    end

    if t_curr > epi_params.timeEnd
        tlist = [tlist; t_curr]; Ilist = [Ilist; I_curr]; 
        break;
    end
end

t_end = t_curr;
I_end = I_curr; 

% shift all currently infected individuals into I_all_matrix, indicating their death time as INF (beyong the time of the simulation):
insert_matrix = [I_matrix_curr Inf*ones(size(I_matrix_curr,1),1) zeros(size(I_matrix_curr,1),1)];

I_all_matrix = [I_all_matrix; insert_matrix];
