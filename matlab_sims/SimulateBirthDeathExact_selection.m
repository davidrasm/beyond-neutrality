function void = SimulateBirthDeathExact_selection(epi_params)

% I_all_matrix: matrix of all individuals who have been infected over the simulation
% ** ONLY CONTAINS INDIVIDUALS WHO HAVE ALREADY RECOVERED FROM INFECTION
% ** WILL BE USED IN CONSTRUCTING THE TRANSMISSION TREE
% column 1: infection number ("ID")
% column 2: start time of infection ("birth time")
% column 3: individual who transmitted the disease to this individual ("parent's ID")
% column 4: number of mutations that occurred during transmission from parent
% column 5: cumulative number of mutations that occurred since index case
% column 6: end time of infection ("death time")
% NOT NEEDED (MAKE 0 FOR ALL): column 7: why infections ended: 1 = sampling (rate psi; with probability prob_s); 2 = recovery (rate mu)

% I_matrix_curr: matrix of I individuals who are currently infected
% column 1: infection number ("ID")
% column 2: start time of infection ("birth time")
% column 3: individual who transmitted the disease to this individual ("parent's ID")
% column 4: number of mutations that occurred during transmission from parent
% column 5: cumulative number of mutations that occurred since the root

I_matrix_curr(1:epi_params.I_init,1) = (1:epi_params.I_init)';  % ID
I_matrix_curr(1:epi_params.I_init,2) = NaN;  % birth time of index case is N/A, infected at start of simulation
I_matrix_curr(1:epi_params.I_init,3) = 0;    % parent is defined as 0
I_matrix_curr(1:epi_params.I_init,4) = 0;    % no mutations (irrelevant)
I_matrix_curr(1:epi_params.I_init,5) = 0;    % no cumulative mutations (irrelevant)

% simulate until you get an effective invasion
while 1
    [tlist, Ilist, I_all_matrix, t_end, I_end] = SimulateStochasticBirthDeathExact_selection(epi_params, [], I_matrix_curr);
    if I_end > 0  % did not go stochastically extinct
        break;
    end
end

% sort I_all_matrix so more amenable to reconstructing the transmission tree later:
[sorted_vals, sorted_indices] = sort(I_all_matrix(:,1));
n_entries = length(sorted_vals); I_ordered_matrix(1:n_entries, 1:7) = NaN;
for i = 1:n_entries
    I_ordered_matrix(i,:) = I_all_matrix(sorted_indices(i),:);
end
% create I_matrix that has columns: birth time, death time, parent ID, sampling status, mutations from parent, cumulative mutations from index case. Row number is implicitly the parent ID number.
I_matrix = [I_ordered_matrix(:,2) I_ordered_matrix(:,6) I_ordered_matrix(:,3) I_ordered_matrix(:,7) I_ordered_matrix(:,4) I_ordered_matrix(:,5)];

epi_params.timeEnd = t_end;

Imatrix_file = strcat(epi_params.filename, '_Imatrix.mat');
timeSeries_file = strcat(epi_params.filename, '_timeSeries.mat');

save(Imatrix_file, 'epi_params', 'I_matrix');
save(timeSeries_file, 'epi_params', 'tlist', 'Ilist')
