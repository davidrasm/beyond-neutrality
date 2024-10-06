function void = BuildBDTransmissionTree(infile, outfile_treeData_pre, outfile_tree_pre, n_samples, t_sampleEnd, t_sampleBegin)

load(infile);

I_matrix(1,3) = NaN;
cumulative_muts = I_matrix(:,6);

indiv_sampled_col_poss = intersect(find(I_matrix(:,2) > t_sampleBegin), find(I_matrix(:,2) < t_sampleEnd));
samp = randperm(length(indiv_sampled_col_poss),n_samples);
indiv_sampled_col = indiv_sampled_col_poss(samp);

seq_times_col = I_matrix(indiv_sampled_col,2);   % determine the times at which these individuals were sampled
indiv_sampled = indiv_sampled_col'; seq_times = seq_times_col';
for i = 1:n_samples
    names{i} = strcat('sample', int2str(i), '_', num2str(seq_times(i)));
    names_numerical(i) = NaN;
end
    
% get each sampled individual's genealogy of who-infected-whom: 
parentLineages = GetIndividualsLineages(n_samples, indiv_sampled, I_matrix(:,3));

n_seqs = n_samples;

% allocating memory to the structures that form the tree:
b = NaN*ones((n_seqs - 1),2);
d = NaN*ones(2*n_seqs - 1, 1);
loc_b = 1;

m = NaN*ones(2*n_seqs - 1, 1);  % number of mutations on branches

curr_indiv = indiv_sampled;

complete_indiv = indiv_sampled;
complete_seq_times = seq_times;

while 1
    
    num_lineages_remaining = n_seqs
    
    % finds which of the remaining lineages should coalesce next:
    [coal_daughters, coal_parent, timeOfCoalescence] = FindMostRecentCoalescence(curr_indiv, parentLineages, I_matrix(:,1), I_matrix(:,2));

    if timeOfCoalescence == -Inf % then no coalescences remain
        break;
    end
    
    indiv1_index = find(complete_indiv == coal_daughters(1));
    indiv2_index = find(complete_indiv == coal_daughters(2));

    d(indiv1_index(end), 1) = complete_seq_times(indiv1_index(end)) - timeOfCoalescence;
    d(indiv2_index(end), 1) = complete_seq_times(indiv2_index(end)) - timeOfCoalescence;

    m(indiv1_index(end), 1) = cumulative_muts(coal_daughters(1)) - cumulative_muts(coal_parent);
    m(indiv2_index(end), 1) = cumulative_muts(coal_daughters(2)) - cumulative_muts(coal_parent);

    b(loc_b, 1:2) = [indiv1_index(end) indiv2_index(end)];
    names{n_samples + loc_b} = strcat('node', int2str(loc_b), '_', num2str(timeOfCoalescence));
    names_numerical(n_samples + loc_b) = loc_b;
    
    loc_b = loc_b + 1;
    
    complete_indiv = [complete_indiv coal_parent];
    complete_seq_times = [complete_seq_times timeOfCoalescence];
    
    % erase coal_daughters from curr_indiv
    loc1_inCurr = find(curr_indiv == coal_daughters(1));
    curr_indiv(loc1_inCurr(1)) = [];
    loc2_inCurr = find(curr_indiv == coal_daughters(2));
    curr_indiv(loc2_inCurr(1)) = [];
    
    % add coal_parent to curr_indiv list
    curr_indiv = [curr_indiv coal_parent];
    
    % redo parentLineages
    n_seqs = length(curr_indiv);
    clear parentLineages
    parentLineages = GetIndividualsLineages(n_seqs, curr_indiv, I_matrix(:,3));
    
end
d(end) = 0;
m(end) = 0;

tree = phytree(b,d,names);

outfile_treeData = strcat(outfile_treeData_pre, int2str(n_samples));
outfile_tree = strcat(outfile_tree_pre, int2str(n_samples), '.tree');

save(outfile_treeData, 'b', 'd', 'm', 'names', 'names_numerical', 'epi_params');

phytreewrite(outfile_tree, tree, 'BranchNames', false)
