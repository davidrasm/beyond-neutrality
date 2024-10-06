function data = GetSequencesDownstream_poisson(data, node_num, init_seq, b, d, names, names_numerical, epi_params)

branchorleaf1 = b(node_num, 1);
branchorleaf2 = b(node_num, 2);
branchLength1 = d(branchorleaf1);
m1 = poissrnd(branchLength1*epi_params.subsPerGenomePerDay);
seq1 = MutateSequence(init_seq, m1);
branchLength2 = d(branchorleaf2);
m2 = poissrnd(branchLength2*epi_params.subsPerGenomePerDay);
seq2 = MutateSequence(init_seq, m2);
    
data(branchorleaf1).Sequence = seq1;
data(branchorleaf1).Header = names{branchorleaf1};
node_num1 = names_numerical(branchorleaf1);

data(branchorleaf2).Sequence = seq2;
data(branchorleaf2).Header = names{branchorleaf2};
node_num2 = names_numerical(branchorleaf2);

if branchorleaf1 > epi_params.n 
    data = GetSequencesDownstream_poisson(data, node_num1, data(branchorleaf1).Sequence, b, d, names, names_numerical, epi_params);
end
if branchorleaf2 > epi_params.n 
    data = GetSequencesDownstream_poisson(data, node_num2, data(branchorleaf2).Sequence, b, d, names, names_numerical, epi_params);
end
