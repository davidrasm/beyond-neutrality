function void = SimulateSequencesFromTree_Poisson(infile, init_seq, outfile_fasta_poisson, n_samples)

load(infile)

epi_params.n = n_samples;

seq_length = length(init_seq);

epi_params.subsPerGenomePerDay = epi_params.m*(1/(1/epi_params.lambda));
plot(d, m, 'k.') % number of mutations on branches of x-axis length
m_poisson = poissrnd(d*epi_params.subsPerGenomePerDay)
hold on; plot(d, m_poisson, 'ro')

init_node = length(d);

data(init_node).Sequence = init_seq;
data(init_node).Header = names{init_node};

curr_node = init_node;
node_num = names_numerical(curr_node);

data = GetSequencesDownstream_poisson(data, node_num, init_seq, b, d, names, names_numerical, epi_params);

for i = 1:epi_params.n
    tip_data(i).Sequence = data(i).Sequence;
    tip_data(i).Header = data(i).Header;
end

fastawrite(outfile_fasta_poisson, tip_data)
