function void = SimulateSequencesFromTree(infile, init_seq, outfile_seqs, n_samples)

load(infile)

epi_params.n = n_samples;

seq_length = length(init_seq);

init_node = length(d);

data(init_node).Sequence = init_seq;
data(init_node).Header = names{init_node};

curr_node = init_node;
node_num = names_numerical(curr_node);

data = GetSequencesDownstream(data, node_num, init_seq, b, d, m, names, names_numerical, epi_params);

for i = 1:epi_params.n
    tip_data(i).Sequence = data(i).Sequence;
    tip_data(i).Header = data(i).Header;
end

fastawrite(outfile_seqs, tip_data)
