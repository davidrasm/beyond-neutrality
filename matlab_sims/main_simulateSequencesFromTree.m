function void = main_simulateSequencesFromTree(void)

clear all; close all; clc;

n_replicates = 1; %10;
n_samples = 200;

init_seq_struct = getgenbank('NC_045512');
init_seq = init_seq_struct.Sequence;
seq_length = length(init_seq)

mlist = 0.3; %[0.3 3.0];
sigmalist = 0; %[0 0.02 0.04 0.08 0.16 0.32 0.64];

t_sampleBegin = 0;
t_sampleEnd = 55;

for m = mlist
    epi_params.m = m;                    % per genome, per transmission mutation rate, from Park, Martin, + Koelle (2023) Nature Communications

    for s = sigmalist
        epi_params.sigma = s;                % deleterious fitness effect
        
        epi_params.filename_pre = strcat('simData_tEnd55_m', int2str(epi_params.m*100), '_sel', int2str(epi_params.sigma*100), '_');

        for i = 1:n_replicates

            infile = strcat(epi_params.filename_pre, 'tStart', int2str(t_sampleBegin), '_', int2str(i), '_treeData_nSeqs', int2str(n_samples));

            outfile_fasta = strcat(infile, '_sequences.fasta');
            outfile_fasta_poisson = strcat(infile, '_sequences_poisson.fasta');
    
            SimulateSequencesFromTree(infile, init_seq, outfile_fasta, n_samples);

            SimulateSequencesFromTree_Poisson(infile, init_seq, outfile_fasta_poisson, n_samples);
        end
    end
end


