function void = main_simulateBDTransmissionTree(void)

clear all; close all;

n_replicates = 1; %10;
n_samples = 200;

mlist = 0.3; %[0.3 3.0];
sigmalist = [0]; %[0 0.02 0.04 0.08 0.16 0.32 0.64];

t_sampleBegin = 0; %50;
t_sampleEnd = 55;

for m = mlist
    epi_params.m = m;                   % per genome, per transmission mutation rate, from Park, Martin, + Koelle (2023) Nature Communications

    for s = sigmalist
        epi_params.sigma = s;               % deleterious fitness effect
        
        epi_params.filename_pre = strcat('simData_tEnd55_m', int2str(epi_params.m*100), '_sel', int2str(epi_params.sigma*100), '_');

        for i = 1:n_replicates

            infile = strcat(epi_params.filename_pre, int2str(i), '_Imatrix');

            outfile_treeData_pre = strcat(epi_params.filename_pre, 'tStart', int2str(t_sampleBegin), '_', int2str(i), '_treeData_nSeqs');            
            outfile_tree_pre = strcat(epi_params.filename_pre, 'tStart', int2str(t_sampleBegin), '_', int2str(i), '_tree_nSeqs');

            BuildBDTransmissionTree(infile, outfile_treeData_pre, outfile_tree_pre, n_samples, t_sampleEnd, t_sampleBegin);
        end
    end

end
