# -*- coding: utf-8 -*-
"""
    Run linear birth-death (exp growth) sims with different mutational fitness effects, mutation rates and sampling intervals
    Note: This is the main script used to run all sims present in manuscript
"""
import subprocess
import tskit
import pyslim
import msprime
import numpy as np
import pandas as pd

def add_neutral_muts(ts,mu):
    
    """
        Add neutral mutations onto tree in TreeSequence ts 
        Mut rate mu is per site per unit time. 
        But need to divide by tic rate to adjust for time discretization step i.e. rate=mu/10
        
        Note: this was not used in paper
    """            
   
    print(f"The old tree sequence had {ts.num_mutations} mutations,\n"
          f"and mean pairwise nucleotide diversity is {ts.diversity():0.3e}.")
    
    next_id = pyslim.next_slim_mutation_id(ts)
    nts = msprime.sim_mutations(
        ts,
        rate=mu,
        model=msprime.SLiMMutationModel(type=0, next_id=next_id),
        keep=True,
    )

    print(f"The new tree sequence now has {nts.num_mutations} mutations,\n"
          f"and mean pairwise nucleotide diversity is {nts.diversity():0.3e}.")
    
def fitness_var(ts,time):
    
    """
        Compute fitness variance in TreeSequence ts at time
        Time t should be in backwards time from present (t=0)
        Note: only works if fitness effects are multiplicative         
    """
    
    tic_rate = 10 # integration dt time steps per unit time
    alive = pyslim.individuals_alive_at(ts, time*tic_rate) # in backwards time
    fit_alive = []
    for i in alive:
        u = ts.individual(i).nodes[0]
        muts = branch_muts[u]
        multi_fit = 1.0 # multiplicative fitness
        for m in muts:
            multi_fit *= 1 + ts.mutation(m).metadata['mutation_list'][0]['selection_coeff']
        fit_alive.append(multi_fit)
    fit_var = np.var(fit_alive)
    
    return fit_var

# Set up rng for seeding SLiM sims
rng = np.random.default_rng()

# Set num of sims to run
n_sims = 2

# Run SLiM simulation and load in the resulting trees
model_script = "./slim/linear_birth_death.slim"

# Init dataframe for sim/tree-level summary stats
stats_cols = ['sim_index',
                'mut_rate',
                'sel_coeff',
                'sampling_start',
                'num_samples',
                'num_muts',
                'total_branch_length',
                'sackin_index',
                'norm_sackin',
                'pi',
                'seg_sites',
                'tajimas_d',
                'num_singletons',
                'external_clock_rate',
                'internal_clock_rate',
                'external_muts',
                'internal_muts',
                'external_length',
                'internal_length',
                ]
stats_df = pd.DataFrame(columns=stats_cols,index=['sim_index'])

# Sim param values
mu_vals = [1.0e-04, 1.0e-05] # per site per generation mutation rates
s_vals = [0.0,-0.02,-0.04,-0.08,-0.16,-0.32,-0.64] # mutational fitness costs (s_d)
sampling_start_times = [0,35,50]

# Set min/max number of allowable samples per simulation
min_samples = 100
max_samples = 10000

end_time = 55 # time to end sims
tic_rate = 10 # integration steps per unit time

for mu in mu_vals:
    
    if mu == 1.0e-04:
        model_script = "./slim/linear_birth_death_high_mu.slim"
        
    if mu == 1.0e-05:
        model_script = "./slim/linear_birth_death_low_mu.slim"

    for s in s_vals:
    
        param_str = 's=' + str(s)
        
        for sampling_start in sampling_start_times:
    
            for sim in range(n_sims):
            
                print(f"Starting sim number {sim} with mu={mu} and s={s}\n")
            
                num_samples = 0
                while num_samples < min_samples or num_samples > max_samples:
                    
                    seed = rng.integers(1e6)
                    subprocess.check_output(["slim", "-s", str(seed), "-d", param_str, model_script])
                    ts = tskit.load("nonWF_haploid.trees")
        
                    print(f"The tree sequence has {ts.num_trees} trees\n"
                          f"on a genome of length {ts.sequence_length},\n"
                          f"{ts.num_individuals} individuals, {ts.num_samples} 'sample' genomes,\n"
                          f"and {ts.num_mutations} mutations.")
                
                    # Simplify ts to only include first "haploid" genome of each indiivudal
                    keep_nodes = []
                    node_times = []
                    for i in range(len(ts.tables.individuals)):
                        keep_nodes.append(ts.individual(i).nodes[0]) # only keep first node or genome
                        node_times.append(ts.nodes_time[ts.individual(i).nodes[0]])
                    keep_nodes = np.unique(keep_nodes)
                    #ts = ts.simplify(keep_nodes, keep_input_roots=True)
                    
                    # Convert node time units and to forward time
                    node_times = end_time - (np.array(node_times) / tic_rate)
                    
                    # Sample exactly 200 indvs after sampling start time
                    keep_nodes = keep_nodes[node_times >= sampling_start]
                    if len(keep_nodes) < 200:
                        print("Number of sampleable nodes too small: " , str(len(keep_nodes)))
                        continue
                    else:
                        keep_nodes = rng.choice(keep_nodes, 200, replace=False)
                    
                    # Simplify again!!!!
                    ts = ts.simplify(keep_nodes, keep_input_roots=True)
                    tree = ts.first()
                    
                    # Plot tree if desired
                    #svg = tree.draw_svg(path='test_nonWF.svg',size=(500, 500),time_scale="rank")
                    
                    num_samples = tree.num_samples() # Total num of samples in tree
                    print("Sim_size: " , str(num_samples))
                
                """
                    Compute pop gen and tree summary stats
                """
                
                num_muts = tree.num_mutations # Total num of mutations across tree
                
                total_branch_length = tree.total_branch_length / tic_rate # Sum of all branch lengths
                
                #colless = tree.colless_index() # does not work for trees with multifurcations
                
                sackin = tree.sackin_index() # Sackin's imbalance statisitc (unnormalized)
                exp_sackin = 2 *  num_samples * np.log(num_samples) # Expected Sackins's index under Yule model
                norm_sackin = (sackin - exp_sackin) / exp_sackin # Normalized Sackin's index
                
                pi = ts.diversity() # Avg pairwise diversity across the genome
                
                seg_sites = ts.segregating_sites() # Density of segregating sites
                
                tajimas_d = ts.Tajimas_D() # Tajima's D
                
                afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
                num_singletons = afs[1] # Num of singletons (muts on external branches)
                
                # Compute internal/external branch stats
                branch_lengths = [tree.branch_length(u) for u in tree.nodes()]
                branch_muts = {u:[] for u in tree.nodes()}
                for mt in tree.mutations():
                     branch_muts[mt.node].append(mt.id)
                external_rates = [] # external branch clock rate
                internal_rates = [] # internal branch clock rates
                external_muts = 0 # mutations on external branches
                internal_muts = 0 # mutations on internal branches
                external_length = 0 # sum of external branch lengths
                internal_length = 0 # sum of internal branch lengths
                for u in tree.nodes():
                    length = tree.branch_length(u) / tic_rate
                    count = len(branch_muts[u])
                    if length > 0: 
                        rate = count / length
                    else:
                        continue
                    if tree.is_leaf(u):
                        external_rates.append(rate)
                        external_muts += count
                        external_length += length
                    else:
                        internal_rates.append(rate)
                        internal_muts += count
                        internal_length += length
                external_clock_rate = np.mean(external_rates) # mean external branch clock rate
                internal_clock_rate = np.mean(internal_rates) # # mean internal branch clock rate
                

                stats_dict = {'sim_index': sim,
                              'mut_rate': mu,
                              'sel_coeff': s,
                              'sampling_start': sampling_start,
                              'num_samples': num_samples,
                              'num_muts': num_muts,
                              'total_branch_length':  total_branch_length,
                              'sackin_index': sackin,
                              'norm_sackin': norm_sackin,
                              'pi': pi,
                              'seg_sites': seg_sites,
                              'tajimas_d': tajimas_d,
                              'num_singletons': num_singletons,
                              'external_clock_rate':external_clock_rate,
                              'internal_clock_rate':internal_clock_rate,
                              'external_muts':external_muts,
                              'internal_muts':internal_muts,
                              'external_length':external_length,
                              'internal_length':internal_length,
                              }
                sim_df = pd.DataFrame(stats_dict, columns=stats_df.columns, index=['sim_index'])
                stats_df = pd.concat([stats_df, sim_df], ignore_index=True)

# Save results to csv
stats_df = stats_df.tail(-1)
summary_file = 'test_summary_stats_varying_mu_sd_ss_55days.csv'
stats_df.to_csv(summary_file,index=False)
