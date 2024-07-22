#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 13:38:17 2024

@author: david
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
#import itertools

df = pd.read_csv('../results/2024-05-13_summary_stats_varying_mu_sd_ss_55days.csv')
#df = pd.read_csv('linear_birth_test_summary_stats_gammaDFE.csv')
#df = pd.read_csv('../results/linear_birth_test_summary_stats_high_mu.csv')
#df = pd.read_csv('linear_birth_test_summary_stats_high_mu.csv')

sns.set_theme(style="darkgrid")
sns.set_palette('Set2',desat=0.5) #, n_colors=None, desat=None

#order = [0.02,0.04,0.08,0.16,0.32,0.64]
df['sel_coeff'] = df['sel_coeff'].abs() # put all s_d on a positve scale
order = [0.0,0.02,0.04,0.08,0.16,0.32,0.64]

text_size = 14
#plt.rcParams.update({'font.size': text_size})

# For adding std dev values
# sd_vals = [0.02]*100
# sd_vals.extend([0.04]*100)
# sd_vals.extend([0.08]*100)
# sd_vals.extend([0.16]*100)
# sd_vals.extend([0.32]*100)
# sd_vals.extend([0.64]*100)
# df['Std dev'] = sd_vals
#df.drop(df.loc[df['sel_coeff']==-0.64].index, inplace=True)

df = df[df['mut_rate']==0.00001]

#sampling_label = 'Sampling interval'
#df = df.rename(columns={'sampling_start': 'Sampling interval'})
#df['Sampling interval'] = df['Sampling interval'].map({0: 'Full', 35: 'Last 20 days', 50:'Last 5 days'})
#df = df[df['Sampling interval']=='Full']
sampling_label = None

"""
    Plot standard pop gen / tree summary stats
"""
x_var = 'sel_coeff'

# Number samples
# fig, ax = plt.subplots(figsize=(8,5))
# sns.violinplot(data=df, x=x_var, y="num_samples", density_norm='width', order=order, hue=sampling_label)
# png_file = 'lbd_varying_s_sampling_window_num-samples.png' # + f'{param:.2f}' +'.png'
# fig.tight_layout()
# fig.savefig(png_file, dpi=200)

# Number mutations
#fig, ax = plt.subplots(figsize=(8,5))
#sns.violinplot(data=df, x=x_var, y="num_muts", density_norm='width', order=order)

# Total branch length
#fig, ax = plt.subplots(figsize=(8,5))
#sns.violinplot(data=df, x=x_var, y="total_branch_length", density_norm='width', order=order)

# Sackin index
# fig, ax = plt.subplots(figsize=(5,3.5))
# sns.violinplot(data=df, x=x_var, y="sackin_index", density_norm='width', order=order, hue=sampling_label)
# ax.set_xlabel('Fitness cost $s_d$',fontsize=text_size)
# ax.set_ylabel('Sackin Index',fontsize=text_size)
# ax.set_ylim([300, 2300])
# png_file = 'lbd_varying_s_sampling_window_sackin.png' # + f'{param:.2f}' +'.png'
# fig.tight_layout()
# fig.savefig(png_file, dpi=200)

# Avg. pair. div. (pi)
# fig, ax = plt.subplots(figsize=(8,5))
# sns.violinplot(data=df, x=x_var, y="pi", density_norm='width', order=order, hue=sampling_label)
# png_file = 'lbd_varying_s_sampling_window_pi.png' # + f'{param:.2f}' +'.png'
# fig.tight_layout()
# fig.savefig(png_file, dpi=200)

# Density of segregating sites
#fig, ax = plt.subplots(figsize=(8,5))
#sns.violinplot(data=df, x=x_var, y="seg_sites", density_norm='width', order=order)
#png_file = 'lbd_varying_s_density-seg-sites.png' # + f'{param:.2f}' +'.png'
#fig.tight_layout()
#fig.savefig(png_file, dpi=200)

# Tajima's D
# fig, ax = plt.subplots(figsize=(8,5))
# sns.violinplot(data=df, x=x_var, y="tajimas_d", density_norm='width', order=order, hue=sampling_label)
# ax.set_xlabel('Fitness effect $s_d$',fontsize=text_size)
# ax.set_ylabel('Tajimas D',fontsize=text_size)
# png_file = 'lbd_varying_s_sampling_window_tajimas.png' # + f'{param:.2f}' +'.png'
# fig.tight_layout()
# fig.savefig(png_file, dpi=200)

# Number singletons (muts on external branches)
# fig, ax = plt.subplots(figsize=(8,5))
# sns.violinplot(data=df, x="sel_coeff", y="num_singletons", density_norm='width', order=order)
# ax.set_xlabel('Fitness effect $s_d$')
# ax.set_ylabel('Number of singletons')

# Propotion singletons (muts on external branches)
# df["prop_singletons"] = df["num_singletons"] / df["num_muts"]
# fig, ax = plt.subplots(figsize=(8,5))
# sns.violinplot(data=df, x=x_var, y="prop_singletons", density_norm='width', order=order, hue=sampling_label)
# ax.set_xlabel('Fitness effect $s_d$',fontsize=text_size)
# ax.set_ylabel('Propotion singletons',fontsize=text_size)
# png_file = 'lbd_varying_s_sampling_window_prop-singletons.png' # + f'{param:.2f}' +'.png'
# fig.tight_layout()
# fig.savefig(png_file, dpi=200)

# External clock rate (muts on external branches)
df["external_internal_clock_ratio"] = np.log10(df["external_clock_rate"] / df["internal_clock_rate"])
df.drop(df.loc[df['external_internal_clock_ratio']==np.Inf].index, inplace=True)
fig, ax = plt.subplots(figsize=(5,3.5))
sns.violinplot(data=df, x=x_var, y="external_internal_clock_ratio", density_norm='width', order=order, hue=sampling_label) 
ax.set_xlabel('Fitness cost $s_d$',fontsize=text_size)
ax.set_ylabel('Log10 External / Internal Clock Ratio',fontsize=text_size)
png_file = 'lbd_varying_s_external-vs-internal-clock.png' # + f'{param:.2f}' +'.png'
fig.tight_layout()
fig.savefig(png_file, dpi=200)

# Proportion of tree length external
df["external_length_prop"] = df["external_length"] / (df["internal_length"] +  df["external_length"])
# fig, ax = plt.subplots(figsize=(5,3.5))
# sns.violinplot(data=df, x=x_var, y="external_length_prop", density_norm='width', order=order, hue=sampling_label)
# ax.set_xlabel('Fitness cost $s_d$',fontsize=text_size)
# ax.set_ylabel('Proportion Tree Length External',fontsize=text_size)
# ax.set_ylim([.6, 1.])
# png_file = 'lbd_varying_s_sampling_window_external_length_proportion.png' # + f'{param:.2f}' +'.png'
# fig.tight_layout()
# fig.savefig(png_file, dpi=200)

# Proportion of mutations external
df["external_muts_prop"] = df["external_muts"] / (df["internal_muts"] +  df["external_muts"])
fig, ax = plt.subplots(figsize=(5,3.5))
sns.violinplot(data=df, x=x_var, y="external_muts_prop", density_norm='width', order=order, hue=sampling_label)
ax.set_xlabel('Fitness cost $s_d$',fontsize=text_size)
ax.set_ylabel('Proportion Mutations External',fontsize=text_size)
png_file = 'lbd_varying_s_external_muts_proportion.png' # + f'{param:.2f}' +'.png'
fig.tight_layout()
fig.savefig(png_file, dpi=200)

# Compare the average proportion of mutations that are external (eta_e / S) to the average proportion of tree length composed of external branches (Ln/Jn).
# Ratio of proportion of external mutations to proportion of external tree lenth
df["external_prop_ratio"] = df["external_muts_prop"] / df["external_length_prop"]
fig, ax = plt.subplots(figsize=(5,3.5))
sns.violinplot(data=df, x=x_var, y="external_prop_ratio", density_norm='width', order=order, hue=sampling_label)
ax.set_xlabel('Fitness cost $s_d$',fontsize=text_size)
ax.set_ylabel('Ratio Mutations / Tree Length External',fontsize=text_size)
png_file = 'lbd_varying_s_ratio_external_legnth_to_muts.png' # + f'{param:.2f}' +'.png'
fig.tight_layout()
fig.savefig(png_file, dpi=200)

# Observed fit variation
# fig, ax = plt.subplots(figsize=(5,3.5))
# sns.violinplot(data=df, x=x_var, y="obsv_fit_var", density_norm='width', order=order, hue=sampling_label)
# png_file = 'lbd_varying_s_sampling_window_obsv_fit_var.png' # + f'{param:.2f}' +'.png'
# ax.set_xlabel('Fitness cost $s_d$',fontsize=text_size)
# ax.set_ylabel('Fitness variation',fontsize=text_size)
# fig.tight_layout()
# fig.savefig(png_file, dpi=200)



