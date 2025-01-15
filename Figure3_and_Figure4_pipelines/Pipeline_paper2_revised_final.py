########################################################################################################################
# SCRIPT CLUSTER, CORRELATION, AND LDA ANALYSIS OF CIRCADIAN RHYTMICITIES IN TRIPLE NEGATIVE BREAST CANCER (TNBC) ######
# Collaboration with The Granada Lab from Charite, Berlin. ###### Code written by Jeff DIDIER, Sebastien DE LANDTSHEER #
########################################################################################################################

# SysBio Group, DLSM, FSTC, University of Luxembourg
# Jeff Didier
# Dr Sebastien De Landtsheer
# Prof Thomas Sauter

# Revision 30.09.2024: new sensitivity data by collaborators (<ctrl-f>, search for '## REVISED')

# Revision 27.11.2024: response to reviewer's major comments (<ctrl-f>, search for '## FINAL REVISION')

# Revision 13.01.2025: data extraction of specific figures for reproduction (Fig 3 and 4)

#########################
# ## LIBRARIES IMPORTS ##

import os
import random
import matplotlib

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from kneed import KneeLocator
from adjustText import adjust_text
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch, Polygon
from sklearn.linear_model import LogisticRegression, Ridge, Lasso
from sklearn.decomposition import PCA
from functions.data_functions import make_mut_table, remove_columns_with_zeros, get_LDA_metrics, \
    LDA_loocv, LogReg_loocv, regressor_loocv, LogReg_on_LDA_loocv, LogRegLDA_on_PCA_loocv
from functions.plot_functions import plot_UMAP, plot_LDA, plot_linear_regressor, plot_regressor_loocv
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

##########################################################################
# ## CREATE FIX VARIABLES, RESULTS FOLDER, AND COLLECT META INFORMATION ##
plt.switch_backend('agg')
formats = ['svg', 'png']
dpi = 300
seed = 42

# Import data files
data_bmal1 = pd.read_excel("../data/circadian_parameters_by_replicate_Bmal1_final_v1.xlsx", sheet_name=None)
data_per2 = pd.read_excel("../data/circadian_parameters_by_replicate_Per2_final_v1.xlsx", sheet_name=None)
data_bmal1_per2 = pd.read_excel("../data/circadian_parameters_by_replicate_Bmal1_Per2_final_v1.xlsx", sheet_name=None)
data_growth = pd.read_excel("../data/growth_and_drug_sensitivity_values_final_v3.xlsx", sheet_name='growthrate')
data_sensitivty = pd.read_excel("../data/growth_and_drug_sensitivity_values_final_v3.xlsx", index_col=0,
                                sheet_name=['GEC50', 'GR50', 'GRAOC', 'GRINF', 'Hill'])  # latest sensitivity data
ccle_drug_data = pd.read_csv('../data/CCLE_NP24.2009_Drug_data_2015.02.24.csv')
data_exp = pd.read_csv("../data/CCLE_expression_from_CCLE_2022_11_22.csv")
data_info = pd.read_csv("../data/sample_info.csv").set_index('CCLE_Name')
short_ccg_list = pd.read_excel("../data/circadian_clock_gene_lists_v2.xlsx", sheet_name="genelist_short",
                               skiprows=1, header=None, names=['Gene', 'AlternativeName'])  # short gene list (16 genes)

# Collect cell line, all subtypes, and unique subtypes
celllines = data_bmal1_per2['classification'].T[0].values[1:]
subtypes = data_bmal1_per2['classification'].columns[1:]
unique_subtypes = pd.unique([x.split('.')[0] for x in subtypes])

# For our convenience, a dict to easily map cell line to its subtype
dict_to_map_line_and_subtype = dict({k: v for (k, v) in zip(data_bmal1_per2['classification'].T[0].values[1:],
                                                            [x.split('.')[0] for x in subtypes])})
sample_tag = []
for k, v in dict_to_map_line_and_subtype.items():
    sample_tag.append(v + '_' + k)

# Get all existing parameters, remove 'classification', and add 'growthrate'
all_circ_params = list(data_bmal1_per2.keys())
all_circ_params.remove('classification')
circ_growth_params = all_circ_params + ['growthrate']

# Create the matrices for each channel that we feed with the mean of the replicates (Skipna=True)
df_matrix_bmal1 = pd.DataFrame(index=sample_tag, columns=circ_growth_params)
df_matrix_per2 = pd.DataFrame(index=sample_tag, columns=circ_growth_params)
df_matrix_bmal1_per2 = pd.DataFrame(index=sample_tag, columns=circ_growth_params)

# do the same for median
medians_bmal1 = pd.DataFrame(index=sample_tag, columns=circ_growth_params)
medians_per2 = pd.DataFrame(index=sample_tag, columns=circ_growth_params)
medians_bmal1_per2 = pd.DataFrame(index=sample_tag, columns=circ_growth_params)

# For loop to feed the matrices with the correct data
for n, channel in enumerate([data_bmal1, data_per2, data_bmal1_per2]):
    for param, df in channel.items():
        for cell in celllines:
            if param in circ_growth_params and param != 'classification' and cell in df.columns:
                this_line = [x for x in sample_tag if cell in x and dict_to_map_line_and_subtype[cell] in x]
                if n == 0:  # BMAL1
                    df_matrix_bmal1.loc[this_line, param] = df[cell].mean()
                    medians_bmal1.loc[this_line, param] = df[cell].median()
                    if 'U2OS' not in cell:   # dont have growthrate for the U2OS cell lines
                        df_matrix_bmal1.loc[this_line, 'growthrate'] = data_growth[cell].values
                        medians_bmal1.loc[this_line, 'growthrate'] = data_growth[cell].values
                if n == 1:  # PER2
                    df_matrix_per2.loc[this_line, param] = df[cell].mean()
                    medians_per2.loc[this_line, param] = df[cell].median()
                    if 'U2OS' not in cell:   # dont have growthrate for the U2OS cell lines
                        df_matrix_per2.loc[this_line, 'growthrate'] = data_growth[cell].values
                        medians_per2.loc[this_line, 'growthrate'] = data_growth[cell].values
                if n == 2:  # BMAL1_PER2
                    df_matrix_bmal1_per2.loc[this_line, param] = df[cell].mean()
                    medians_bmal1_per2.loc[this_line, param] = df[cell].median()
                    if 'U2OS' not in cell:   # dont have growthrate for the U2OS cell lines
                        df_matrix_bmal1_per2.loc[this_line, 'growthrate'] = data_growth[cell].values
                        medians_bmal1_per2.loc[this_line, 'growthrate'] = data_growth[cell].values

# list of features to be removed for reduced data frame
feature_to_remove = \
    ['ridgelength', 'period_coeffvar', 'expdecay_envelope', 'mra_noise', 'mra_ultradian', 'phdiff_median']

# cleaning CCLE expression ##
# filter data: only cell lines we need
DepMapId = [data_info[data_info['stripped_cell_line_name'] == x.split('_')[0]].loc[
            :, 'DepMap_ID'].to_numpy()[0] for x in celllines[:15]]  # excluding osteosarcoma, [:16] to include
data_exp_filtered = data_exp.loc[data_exp['Unnamed: 0'].isin(DepMapId)].set_index('Unnamed: 0')
# remove genes for which there is more  than one "0.00"
data_exp_filtered = remove_columns_with_zeros(data_exp_filtered, 2)
# remove the number from the gene name
colnames = data_exp_filtered.columns
colnames_correct = [x.split()[0] for x in colnames]
data_exp_filtered.columns = colnames_correct

# expression data - normalizing the expression (min-max) [ xscaled = (x - min(x)) / (max(x) - min(x)) ]
data_exp_norm = (data_exp_filtered-data_exp_filtered.min())/(data_exp_filtered.max()-data_exp_filtered.min())
data_exp_norm.index = \
    [data_info['stripped_cell_line_name'][data_info['DepMap_ID'] == x][0] for x in data_exp_norm.index]

# short and ccg dataframes (all genes from short list are available in CCLE)
short_ccg_norm = data_exp_norm[short_ccg_list['Gene']]
X_short_unscaled = data_exp_filtered[short_ccg_list['Gene']]
X_short = short_ccg_norm  # final short core clock genes expression (16 genes)

# ## formatting the CCLE drug data set
ccle_drug_data['Primary Cell Line Name'] = ccle_drug_data['Primary Cell Line Name'].replace('-', '', regex=True)
ccle_drug_data['Compound'] = ccle_drug_data['Compound'].replace('-', '', regex=True)

# find the available cell lines in CCLE drug data
remaining_ccle_drugs = pd.DataFrame(columns=ccle_drug_data.columns)
for line in celllines:
    if line in set(ccle_drug_data['Primary Cell Line Name']):
        remaining_ccle_drugs = \
            pd.concat([remaining_ccle_drugs, ccle_drug_data.loc[ccle_drug_data['Primary Cell Line Name'] == line]])

remaining_drugs = remaining_ccle_drugs['Compound'].unique()  # no Cisplatin remaining...

#########################################################
# ## RESULTS FOLDER CREATIONS FOR THE VARIOUS VERSIONS ##

# general results folder
if os.path.isdir('./results') is False:
    os.mkdir('./results')

# general results folder revision 30.09.2024
if os.path.isdir('./results_revised') is False:
    os.mkdir('./results_revised')

# general results folder revision final 27.11.2024
if os.path.isdir('./results_revised_review') is False:
    os.mkdir('./results_revised_review')

# since the revision requests some more plots, we will organize them by comment
for sub_fold in ['comment_2', 'comment_3', 'comment_5', 'comment_5_alternatives']:
    if os.path.isdir(f'./results_revised_review/{sub_fold}') is False:
        os.mkdir(f'./results_revised_review/{sub_fold}')

# comment 5 leads to a lot of figures, and will be subdivided into the input data and output sensitivity sources
for comm_5_sub_fold in ['CCLE_in_Cisplatin_out', 'Circadian_in_Cisplatin_out',
                        'CCLE_in_CCLE_out', 'Circadian_in_CCLE_out']:
    if os.path.isdir(f'./results_revised_review/comment_5/{comm_5_sub_fold}') is False:
        os.mkdir(f'./results_revised_review/comment_5/{comm_5_sub_fold}')
    # sub folders for the binarized versions
    for comm_5_algo_path in ['lda', 'logreg', 'logreg_on_lda']:
        if os.path.isdir(f'./results_revised_review/comment_5/{comm_5_sub_fold}/{comm_5_algo_path}') is False:
            os.mkdir(f'./results_revised_review/comment_5/{comm_5_sub_fold}/{comm_5_algo_path}')
        # sub folder of lasso and ridge for cisplatin data only
        if 'Cisplatin' in comm_5_sub_fold:
            for comm_5_algo_cont in ['ridge', 'lasso']:
                if os.path.isdir(f'./results_revised_review/comment_5/{comm_5_sub_fold}/{comm_5_algo_cont}') is False:
                    os.mkdir(f'./results_revised_review/comment_5/{comm_5_sub_fold}/{comm_5_algo_cont}')

#####################################
# ## Figure 3 a), and Figure S4 a) ##

# ## REVISED

data_here = df_matrix_bmal1_per2
figure_panels = ['3_a', 's4_a']

for elem in figure_panels:
    df = data_here.copy()
    if elem == '3_a':  # reduced feature set
        df.drop(columns=[f for f in feature_to_remove if f in df.columns], inplace=True)
    df = df.dropna(axis=1, how='all')
    df = df.dropna(axis=0, how='all')
    # save corr and pvals as excel
    corr_matrix = df.corr(method='pearson', numeric_only=False)
    pvals_matrix = df.corr(method=lambda x, y: pearsonr(x, y)[1], numeric_only=False) - np.eye(len(df.columns))
    with pd.ExcelWriter(f'./results_revised/{elem}.xlsx', mode='w') as writer:
        corr_matrix.to_excel(writer, sheet_name='correlations')
    with pd.ExcelWriter(f'./results_revised/{elem}.xlsx', mode='a') as writer:
        pvals_matrix.to_excel(writer, sheet_name='pvalues')
    # plotting
    fig = plt.figure(figsize=(12, 12))
    ax = plt.gca()
    im = ax.matshow(corr_matrix, cmap=plt.get_cmap('coolwarm'), interpolation='none', vmin=-1, vmax=1)
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(range(len(corr_matrix.columns)), corr_matrix.columns, fontsize=10, rotation=45, ha='right')
    ax.set_yticks(range(len(corr_matrix.columns)), corr_matrix.columns, fontsize=10)
    ax.set_title(f"Circadian parameter relationships")
    ax.tick_params(axis=u'both', which=u'both', length=0)
    for (i, j), z in np.ndenumerate(np.array(pvals_matrix)):
        if (i != j or z != 0) and z < 0.05:
            star_multiplier = int(('%.2e' % z)[-1]) - 1
            if star_multiplier > 4 or ('%.2e' % z)[-2:] == 10:
                star_multiplier = 4  # block if 4 or above
            z = '*' * star_multiplier
            ax.text(j, i, "%s\n%.1f" % (z, corr_matrix.iloc[j, i]), ha="center", va="center",
                    fontsize=18 if elem == '3_a' else 16,
                    color='w' if np.abs(float('%.1f' % corr_matrix.iloc[j, i])) >= 0.7 else 'k')
    cbar = fig.colorbar(im, shrink=.6)
    cbar.ax.set_ylabel(r"Pearson $r$", rotation=90)
    cbar.ax.get_yaxis().set_ticks([-1, 0, 1])
    cbar.ax.yaxis.set_label_position('left')
    plt.tight_layout()
    # save figs
    for form in formats:
        plt.savefig(f'./results_revised/{elem}.{form}', bbox_inches='tight', dpi=dpi)
    plt.close()

#####################################
# ## Figure 3 b), and Figure S4 b) ##

# ## REVISED

width = 0.3  # for the barplots
best_components = 5  # enough to describe 95% of the variance

random.seed(seed)
np.random.seed(seed)

data_here = df_matrix_bmal1_per2
matrix = data_here.copy()
df = matrix.dropna(axis=1, how='all')
df = df.dropna(axis=0)  # need to have all possible nans removed
df.drop(columns=[f for f in feature_to_remove if f in df.columns], inplace=True)
# normalize by min-max
df_norm = (df - df.min()) / (df.max() - df.min())
# Get pca of the data
pca = PCA(n_components=best_components, random_state=seed)
pca.fit(df_norm)
PCs = pca.fit_transform(df_norm)
PCdf = pd.DataFrame(data=PCs, columns=["PC" + str(i) for i in range(1, PCs.shape[1] + 1)])
PCdf.index = df_norm.index
# Match PC names to loadings
pc_loadings = dict(zip(PCdf.columns, pca.components_))
# Matrix of corr coefficients between pcs and features
loadings_df = pd.DataFrame.from_dict(pc_loadings)
loadings_df['feature_names'] = df_norm.columns
loadings_df = loadings_df.set_index('feature_names')

# plot the loading plot with the scatterplot
target_markers = {'TNBC-BL1': 'o',
                  'TNBC-BL2': '<',
                  'TNBC-M': 'd',
                  'BC-LumA': '^',
                  'Epithelial': 's'}
target_colors = {'TNBC-BL1': 'brown',
                 'TNBC-BL2': 'orange',
                 'TNBC-M': 'green',
                 'BC-LumA': 'blue',
                 'Epithelial': 'black'}
xs = pca.components_[0]
ys = pca.components_[1]
targets = pd.DataFrame(index=PCdf.index, columns=['subtype'], data=[x.split('_')[0] for x in df_norm.index])
g = sns.lmplot(x='PC1', y='PC2', data=pd.concat([PCdf, targets], axis=1), fit_reg=False, hue='subtype',
               palette=target_colors, markers=list(target_markers.values()),
               scatter_kws={'linewidths': 1, 'edgecolor': 'k', "s": 75})
sns.move_legend(g, "upper right", bbox_to_anchor=(1.02, 0.7), frameon=True)
text_points = []
text_arrows = []
# label points
for ind in PCdf.index:
    text_points.append(plt.annotate(ind.split('_')[1], (PCdf.loc[ind]['PC1'], PCdf.loc[ind]['PC2']), fontsize=9,
                       ha='center', va='center', weight='bold'))
# label arrows
for i, varnames in enumerate(df_norm.columns):
    plt.arrow(0, 0, xs[i], ys[i], color='k', head_width=0.03, alpha=0.3, linestyle='-', linewidth=2)
    text_arrows.append(plt.text(xs[i], ys[i], varnames, fontsize=7))
plt.xlabel(f'PC1 ({"{:.1f}".format(round(pca.explained_variance_ratio_[0] * 100, 2))}%)')
plt.ylabel(f'PC2 ({"{:.1f}".format(round(pca.explained_variance_ratio_[1] * 100, 2))}%)')
plt.title(f'Circadian cell model mapping')
plt.tight_layout()
adjust_text(text_arrows)
adjust_text(text_points, force_text=(1.1, 2), expand=(1.1, 1.9))
for form in formats:
    plt.savefig(f'./results_revised/3_b.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

# PCA bars, here seed reset is not required as we continue to use transformed data from the biplot
original_hatch_width = 1.0
matplotlib.rcParams['hatch.linewidth'] = 8.0
fig, ax = plt.subplots(figsize=(8, 6))
loadings_df = loadings_df.reindex(loadings_df['PC1'].abs().sort_values(ascending=False).index)  # sort by absolute PC1
rect1 = ax.bar(np.arange(len(loadings_df.index)) - width/2, loadings_df["PC1"].abs(), width=width, label="PC1",
               color='steelblue', edgecolor="k", hatch=['\\' if i < 0 else None for i in loadings_df["PC1"]])
rect2 = ax.bar(np.arange(len(loadings_df.index)) + width/2,
               loadings_df["PC2"].abs(), width=width, label="PC2", edgecolor="k",
               color='lemonchiffon', hatch=['\\' if i < 0 else None for i in loadings_df["PC2"]])
ax.set_xticks(np.arange(len(loadings_df.index)))
ax.set_xticklabels(loadings_df.index)
ax.tick_params(axis="x", labelsize=10, labelrotation=45)
lgd = ax.legend(labels=['PC1', 'PC2'], fontsize=12, ncol=3)
handles, labs = lgd.axes.get_legend_handles_labels()
handles.append(Patch(facecolor='lightgrey', edgecolor='k', hatch='\\'))
labs.append('Negative')
lgd._legend_box = None
lgd._init_legend_box(handles, labs)
lgd._set_loc(lgd._loc)
lgd.set_title(lgd.get_title().get_text())
for num, ha in enumerate(lgd.legendHandles):
    if num < len(lgd.legendHandles) - 1:
        ha.set_hatch(None)
plt.xlabel('')
plt.ylabel('Absolute PC loading')
plt.title('Ranking of PC loadings')
plt.tight_layout()
for form in formats:
    plt.savefig(f'./results_revised/s4_b.{form}', bbox_inches='tight', dpi=300)
plt.close()
matplotlib.rcParams['hatch.linewidth'] = original_hatch_width

# write the loadings of PC1, PC2, PC3, PC4, PC5 in separate worksheets
with pd.ExcelWriter(f'./results_revised/s4_b.xlsx') as writer:
    for pc in ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']:
        loadings_df[pc][loadings_df[pc].abs().sort_values(ascending=False).index].to_excel(writer, sheet_name=pc)

# ## FINAL REVISION (Comment 3)
# - show box(en) plot distributions of original data, min-max data, standard-scaled data across cell lines

data_here = df_matrix_bmal1_per2
matrix = data_here.copy()
df = matrix.dropna(axis=1, how='all')
df = df.dropna(axis=0)  # need to have all possible nans removed
df.drop(columns=[f for f in feature_to_remove if f in df.columns], inplace=True)
# unscaled
df_unscaled = df
# normalize by standard scaling
df_standard = (df - df.mean()) / df.std()
# normalize by min-max
df_minmax = (df - df.min()) / (df.max() - df.min())
# boxen plot of features unscaled vs minmax vs standard
df_unscaled['scaling'] = 'raw data'
df_standard['scaling'] = 'standard scaling'
df_minmax['scaling'] = 'min-max scaling'
# combine all
combined_df = pd.concat([df_unscaled, df_standard, df_minmax])
# melt for boxenplot
melted_df = pd.melt(combined_df.reset_index(), id_vars=['index', 'scaling'], var_name='feature', value_name='value')
hue_order = melted_df['feature'].unique()
palette = sns.color_palette("deep", len(hue_order))
hue_colors = dict(zip(hue_order, palette))

# plot
fig, axes = plt.subplots(figsize=(14, 6), nrows=1, ncols=3)
for ax, scale in zip(axes, ['raw data', 'standard scaling', 'min-max scaling']):
    ax.yaxis.grid(True)
    data_here = melted_df[melted_df['scaling'] == scale]
    # Box plot to get the error bars, without fliers
    g = sns.boxplot(ax=ax, data=data_here, x='feature', y='value', hue='feature', saturation=1, linewidth=.8,
                    linecolor='dimgrey', legend=True if scale == 'raw data' else False, zorder=2, showfliers=False,
                    palette=palette)
    for hue_value, color in hue_colors.items():
        feature_here = data_here[data_here['feature'] == hue_value]
        # Boxen plot to have the fliers, avoid overlapping of out-of-error-bar flier, no legend needed here
        sns.boxenplot(ax=ax, data=feature_here, x='feature', y='value', hue='feature', saturation=1, linewidth=.8,
                      ec='dimgrey', legend=False, zorder=2, flier_kws=dict(facecolor=color, edgecolor='dimgrey'),
                      facecolor=color)
        if hue_value == 'amplitude_median' and scale == 'raw data':
            outliers = feature_here[feature_here['value'] >= 100]
            if len(outliers) >= 1:
                for out in outliers.iloc:
                    ax.annotate(out['index'].split('_')[1], (out['feature'], out['value']))
    if scale == 'raw data':
        g.legend_.set_title('$\\bf{Circadian\\ features}$')
        plt.setp(g.get_legend().get_texts(), fontsize='8')  # for legend text
        plt.setp(g.get_legend().get_title(), fontsize='8')  # for legend title
    ax.set_title(f'{scale}', fontsize=14)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.set_xticks([])
    ax.set_yticks(ticks=ax.get_yticks(), labels=ax.get_yticklabels(), fontsize=12)
    if scale != 'raw data':
        ax.set_ylim((melted_df[melted_df['scaling'] == scale].min().value - 0.1,
                     melted_df[melted_df['scaling'] == scale].max().value + 0.1))
fig.suptitle(f'Raw, min-max, and standard scaled circadian features (n={len(df_unscaled)})', fontsize=16)
fig.tight_layout()
for form in formats:
    plt.savefig(f'./results_revised_review/comment_3/feature_scaling.{form}', bbox_inches='tight', dpi=300)
plt.close()

# --------------------
# Alternative, showing all 9 features on their individual scale
fig = plt.figure(figsize=(16, 10))
second_row_axes = []
gs0 = fig.add_gridspec(2, 18)
second_row_axes.append(fig.add_subplot(gs0[1, :9]))
second_row_axes.append(fig.add_subplot(gs0[1, 9:]))
for scale in ['raw data', 'standard scaling', 'min-max scaling']:
    data_here = melted_df[melted_df['scaling'] == scale]
    if scale == 'raw data':
        for i, feat in enumerate(data_here['feature'].unique()):  # first row, all individual raw features
            ax = fig.add_subplot(gs0[0, i * 2:(i + 1) * 2])
            ax.yaxis.grid(True)
            # Box plot to get the error bars, without fliers
            sns.boxplot(ax=ax, data=data_here[data_here['feature'] == feat], x='feature', y='value',
                        color=hue_colors[feat], saturation=1, linewidth=.8, linecolor='dimgrey', legend=False, zorder=2,
                        showfliers=False)
            # Boxen plot to have the fliers, avoid overlapping of out-of-error-bar flier, no legend needed here
            sns.boxenplot(ax=ax, data=data_here[data_here['feature'] == feat], x='feature', y='value',
                          color=hue_colors[feat], saturation=1, linewidth=.8, ec='dimgrey', legend=False, zorder=2,
                          flier_kws=dict(facecolor=hue_colors[feat], edgecolor='dimgrey'))
            ax.set_xlabel(None)
            ax.set_ylabel(None)
            ax.set_xticks([])
            ax.set_yticks(ticks=ax.get_yticks(), labels=ax.get_yticklabels(), fontsize=12)
            if feat == 'amplitude_median':
                feature_here = data_here[data_here['feature'] == feat]
                outliers = feature_here[feature_here['value'] >= 100]
                if len(outliers) >= 1:
                    for out in outliers.iloc:
                        ax.annotate(out['index'].split('_')[1], (out['feature'], out['value']), fontsize=10)
        fig.suptitle(f'{scale}', fontsize=14, y=0.92)
    if scale != 'raw data':
        # second row, standard vs min-max scaled
        ax = second_row_axes[0] if scale == 'standard scaling' else second_row_axes[1]
        ax.yaxis.grid(True)
        # Box plot to get the error bars, without fliers
        g = sns.boxplot(ax=ax, data=data_here, x='feature', y='value', hue='feature', saturation=1, linewidth=.8,
                        linecolor='dimgrey', legend=True if scale == 'standard scaling' else False, zorder=2,
                        showfliers=False, palette=palette)
        for hue_value, color in hue_colors.items():
            feature_here = data_here[data_here['feature'] == hue_value]
            # Boxen plot to have the fliers, avoid overlapping of out-of-error-bar flier, no legend needed here
            sns.boxenplot(ax=ax, data=feature_here, x='feature', y='value', hue='feature', saturation=1,
                          linewidth=.8, ec='dimgrey', legend=False, zorder=2, facecolor=color,
                          flier_kws=dict(facecolor=color, edgecolor='dimgrey'))
        ax.set_title(f'{scale}', fontsize=14)
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticks([])
        ax.set_yticks(ticks=ax.get_yticks(), labels=ax.get_yticklabels(), fontsize=12)
        if scale == 'standard scaling':
            g.legend_.set_title('$\\bf{Circadian\\ features}$')
            plt.setp(g.get_legend().get_texts(), fontsize='6')  # for legend text
            plt.setp(g.get_legend().get_title(), fontsize='6')  # for legend title
            g.legend(loc='lower right', ncol=3)
fig.suptitle(f'Raw, min-max, and standard scaled circadian features (n={len(df_unscaled)})', fontsize=16)
fig.tight_layout()
for form in formats:
    plt.savefig(f'./results_revised_review/comment_3/feature_scaling_individual.{form}', bbox_inches='tight', dpi=300)
plt.close()
# --------------------

#  - rerun PCA on standard scaled data (careful with arrow scaling) (+ stats if possible, if so also for the min-max)

width = 0.3  # for the barplots
best_components = 5  # enough to describe 95% of the variance

random.seed(seed)
np.random.seed(seed)

data_here = df_matrix_bmal1_per2
matrix = data_here.copy()
df = matrix.dropna(axis=1, how='all')
df = df.dropna(axis=0)  # need to have all possible nans removed
df.drop(columns=[f for f in feature_to_remove if f in df.columns], inplace=True)
# normalize by standard scaling
df_norm = (df - df.mean()) / df.std()
# Get pca of the data
pca = PCA(n_components=best_components, random_state=seed)
pca.fit(df_norm)
PCs = pca.fit_transform(df_norm)
PCdf = pd.DataFrame(data=PCs, columns=["PC" + str(i) for i in range(1, PCs.shape[1] + 1)])
PCdf.index = df_norm.index
# Match PC names to loadings
pc_loadings = dict(zip(PCdf.columns, pca.components_))
# Matrix of corr coefficients between pcs and features
loadings_df = pd.DataFrame.from_dict(pc_loadings)
loadings_df['feature_names'] = df_norm.columns
loadings_df = loadings_df.set_index('feature_names')

# plot the loading plot with the scatterplot
target_markers = {'TNBC-BL1': 'o',
                  'TNBC-BL2': '<',
                  'TNBC-M': 'd',
                  'BC-LumA': '^',
                  'Epithelial': 's'}
target_colors = {'TNBC-BL1': 'brown',
                 'TNBC-BL2': 'orange',
                 'TNBC-M': 'green',
                 'BC-LumA': 'blue',
                 'Epithelial': 'black'}
xs = pca.components_[0] * PCdf['PC1'].max()  # scale loadings by max PC1 value
ys = pca.components_[1] * PCdf['PC2'].max()  # scale loadings by max PC2 value
targets = pd.DataFrame(index=PCdf.index, columns=['subtype'], data=[x.split('_')[0] for x in df_norm.index])

g = sns.lmplot(x='PC1', y='PC2', data=pd.concat([PCdf, targets], axis=1), fit_reg=False, hue='subtype',
               palette=target_colors, markers=list(target_markers.values()),
               scatter_kws={'linewidths': 1, 'edgecolor': 'k', "s": 75})
sns.move_legend(g, "upper right", bbox_to_anchor=(1.02, 0.98), frameon=True)
text_points = []
text_arrows = []
# label points
for ind in PCdf.index:
    text_points.append(plt.annotate(ind.split('_')[1], (PCdf.loc[ind]['PC1'], PCdf.loc[ind]['PC2']), fontsize=8,
                       ha='left', va='center', weight='bold'))
# label arrows
for i, varnames in enumerate(df_norm.columns):
    plt.arrow(0, 0, xs[i], ys[i], color='k', head_width=0.08, alpha=0.3, linestyle='-', linewidth=2)
    text_arrows.append(plt.text(xs[i], ys[i], varnames, fontsize=8))
plt.xlabel(f'PC1 ({"{:.1f}".format(round(pca.explained_variance_ratio_[0] * 100, 2))}%)')
plt.ylabel(f'PC2 ({"{:.1f}".format(round(pca.explained_variance_ratio_[1] * 100, 2))}%)')
plt.title(f'Circadian cell model mapping\n(component-wise scaling of loadings)')
plt.tight_layout()
adjust_text(text_arrows + text_points)
for form in formats:
    plt.savefig(f'./results_revised_review/comment_3/pca_standardscaled.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

# PCA bars, here seed reset is not required as we continue to use transformed data from the biplot
original_hatch_width = 1.0
matplotlib.rcParams['hatch.linewidth'] = 8.0
fig, ax = plt.subplots(figsize=(8, 6))
loadings_df = loadings_df.reindex(loadings_df['PC1'].abs().sort_values(ascending=False).index)  # sort by absolute PC1
rect1 = ax.bar(np.arange(len(loadings_df.index)) - width/2, loadings_df["PC1"].abs(), width=width, label="PC1",
               color='steelblue', edgecolor="k", hatch=['\\' if i < 0 else None for i in loadings_df["PC1"]])
rect2 = ax.bar(np.arange(len(loadings_df.index)) + width/2,
               loadings_df["PC2"].abs(), width=width, label="PC2", edgecolor="k",
               color='lemonchiffon', hatch=['\\' if i < 0 else None for i in loadings_df["PC2"]])
ax.set_xticks(np.arange(len(loadings_df.index)))
ax.set_xticklabels(loadings_df.index)
ax.tick_params(axis="x", labelsize=10, labelrotation=45)
lgd = ax.legend(labels=['PC1', 'PC2'], fontsize=12, ncol=3)
handles, labs = lgd.axes.get_legend_handles_labels()
handles.append(Patch(facecolor='lightgrey', edgecolor='k', hatch='\\'))
labs.append('Negative')
lgd._legend_box = None
lgd._init_legend_box(handles, labs)
lgd._set_loc(lgd._loc)
lgd.set_title(lgd.get_title().get_text())
for num, ha in enumerate(lgd.legendHandles):
    if num < len(lgd.legendHandles) - 1:
        ha.set_hatch(None)
plt.xlabel('')
plt.ylabel('Absolute PC loading')
plt.title('Ranking of PC loadings')
plt.tight_layout()
for form in formats:
    plt.savefig(f'./results_revised_review/comment_3/pca_standardscaled_bars.{form}', bbox_inches='tight', dpi=300)
plt.close()
matplotlib.rcParams['hatch.linewidth'] = original_hatch_width

# write the loadings of PC1, PC2, PC3, PC4, PC5 in separate worksheets
with pd.ExcelWriter(f'./results_revised_review/comment_3/pca_standardscaled.xlsx') as writer:
    for pc in ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']:
        loadings_df[pc][loadings_df[pc].abs().sort_values(ascending=False).index].to_excel(writer, sheet_name=pc)

###################
# ## Figure 3 c) ##

# ## REVISED

data_here = df_matrix_bmal1_per2

matrix = data_here.copy()
df = matrix.dropna(axis=1, how='all')
df = df.dropna(axis=0)  # to keep it consistent with exercise 2
df_norm = (df - df.min()) / (df.max() - df.min())

# obvious candidates in BMAL1_PER2
obvious_candidates_comb = ['period_median', 'phdiff_coeffvar', 'mra_circadian']

# Export scaled obvious candidates as xlsx
df_norm[obvious_candidates_comb].to_excel('./results_revised_review/MSB-2024-12581_SourceData_Fig3C.xlsx',
                                          sheet_name='sheet1', index=True)

f = plt.figure(figsize=(8, 8))
ax = f.add_subplot(1, 1, 1, projection='3d')
# plotting normalized would visually be the same as the raw data, so let's stay with the normalized
xx, yy, zz = df_norm[obvious_candidates_comb[0]], df_norm[obvious_candidates_comb[1]], \
             df_norm[obvious_candidates_comb[2]]
shapes = []
colors = []
for num, ind in enumerate(df_norm.index):
    shapes.append('^' if 'BC-LumA' in ind else 's' if 'Epithelial' in ind else
                  'o' if 'TNBC-BL1' in ind else '<' if 'TNBC-BL2' in ind else 'd')
    colors.append('lightseagreen' if ind.split('_')[-1] in ('MCF7', 'MCF10A', 'MDAMB468',
                                                            'T47D', 'HCC1937', 'HCC1143') else
                  'purple' if ind.split('_')[-1] in ('HCC1806', 'MDAMB231', 'CAL51') else
                  'violet' if ind.split('_')[-1] in ('HCC1428', 'MDAMB436', 'BT549', 'HCC70') else
                  'yellow')
    # texts
    txt = ind.split('_')[1]
    ax.text3D(xx[num], yy[num] * 1.05, zz[num] * 1.075 if num // 2 == 0 else zz[num] - (0.075 * zz[num]),
              txt, fontsize=8, fontweight='bold', ha='center', va='bottom', color='k', zorder=20)
    ax.plot([xx[num], xx[num]], [yy[num], yy[num] * 1.05],
            [zz[num], zz[num] * 1.075] if num // 2 == 0 else [zz[num], zz[num] - (0.075 * zz[num])],
            color='grey', linestyle='-', alpha=.75, lw=1, zorder=2)
for num, _ in enumerate(xx.values):
    ax.scatter(xx[num], yy[num], zz[num], s=100, color=colors[num], marker=shapes[num], edgecolor='k', zorder=0)
# lines
z2 = np.ones(shape=xx.shape) * min(df_norm[df_norm.columns[2]])
for i, j, k, h in zip(xx, yy, zz, z2):
    ax.plot([i, i], [j, j], [k, h], color='grey', linestyle='dashed', alpha=.75, lw=1)
ax.set_title('Defining circadian-based phenotypes', fontsize=14)
ax.set_xlabel('Period scale')
ax.set_ylabel('CV phase diff.\nscale')
ax.set_zlabel('Circadian band\nscale')
ax.set_xlim(0, 1)
ax.set_xticks(ticks=[0, 0.25, 0.5, 0.75, 1], labels=[0, 0.25, 0.5, 0.75, 1])
ax.set_ylim(0, 1)
ax.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1], labels=[0, 0.25, 0.5, 0.75, 1])
ax.set_zlim(0, 1)
ax.set_zticks(ticks=[0, 0.25, 0.5, 0.75, 1], labels=[0, 0.25, 0.5, 0.75, 1])
ax.set_box_aspect((1, 1, 1))
# doing the legend
legend_elements = []
for c, typ in zip(np.unique(colors), ['Functional', 'Unstable', 'Weak', 'Dysfunctional']):
    legend_elements.append(Patch(color=c, label=typ))
ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(-0.11, 1))

for form in formats:
    plt.savefig(f'./results_revised/3_c.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

# ## FINAL REVISION (Comment 2)

# - UMAP on obvious features with n_neighbors 2, 3, 4

seed = 42

data_here = df_matrix_bmal1_per2
matrix = data_here.copy()
df = matrix.dropna(axis=1, how='all')
df = df.dropna(axis=0)  # to keep it consistent with exercise 2
df_norm = (df - df.min()) / (df.max() - df.min())

# obvious candidates in BMAL1_PER2
obvious_candidates_comb = ['period_median', 'phdiff_coeffvar', 'mra_circadian']

targets = pd.DataFrame([x.split('_')[0] for x in df.index], index=df.index, columns=['subtype'])

# plot UMAP (neighbors=2, 3, 4)
for n in [2, 3, 4]:
    random.seed(seed)
    np.random.seed(seed)
    plot_UMAP(df_norm[obvious_candidates_comb], pd.concat([df_norm, targets], axis=1), 'subtype', n_neighbors=n,
              seed=seed, title=f'Circadian clock feature mapping\n(UMAP n_neighbors={n})', custom_color_marker=True,
              force_text=(1.2, 1.2))
    for form in formats:
        plt.savefig(f'./results_revised_review/comment_2/UMAP_on_obvious_candidates_n_{n}.{form}',
                    bbox_inches='tight', dpi=dpi)
    plt.close()

# - UMAP on PC1&2 (&3) with n_neighbors 2, 3, 4

seed = 42

random.seed(seed)
np.random.seed(seed)

data_here = df_matrix_bmal1_per2
matrix = data_here.copy()
df = matrix.dropna(axis=1, how='all')
df = df.dropna(axis=0)  # to keep it consistent with exercise 2
df_norm = (df - df.min()) / (df.max() - df.min())

pca = PCA(n_components=best_components, random_state=seed)
pca.fit(df_norm)
PCs = pca.fit_transform(df_norm)
PCdf = pd.DataFrame(data=PCs, columns=["PC" + str(i) for i in range(1, PCs.shape[1] + 1)])
PCdf.index = df_norm.index

targets = pd.DataFrame([x.split('_')[0] for x in df.index], index=df.index, columns=['subtype'])

# plot UMAP (neighbors=2, 3, 4)
for number_pcs in [2, 3]:
    for n in [2, 3, 4]:
        random.seed(seed)
        np.random.seed(seed)
        plot_UMAP(PCdf.iloc[:, :number_pcs], pd.concat([PCdf.iloc[:, :number_pcs], targets], axis=1), 'subtype',
                  n_neighbors=n, seed=seed, title=f'Circadian clock feature mapping\n'
                                                  f'(UMAP n_neighbors={n} on PC 1-{number_pcs})',
                  custom_color_marker=True, force_text=(1.2, 1.2))
        for form in formats:
            plt.savefig(f'./results_revised_review/comment_2/UMAP_on_PCA1-{number_pcs}_n_{n}.{form}',
                        bbox_inches='tight', dpi=dpi)
        plt.close()

# - k-means clustering, silhouette and WCSS on the obvious features with k 2, 3, 4, 5, 6

seed = 42

n_clusters = [2, 3, 4, 5, 6, 7]

data_here = df_matrix_bmal1_per2
matrix = data_here.copy()
df = matrix.dropna(axis=1, how='all')
df = df.dropna(axis=0)  # to keep it consistent with exercise 2
df_norm = (df - df.min()) / (df.max() - df.min())
df_obv = df_norm[obvious_candidates_comb]

scores = pd.DataFrame(columns=n_clusters, index=['silhouette_avg', 'wcss'])
for clusters in n_clusters:
    random.seed(seed)
    np.random.seed(seed)
    kmeans = KMeans(n_clusters=clusters, init='k-means++', random_state=seed)
    kmeans.fit(df_obv)
    clustered = kmeans.fit_predict(df_obv)
    # save silhouette and wcss scores
    scores.loc['wcss', clusters] = kmeans.inertia_  # wcss is inertia in kmeans
    scores.loc['silhouette_avg', clusters] = silhouette_score(df_obv, clustered)
    # plot k means clusters
    f = plt.figure(figsize=(8, 8))
    ax = f.add_subplot(1, 1, 1, projection='3d')
    # plotting normalized would visually be the same as the raw data, so let's stay with the normalized
    xx, yy, zz = \
        df_norm[obvious_candidates_comb[0]], df_norm[obvious_candidates_comb[1]], df_norm[obvious_candidates_comb[2]]
    shapes = []
    colors = sns.color_palette("deep", len(np.unique(clustered)))
    colors_this_cluster = [colors[n] for n in clustered]
    for num, ind in enumerate(df_norm.index):
        shapes.append('^' if 'BC-LumA' in ind else 's' if 'Epithelial' in ind else 'o' if 'TNBC-BL1' in ind else
                      '<' if 'TNBC-BL2' in ind else 'd')
        # texts
        txt = ind.split('_')[1]
        ax.text3D(xx[num], yy[num] * 1.05, zz[num] * 1.075 if num // 2 == 0 else zz[num] - (0.075 * zz[num]),
                  txt, fontsize=8, fontweight='bold', ha='center', va='bottom', color='k', zorder=20)
        ax.plot([xx[num], xx[num]], [yy[num], yy[num] * 1.05],
                [zz[num], zz[num] * 1.075] if num // 2 == 0 else [zz[num], zz[num] - (0.075 * zz[num])],
                color='grey', linestyle='-', alpha=.75, lw=1, zorder=2)
    for num, _ in enumerate(xx.values):
        ax.scatter(xx[num], yy[num], zz[num], s=100, color=colors_this_cluster[num], marker=shapes[num],
                   edgecolor='k', zorder=0)
    # lines
    z2 = np.ones(shape=xx.shape) * min(df_norm[df_norm.columns[2]])
    for i, j, k, h in zip(xx, yy, zz, z2):
        ax.plot([i, i], [j, j], [k, h], color='grey', linestyle='dashed', alpha=.75, lw=1)
    ax.set_title(f'Defining circadian-based phenotypes, clusters={clusters}\n'
                 f'(Silh.: {"%.2f" % silhouette_score(df_obv, clustered)}, WCSS: {"%.2f" % kmeans.inertia_})',
                 fontsize=14)
    ax.set_xlabel('Period scale')
    ax.set_ylabel('CV phase diff.\nscale')
    ax.set_zlabel('Circadian band\nscale')
    ax.set_xlim(0, 1)
    ax.set_xticks(ticks=[0, 0.25, 0.5, 0.75, 1], labels=[0, 0.25, 0.5, 0.75, 1])
    ax.set_ylim(0, 1)
    ax.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1], labels=[0, 0.25, 0.5, 0.75, 1])
    ax.set_zlim(0, 1)
    ax.set_zticks(ticks=[0, 0.25, 0.5, 0.75, 1], labels=[0, 0.25, 0.5, 0.75, 1])
    ax.set_box_aspect((1, 1, 1))
    # doing the legend
    legend_elements = []
    for c, typ in zip(colors, [f'cluster {x + 1}' for x in np.unique(clustered)]):
        legend_elements.append(Patch(color=c, label=typ))
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(-0.11, 1))
    for form in formats:
        plt.savefig(f'./results_revised_review/comment_2/kmeans_clustering_n_{clusters}.{form}',
                    bbox_inches='tight', dpi=dpi)
    plt.close()

# find elbow point
knee_locator = KneeLocator(n_clusters, scores.loc['wcss'], curve="convex", direction="decreasing", S=1)
elbow_point = knee_locator.knee
if elbow_point is not None:
    print('Elbow point found at n_clusters = ', elbow_point)
else:
    print('No elbow point found by KneeLocator, consider to change sensitivity.')
# plot the elbow curve
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(n_clusters, scores.loc['wcss'], zorder=0)
ax.vlines(elbow_point, ax.get_ylim()[0], ax.get_ylim()[1], ls='--', color='red', label='elbow location', zorder=1)
ax.set_xlabel('k-means clusters')
ax.set_ylabel('Within-cluster sum of squares')
ax.set_xticks(n_clusters)
ax.set_xticklabels(n_clusters)
ax.set_title('Elbow plot for 2-7 k-means cluster on the top 3 PCA drivers')
for p in n_clusters:
    ax.text(p, scores.loc['wcss', p], 'Silh.: ' + str("%.3f" % scores.loc['silhouette_avg', p]), zorder=2,
            ha='center', va='center', fontsize=8, bbox=dict(facecolor='w', edgecolor='k', boxstyle='round', alpha=0.5))
plt.tight_layout()
for form in formats:
    plt.savefig(f'./results_revised_review/comment_2/kmeans_elbow_plot.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

###################
# ## Figure 3 d) ##

# ## importing one of the files to get the cell lines names ##
data_temp = pd.read_excel("../data/circadian_parameters_by_replicate_Bmal1_final_v1.xlsx",
                          sheet_name="autocorrelation_lag")
cell_lines_names = list(data_temp.columns)[1:17]
cell_lines_names = [s.split("_")[0] for s in cell_lines_names]

# making the table using the function
mut_table23q2 = make_mut_table(cell_lines_names, "../data/OmicsSomaticMutations.csv")
mut_table = mut_table23q2.copy()
mut_table.columns = [k.replace('-', '').replace(' ', '') for k in mut_table.columns]

# merge missense_silent and missense
mut_table[mut_table.values == 'MISSENSE + SILENT'] = 'MISSENSE'

# remove U2OS cell line and RORA gene
mut_table = mut_table.drop('U2OS', axis=1)
mut_table = mut_table.drop('RORA', axis=0)

# add the other cell lines
missing_lines = [k for k in cell_lines_names if k not in mut_table.columns]
missing_lines.remove('MCF10A')
missing_lines.remove('U2OS')
mut_table[missing_lines] = np.nan

# transforming the table for output: mapping values to integers
unique_values = [np.nan, 'SILENT', 'MISSENSE', 'FRAME_SHIFT_DEL', 'NONSENSE']
value_to_num = {ni: indi for indi, ni in enumerate(unique_values, start=1)}
numerical_df = mut_table.replace(value_to_num).astype(int)
n = len(unique_values)

# Define the RGB values for the colors
light_grey = (211/255, 211/255, 211/255)
wheat = (253/255, 245/255, 230/255)
salmon = (250/255, 128/255, 114/255)
red = (210/255, 0, 0)
dark_red = (139/255, 0, 0)

# Create an array of RGBA values transitioning through the colors
vals = np.ones((256, 4))
vals[:64, 0] = np.linspace(light_grey[0], wheat[0], 64)  # Light grey to beige
vals[:64, 1] = np.linspace(light_grey[1], wheat[1], 64)
vals[:64, 2] = np.linspace(light_grey[2], wheat[2], 64)
vals[64:128, 0] = np.linspace(wheat[0], salmon[0], 64)  # Beige to salmon
vals[64:128, 1] = np.linspace(wheat[1], salmon[1], 64)
vals[64:128, 2] = np.linspace(wheat[2], salmon[2], 64)
vals[128:192, 0] = np.linspace(salmon[0], red[0], 64)  # Salmon to red
vals[128:192, 1] = np.linspace(salmon[1], red[1], 64)
vals[128:192, 2] = np.linspace(salmon[2], red[2], 64)
vals[192:, 0] = np.linspace(red[0], dark_red[0], 64)  # Red to dark red
vals[192:, 1] = np.linspace(red[1], dark_red[1], 64)
vals[192:, 2] = np.linspace(red[2], dark_red[2], 64)

# Create a ListedColormap from the RGBA values
my_cmp = ListedColormap(vals)

# plotting
fig = sns.clustermap(numerical_df, annot=False, linewidth=1, linecolor='w', cmap=my_cmp, cbar=False,
                     metric='hamming', method='single', figsize=(12, 6), dendrogram_ratio=(0, 0), cbar_pos=None,
                     row_cluster=False)
ax = plt.gca()
for xpos, cluster in zip([1.5, 6, 11.5], ['≥2 mutations', '1 mutation', 'no mutation']):
    ax.text(xpos, -0.2, cluster, ha='center', fontsize=14)
ax.plot([0.1, 2.9], [-0.1, -0.1], linewidth=2, color='k', clip_on=False)
ax.plot([3.1, 8.9], [-0.1, -0.1], linewidth=2, color='k', clip_on=False)
ax.plot([9.1, 13.9], [-0.1, -0.1], linewidth=2, color='k', clip_on=False)
ax.tick_params(axis=u'both', which=u'both', length=0)
# get the colors to add to our custom legend
im = fig.ax_heatmap.collections[0]
rgba_values = im.cmap(im.norm(im.get_array()))
cols, find_unique_val_location = np.unique(rgba_values, axis=0, return_index=True)
labels = im.get_array()[find_unique_val_location]

map_number_to_color = {unique_values[val-1]: col for (val, col) in zip(labels[labels-1], cols[labels-1])}
map_number_to_color = {key if not pd.isna(key) else 'NOT MUTATED': value for key, value in map_number_to_color.items()}

correct_order = ['NOT MUTATED', 'SILENT', 'MISSENSE', 'FRAME_SHIFT_DEL', 'NONSENSE']
ordered_color_map = {k: map_number_to_color[k] for k in correct_order}
dict_correct_name = {'NOT MUTATED': 'Not mutated', 'SILENT': 'Silent', 'MISSENSE': 'Missense',
                     'FRAME_SHIFT_DEL': 'Frameshift deletion', 'NONSENSE': 'Nonsense'}
# Create a list to hold the patch objects
patches = [Patch(facecolor=value, edgecolor='dimgrey',
                 label=dict_correct_name[label]) for label, value in ordered_color_map.items()]
# Add legend to the plot
plt.legend(handles=patches, title="Mutation", bbox_to_anchor=(1.01, 1.01), loc='upper left', alignment='left',
           title_fontproperties={'weight': 'bold'}, frameon=False)
plt.setp(fig.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
fig.ax_heatmap.axvline(3, color='white', lw=5)
fig.ax_heatmap.axvline(9, color='white', lw=5)
fig.ax_heatmap.yaxis.tick_left()
fig.ax_heatmap.yaxis.label_position = 'left'
fig.ax_col_dendrogram.set_title('Mutations of circadian clock genes across cell models')
plt.ylabel('Circadian clock genes', va='bottom')
plt.xlabel('')
plt.tight_layout()
for form in formats:
    plt.savefig(f'./results/3_d.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

####################
# ## Figure S4 c) ##

# ## REVISED

seed_here = 40

random.seed(seed_here)
np.random.seed(seed_here)

# S4 c circadian clock features
data_here = df_matrix_bmal1_per2

matrix = data_here.copy()
df = matrix.dropna(axis=1, how='all')
df = df.dropna(axis=0)
df_norm = (df - df.min()) / (df.max() - df.min())

targets = pd.DataFrame([x.split('_')[0] for x in df.index], index=df.index, columns=['subtype'])

# plot UMAP (neighbors=3)
plot_UMAP(df_norm, pd.concat([df_norm, targets], axis=1), 'subtype', n_neighbors=3, seed=seed_here,
          title='Circadian clock feature mapping', custom_color_marker=True, force_text=(2, 2))
for form in formats:
    plt.savefig(f'./results_revised/s4_c.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

# ## FINAL REVISION (Comment 2)

# - UMAP on all 9 raw features with n_neighbors 2, 3, 4 (+ stats if possible), currently only done for n=3

seed = 42

# S4 c circadian clock features
data_here = df_matrix_bmal1_per2

matrix = data_here.copy()
df = matrix.dropna(axis=1, how='all')
df = df.dropna(axis=0)
df_norm = (df - df.min()) / (df.max() - df.min())

targets = pd.DataFrame([x.split('_')[0] for x in df.index], index=df.index, columns=['subtype'])

n_neighbors = [2, 3, 4]

# plot UMAP (neighbors=2, 3, 4)
for n_neighbor in n_neighbors:
    random.seed(seed)
    np.random.seed(seed)
    plot_UMAP(df_norm, pd.concat([df_norm, targets], axis=1), 'subtype', n_neighbors=n_neighbor, seed=seed,
              title=f'Circadian clock feature mapping\n(UMAP n_neighbors={n_neighbor})',
              custom_color_marker=True, force_text=(2, 2))
    for form in formats:
        plt.savefig(f'./results_revised_review/comment_2/UMAP_all_nine_features_n_{n_neighbor}.{form}',
                    bbox_inches='tight', dpi=dpi)
    plt.close()

###########################
# ## Figure S4 d) and e) ##

# core clock gene expression across cell models
fig = sns.clustermap(X_short.T, cmap=plt.get_cmap('coolwarm'), figsize=(18, 10),
                     metric='euclidean', method='complete', dendrogram_ratio=(0.1, 0.1),
                     cbar_pos=(1.01, 0.25, 0.02, 0.425), annot=False, vmin=0, center=0.5, vmax=1, linewidth=.003)
fig.ax_heatmap.tick_params(axis=u'both', which=u'both', length=0)
fig.ax_heatmap.set_ylabel(None)
# add frame lines around heatmap
fig.ax_heatmap.axhline(y=0, color='k', linewidth=2)
fig.ax_heatmap.axhline(y=X_short.T.shape[0], color='k', linewidth=2)
fig.ax_heatmap.axvline(x=0, color='k', linewidth=2)
fig.ax_heatmap.axvline(x=X_short.T.shape[1], color='k', linewidth=2)
fig.ax_col_dendrogram.set_title("Core clock gene expression across cell models")
# adapt colorbar ticks
cbar = fig.ax_heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=8)
cbar.ax.set_yticks([0, 0.5, 1], [0, 0.5, 1])
cbar.ax.set_ylabel("Normalized expression", rotation=90)
cbar.ax.yaxis.set_label_position('left')
for form in formats:
    plt.savefig(f'./results/s4_d.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

# S4e)

seed_here = 40

random.seed(seed_here)
np.random.seed(seed_here)

X_short_this = X_short.copy()
X_short_this.index = [dict_to_map_line_and_subtype[x] + '_' + x for x in X_short.index]

targets = pd.DataFrame([x.split('_')[0] for x in X_short_this.index], index=X_short_this.index, columns=['subtype'])

# plot UMAP on ccle
plot_UMAP(X_short_this, pd.concat([X_short_this, targets], axis=1), 'subtype', n_neighbors=3, seed=seed_here,
          title="Circadian gene expression mapping", custom_color_marker=True, force_text=(0.5, 0.5))
for form in formats:
    plt.savefig(f'./results/s4_e.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

###################
# ## Figure 4 d) ##

# ## REVISED

random.seed(seed)
np.random.seed(seed)

drug = 'Cisplatin'
data_here = df_matrix_bmal1_per2.copy()
grinf_data = data_sensitivty['GRINF']
grinf_clean = grinf_data.drop('MCF10A', axis=1)

# transpose to get cisplatin data
drug_here = grinf_clean.T[drug]
this_drug_end = drug_here.dropna()
# binarize
this_drug_binarized = this_drug_end > this_drug_end.median()
# collect subtypes here
subtypes_here = pd.DataFrame([dict_to_map_line_and_subtype[x] for x in this_drug_binarized.index],
                             index=this_drug_binarized.index, columns=['subtype'])
# reshape circadian data
matrix_circ = data_here.dropna(axis=1, how='all')
matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
matrix_circ_norm_final = matrix_circ_norm.loc[subtypes_here.index]
# Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
# sample from the drug to keep consistent data sets
all_together = pd.concat([matrix_circ_norm_final, this_drug_binarized, subtypes_here], axis=1)

# plot LDA
plot_LDA(matrix_circ_norm_final, all_together, label=drug,
         title='$\\bf{Own\\ drug\\ sensitivity\\ dataset}$\nCisplatin | ' + '$GR_{inf}$ | ' + 'Circadian parameters',
         lda_output_target=drug, seed=seed, annot_subtype=True, annot_cellline=True, data_info=data_info,
         custom_color_marker=True, leg_labels=['$\\bf{Drug\\ effect}$', 'cytotoxic', 'cytostatic'], expand=(1.1, 1.5))

for form in formats:
    plt.savefig(f'./results_revised/4_d.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

# Export this LDA data and contributions
LDA = LinearDiscriminantAnalysis(n_components=1)
LDAs = LDA.fit_transform(matrix_circ_norm_final, all_together[drug])
LDAdf = pd.DataFrame(data=LDAs, columns=[f'LD{number + 1}' for number in np.arange(1)],
                     index=matrix_circ_norm_final.index)
LDAloadings = LDA.coef_
normalized_lda_loading = abs(LDAloadings[0]) / sum(abs(LDAloadings[0]))
sorted_normalized_idx = np.argsort(normalized_lda_loading)[::-1]
loadingsDF = pd.DataFrame(data=normalized_lda_loading[sorted_normalized_idx] * 100,
                          columns=['Distinctive contribution (%)'],
                          index=matrix_circ_norm_final.columns[sorted_normalized_idx])
LDAdf.to_excel('./results_revised_review/MSB-2024-12581_SourceData_Fig4C.xlsx', sheet_name='sheet1', index=True)
loadingsDF.to_excel('./results_revised_review/MSB-2024-12581_SourceData_Fig4D.xlsx', sheet_name='sheet1', index=True)

###################
# ## Figure 4 e) ##

# ## REVISED

random.seed(seed)
np.random.seed(seed)

drug = '17AAG'
data_here = df_matrix_bmal1_per2.copy()
ic50_data = ccle_drug_data[ccle_drug_data['Compound'] == drug]
ic50_clean = ic50_data[['Primary Cell Line Name', 'Compound', 'IC50 (uM)']]

# get this drug only
this_drug = ic50_clean.loc[ic50_clean['Compound'] == drug]
this_drug = this_drug[this_drug['Primary Cell Line Name'].isin([k.split('_')[1] for k in df_matrix_bmal1_per2.index])]
this_drug_end = this_drug.dropna()
this_drug_binarized = this_drug_end['IC50 (uM)'] > this_drug_end['IC50 (uM)'].median()
subtypes_here = pd.DataFrame([dict_to_map_line_and_subtype[x] for x in this_drug_end['Primary Cell Line Name']],
                             index=this_drug_end['Primary Cell Line Name'], columns=['subtype'])
this_drug_binarized.index = subtypes_here.index

# circadian data
matrix_circ = data_here.dropna(axis=1, how='all')
matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
matrix_circ_norm_final = matrix_circ_norm.loc[subtypes_here.index]

all_together = pd.concat([matrix_circ_norm_final, this_drug_binarized, subtypes_here], axis=1)

# plot LDA
plot_LDA(matrix_circ_norm_final, all_together, label='IC50 (uM)',
         title='$\\bf{Public\\ CCLE\\ dataset}$\n17-AAG | ' + '$IC_{50}$ | ' + 'Circadian parameters',
         lda_output_target='IC50 (uM)', seed=seed, annot_subtype=True, annot_cellline=True, data_info=data_info,
         custom_color_marker=True, leg_labels=['$\\bf{Sensitivity}$', 'high', 'low'], expand=(1.1, 1.5))

for form in formats:
    plt.savefig(f'./results_revised/4_e.{form}', bbox_inches='tight', dpi=dpi)
plt.close()
###################
# ## Figure 4 g) ##

# ## REVISED

random.seed(seed)
np.random.seed(seed)

# score dataframes
scores_lda_performance_circadian = pd.DataFrame(index=remaining_ccle_drugs['Compound'].unique(), columns=['IC50 (uM)'])
scores_lda_performance_ccle = pd.DataFrame(index=remaining_ccle_drugs['Compound'].unique(), columns=['IC50 (uM)'])

df = remaining_ccle_drugs.copy()
df = df[['Primary Cell Line Name', 'Compound', 'IC50 (uM)']]

# collect the bcd, wcd, and bcd-wcd ratio in bulk for both circadian parameters and gene expression sets, targeting IC50
for drug in remaining_drugs:
    print(f'\n* Analysing for drug {drug}...')
    # get this drug only
    this_drug = df.loc[df['Compound'] == drug]
    this_drug_end = this_drug.dropna()
    # the number of samples must be more than the number of classes, so we always need at least 3 samples
    if len(this_drug_end) <= 2:
        # somtimes, IC50 has mostly the value 8, resulting in all classes being True / False
        text = f'only {len(this_drug_end)} cell line was' if len(this_drug_end) == 1 else \
            f'only {len(this_drug_end)} cell line were' if len(this_drug_end) > 1 else 'no cell lines were'
        print(f'For the drug {drug}, {text} measured for IC50 (uM). This example is skipped for the LDA analysis.')
        continue  # skip if we only remain with 1 or 0 samples (most cases are EC50 (uM) not being measured)
    print(f'** Discretize drug sensitivity by median...')
    this_drug_binarized = this_drug_end['IC50 (uM)'] > this_drug_end['IC50 (uM)'].median()
    if len(this_drug_binarized.unique()) == 1:  # apply the -0.001 trick only if necessary!
        this_drug_binarized = this_drug_end['IC50 (uM)'] > this_drug_end['IC50 (uM)'].median() - 0.001
        print('## Trying "median - 0.001" for discretiztation cutoff ##')
    if len(this_drug_binarized.unique()) == 1:
        print(f'For the drug {drug} measured for IC50 (uM), the binarization by median resulted in one '
              f'single class. This example is skipped for the LDA analysis.')
        continue
    subtypes_here = pd.DataFrame(
        [dict_to_map_line_and_subtype[x] for x in this_drug_end['Primary Cell Line Name']],
        index=this_drug_end['Primary Cell Line Name'], columns=['subtype'])
    this_drug_binarized.index = subtypes_here.index
    # •	By circadian clock gene expression levels
    print('** LDA on CCLE expression...')  # without 2UOS
    ccle_expression = X_short.loc[this_drug_end['Primary Cell Line Name']]
    # LDA machinery
    LDA = LinearDiscriminantAnalysis(n_components=1)
    LDAs = LDA.fit_transform(ccle_expression, this_drug_binarized)
    LDAdf = pd.DataFrame(data=LDAs, columns=['LD1'])
    LDAloadings = LDA.coef_
    normalized_lda_loading = abs(LDAloadings[0]) / sum(abs(LDAloadings[0]))
    sorted_normalized_idx = np.argsort(normalized_lda_loading)[::-1]
    # collect scores here in database1 and excel outputs
    LDAdf.index = this_drug_binarized.index
    loadings_to_export = pd.DataFrame([normalized_lda_loading[sorted_normalized_idx],
                                       LDAloadings.flatten()[sorted_normalized_idx]]).T
    loadings_to_export.columns = ['Loadings (normalized)', 'Loadings (raw)']
    loadings_to_export.index = ccle_expression.columns[sorted_normalized_idx]
    with pd.ExcelWriter(f'./results_revised/4g_gene_expression_{drug}.xlsx') as writer:
        LDAdf.to_excel(writer, sheet_name='Discriminants')
        loadings_to_export.to_excel(writer, sheet_name='Loadings')
    scores_lda_performance_ccle['IC50 (uM)'][drug] = get_LDA_metrics(classes=this_drug_binarized, LDs=LDAdf)
    # •	By circadian features mean
    mat = df_matrix_bmal1_per2.copy()
    print(f'** LDA on BMAL1-PER2 parameters...')
    matrix_circ = mat.dropna(axis=1, how='all')
    matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
    matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
    matrix_circ_norm = matrix_circ_norm.loc[subtypes_here.index]
    # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
    # sample from the drug to keep consistent data sets
    all_together = pd.concat([matrix_circ_norm, this_drug_binarized, subtypes_here], axis=1)
    all_together_final = all_together.dropna(axis=0)
    matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
    this_drug_binarized_consistent = all_together_final['IC50 (uM)']
    if this_drug_binarized_consistent.nunique() == 1 or len(all_together_final) <= 2:
        print(f'/!\\ LDA on circadian parameters targeting {drug} is only possible '
              f'for\n{len(all_together_final)} samples where parameters and IC50 (uM) are available,'
              f'\nbut the binarization of these remaining samples is left with 1 class only using median.'
              f'\nLDA requires at least 2 classes and 3 samples. /!\\')
        continue
    # LDA machinery
    LDA = LinearDiscriminantAnalysis(n_components=1)
    LDAs = LDA.fit_transform(matrix_circ_norm_final, all_together_final['IC50 (uM)'])
    LDAdf = pd.DataFrame(data=LDAs, columns=['LD1'])
    LDAloadings = LDA.coef_
    normalized_lda_loading = abs(LDAloadings[0]) / sum(abs(LDAloadings[0]))
    sorted_normalized_idx = np.argsort(normalized_lda_loading)[::-1]
    # collect scores here in database2 and excel outputs
    LDAdf.index = this_drug_binarized_consistent.index
    loadings_to_export = pd.DataFrame(
        [normalized_lda_loading[sorted_normalized_idx], LDAloadings.flatten()[sorted_normalized_idx]]).T
    loadings_to_export.columns = ['Loadings (normalized)', 'Loadings (raw)']
    loadings_to_export.index = matrix_circ_norm_final.columns[sorted_normalized_idx]
    with pd.ExcelWriter(f'./results_revised/4g_circadian_parameters_{drug}.xlsx') as writer:
        LDAdf.to_excel(writer, sheet_name='Discriminants')
        loadings_to_export.to_excel(writer, sheet_name='Loadings')
    scores_lda_performance_circadian['IC50 (uM)'][drug] = get_LDA_metrics(classes=this_drug_binarized_consistent,
                                                                          LDs=LDAdf)
# NaNs for Nutlin3, LBW242, PD0332991, AZD6244, and PLX4720

# circadian results
scores_by_meas_circadian = scores_lda_performance_circadian['IC50 (uM)'].dropna()
bcd_circadian = scores_by_meas_circadian.apply(lambda x: x[0])
wcd_circadian = scores_by_meas_circadian.apply(lambda x: x[1])
ratios_circadian = scores_by_meas_circadian.apply(lambda x: x[-1])

# ccle results
scores_by_meas_ccle = scores_lda_performance_ccle['IC50 (uM)'].dropna()
bcd_ccle = scores_by_meas_ccle.apply(lambda x: x[0])
wcd_ccle = scores_by_meas_ccle.apply(lambda x: x[1])
ratios_ccle = scores_by_meas_ccle.apply(lambda x: x[-1])

# custom_legend
legend_elements = []
for mark, typ in zip(['x', '.', '--'],
                     ['Between-Cluster-Distance', 'Within-Cluster-Distance', 'Discrimination threshold']):
    if mark != '--':
        legend_elements.append(Line2D([0], [0], marker=mark, color='w', markerfacecolor='none' if mark != 'x' else 'k',
                                      markeredgecolor='k', markersize=8, label=typ))
    else:
        legend_elements.append(Line2D([0], [0], ls=mark, lw=1, color='k', label=typ))

# plot both
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), layout='constrained')
df_ax1 = ratios_circadian.sort_values(ascending=False)
sns.barplot(x=df_ax1.index, y=df_ax1.values, color='lightblue', edgecolor='k', ax=ax1)
# add bcd as black crosses and wcd as black dots (sorted by the ratio!)
sns.scatterplot(x=df_ax1.index, y=bcd_circadian[df_ax1.index].values, color='k', marker='x', s=25, ax=ax1)
sns.scatterplot(x=df_ax1.index, y=wcd_circadian[df_ax1.index].values, facecolor='none', marker='.', s=50,
                edgecolor='k', ax=ax1)
ax1.hlines(y=2, xmin=0-1, xmax=len(df_ax1.index), color='k', linestyle='--')
ax1.set_xticks(ticks=np.arange(len(df_ax1.index)),
               labels=['17-AAG' if x == '17AAG' else x for x in df_ax1.index], rotation=90, ha='center')
ax1.set_xlim(-1, len(df_ax1.index))
ax1.set_xlabel('')
ax1.set_ylabel('Chronosensitivity index')
ax1.set_title('Dependency on circadian parameters')
df_ax2 = ratios_ccle.sort_values(ascending=False)
sns.barplot(x=df_ax2.index, y=df_ax2.values, color='lightblue', edgecolor='k', ax=ax2)
# add bcd as black crosses and wcd as black dots (sorted by the ratio!)
sns.scatterplot(x=df_ax2.index, y=bcd_ccle[df_ax2.index].values, color='k', marker='x', s=25, ax=ax2)
sns.scatterplot(x=df_ax2.index, y=wcd_ccle[df_ax2.index].values, facecolor='none', marker='.', s=50,
                edgecolor='k', ax=ax2)
ax2.hlines(y=2, xmin=0-1, xmax=len(df_ax2.index), color='k', linestyle='--')
ax2.set_xticks(ticks=np.arange(len(df_ax2.index)),
               labels=['17-AAG' if x == '17AAG' else x for x in df_ax2.index], rotation=90, ha='center')
ax2.set_xlim(-1, len(df_ax2.index))
ax2.set_xlabel('')
ax2.set_ylabel('')
ax2.set_title('Dependency on circadian gene expression')
ax2.legend(handles=legend_elements)
for ax in (ax1, ax2):
    triangle = Polygon([[0.03, -0.45], [0.97, -0.45], [0.03, -0.35]], facecolor='lightblue',
                       transform=ax.transAxes, clip_on=False, closed=True)
    ax.add_patch(triangle)
    ax.text(0.5, -0.5, 'Circadian clock dependency', va='center', ha='center', transform=ax.transAxes)
fig.suptitle('$\\bf{Public\\ CCLE\\ dataset}$\nChronosensitivity indices across drugs | ' + '$IC_{50}$')
for txt, xpos, ypos in zip(['> 2 good discrimination', '< 2 poor discrimination'], [16, 16], [2.2, 1.8]):
    ax2.annotate(txt, xy=(xpos, ypos), va='center', ha='center')

for form in formats:
    plt.savefig(f'./results_revised/4_g.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

# Export data into xlsx
ratios_ccle.to_excel('./results_revised_review/MSB-2024-12581_SourceData_Fig4H.xlsx', sheet_name='sheet1', index=True)
ratios_circadian.to_excel('./results_revised_review/MSB-2024-12581_SourceData_Fig4I.xlsx', sheet_name='sheet1', index=True)

########################################
# ## Figure S6 a), b), c), d), and e) ##

# ## REVISED

random.seed(seed)
np.random.seed(seed)

drug = 'Cisplatin'
this_drug_data = data_sensitivty

elems = ['s6_a', 's6_b', 's6_c', 's6_d', 's6_e']

for fig_elem, meas, data_here in zip(elems,
                                     ['GRINF', 'GR50', 'GR50', 'Hill', 'Hill'],
                                     [X_short, df_matrix_bmal1_per2, X_short, df_matrix_bmal1_per2, X_short]):
    # get the correct binarized target
    sensitivity = this_drug_data.copy()
    meas_data = sensitivity[meas]
    meas_clean = meas_data.drop('MCF10A', axis=1)
    # transpose to get cisplatin data
    drug_here = meas_clean.T[drug]
    this_drug_end = drug_here.dropna()
    # binarize
    this_drug_binarized = this_drug_end > this_drug_end.median()
    # collect subtypes here
    subtypes_here = pd.DataFrame([dict_to_map_line_and_subtype[x] for x in this_drug_binarized.index],
                                 index=this_drug_binarized.index, columns=['subtype'])
    this_drug_binarized.index = subtypes_here.index
    # put the input data in shape
    df = data_here.copy()
    if fig_elem in ('s6_b', 's6_d'):  # in case of circadian parameters
        matrix_circ = df.dropna(axis=1, how='all')
        matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
        matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
        input_final = matrix_circ_norm.loc[subtypes_here.index]
    else:
        input_final = df.loc[this_drug_end.index]
    # combine input data with subtypes and drug sensitivity for plotting
    all_together = pd.concat([input_final, this_drug_binarized, subtypes_here], axis=1)
    all_together_final = all_together.dropna(axis=0)
    input_final = input_final.dropna(axis=0)
    # plot LDA
    sub_meas = '$GR_{inf}$' if meas == 'GRINF' else '$GR_{50}$' if meas == 'GR50' else 'Hill coefficient'
    plot_LDA(input_final, all_together_final, label=drug,
             title=f'{sub_meas} | Circadian parameters' if fig_elem in ('s6_b', 's6_d') else
             f'{sub_meas} | Circadian gene expression',
             lda_output_target=drug, seed=seed, annot_subtype=True, annot_cellline=True, data_info=data_info,
             custom_color_marker=True,
             leg_labels=['$\\bf{Drug\\ effect}$', 'cytotoxic', 'cytostatic'] if fig_elem == 's6_a' else
             ['$\\bf{Sensitivity}$', 'high', 'low'], expand=(1.1, 1.5))
    # save each of them
    for form in formats:
        plt.savefig(f'./results_revised/{fig_elem}.{form}', bbox_inches='tight', dpi=dpi)
    plt.close()

####################
# ## Figure S6 f) ##

# ## REVISED

random.seed(seed)
np.random.seed(seed)

measure = 'IC50 (uM)'
drugs = ['17AAG', 'AEW541', 'Paclitaxel', 'Topotecan']
data_here = df_matrix_bmal1_per2.copy()
this_drug_data = remaining_ccle_drugs.copy()
df = this_drug_data[['Primary Cell Line Name', 'Compound', measure]]

matrix_circ = data_here.dropna(axis=1, how='all')
matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
input_final = matrix_circ_norm.loc[subtypes_here.index]

fig, axes = plt.subplots(1, 4, figsize=(16, 4))
txt_by_axis = [[], [], [], []]
for num, (drug, ax) in enumerate(zip(drugs, axes.ravel())):
    # get this drug only
    this_drug = df.loc[df['Compound'] == drug]
    this_drug_end = this_drug.dropna()
    this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median()
    subtypes_here = pd.DataFrame([dict_to_map_line_and_subtype[x] for x in this_drug_end['Primary Cell Line Name']],
                                 index=this_drug_end['Primary Cell Line Name'], columns=['subtype'])
    this_drug_binarized.index = subtypes_here.index
    LDA = LinearDiscriminantAnalysis(n_components=1)
    LDAs = LDA.fit_transform(input_final, this_drug_binarized)
    LDAdf = pd.DataFrame(data=LDAs, columns=['LD1'])
    np.random.seed(seed)
    jittered_y = pd.DataFrame(np.zeros(len(LDAdf)) + 0.1 * np.random.rand(len(LDAdf)) - 0.05, columns=['LD1'])
    ax.scatter(LDAdf['LD1'], jittered_y['LD1'],
               color=['teal' if not k else 'orange' for k in this_drug_binarized.values], s=50, alpha=0.75, marker='o')
    for name, x, y in zip(this_drug_binarized.index, LDAdf['LD1'], jittered_y['LD1']):
        txt_by_axis[num].append(ax.text(x, y, name, fontsize=8, fontweight='bold', ha='center', va='center'))
    _, _, ratio = get_LDA_metrics(this_drug_binarized, LDAdf)
    sign_and_ratio = '= %.1f' % ratio if ratio >= 2 and drug == '17AAG' else \
        '= %.0f' % ratio if ratio >= 2 and drug != '17AAG' else '< 2'
    ax.set_title(f'{"17-AGG" if drug == "17AAG" else drug} | Chronosensitivity index {sign_and_ratio}', fontsize=10)
    ax.set_xlabel('LD1')
    ax.set_ylabel('')
    ax.set_ylim([-0.1, 0.1])
    ax.set_yticks([])
    if num == 0:
        # legend in first plot
        legend_elements = [plt.plot([], marker="", ls="")[0]]
        for c, typ in zip(['teal', 'orange'], ['high', 'low']):
            legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=8, label=typ))
        leg = ax.legend(handles=legend_elements, labels=['$\\bf{Sensitivity}$', 'high', 'low'],
                        loc='lower left', ncol=3, fontsize=9, frameon=True, columnspacing=0.3, handletextpad=0.1)
        for vpack in leg._legend_handle_box.get_children()[:1]:
            for hpack in vpack.get_children():
                hpack.get_children()[0].set_width(0)
fig.suptitle('$\\bf{Public\\ CCLE\\ dataset}$\nCircadian parameters | ' + '$IC_{50}$')
fig.tight_layout()
for txt, ax in zip(txt_by_axis, axes.ravel()):
    adjust_text(txt, ax=ax, expand=(1.1, 1.5), arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
# save figure
for form in formats:
    plt.savefig(f'./results_revised/s6_f.{form}', bbox_inches='tight', dpi=dpi)
plt.close()

###########################
# ## Additional analysis ##

# ## FINAL REVISION (Comment 5 --> 6, 7)

# - comment 5: LOOCV LDA on data of manuscript (own cisplatin)
# - comment 5: include GRAOC and GEC50 for cisplatin, add LOOCV

drug = 'Cisplatin'
tod_scores_for_loocv = list(data_sensitivty.keys())

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                       index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_circadian = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                            index=metrics, columns=tod_scores_for_loocv)

# Without sample MCF10A
for measure in list(data_sensitivty.keys()):
    print(f'* Analysing stat {measure}...')
    df = data_sensitivty[measure].T.copy()
    df = df.drop('MCF10A', axis=0)
    this_drug = df[drug]
    this_drug_end = this_drug.dropna()
    print(f'** Discretize drug sensitivity by median...')
    this_drug_binarized = this_drug_end > this_drug_end.median()
    # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
    print('**** LDA on CCLE expression...')
    ccle_expression = ccle_data.loc[this_drug_end.index]
    # LDA LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LDA_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                  target=drug)
    loocv_results_grid_ccle[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # •	By circadian features mean
    print(f'**** LDA on BMAL1_PER2 mean parameters...\n')
    mat = circadian_data.copy()
    matrix_circ = mat.dropna(axis=1, how='all')
    matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
    matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
    matrix_circ_norm = matrix_circ_norm.loc[this_drug_end.index]
    # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
    # sample from the drug to keep consistent data sets
    all_together = pd.concat([matrix_circ_norm, this_drug_binarized], axis=1)
    all_together_final = all_together.dropna(axis=0)
    matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
    this_drug_binarized_consistent = all_together_final[drug]
    # LDA LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LDA_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug)
    loocv_results_grid_circadian[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]

# plot at this level for all measures (target) and only Cisplatin, 2 per measure (CCLE and circ (including channels))

# Plot for CCLE input
loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
columns = set()
for col in loocv_samples_ccle.columns:
    for index_object in loocv_samples_ccle[col]:
        if isinstance(index_object, pd.Index):
            columns.update(index_object)
columns = list(columns)
indeces = loocv_samples_ccle.index
this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
for drug in indeces:
    for sample in columns:
        # case of correct samples and wrong samples
        if sample in loocv_samples_ccle.loc[drug]['corr_samples'].__str__():
            this_table_ccle.loc[drug][sample] += 1
        # case of the N/A
        if sample not in loocv_samples_ccle.loc[drug]['corr_samples'].__str__() and sample not in \
                loocv_samples_ccle.loc[drug]['not_corr_samples'].__str__():
            this_table_ccle.loc[drug][sample] = np.nan
# add accuracy
this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
sorted_table = sorted_table.dropna(axis=0, how='all')
# table figure
# Create a mask for the 10th column and nans
nanmask = sorted_table.isna()
mask = np.zeros_like(sorted_table, dtype=bool)
mask[:, -1] = ~nanmask.iloc[:, -1]
combinedmask = mask | nanmask
fig = sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                  annot=False, vmin=0, center=.5, vmax=1, linewidth=.003)
# set the np.nan tiles to grey
fig.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
           extent=fig.get_xlim() + fig.get_ylim(), alpha=0.5)
for i in range(sorted_table.shape[0]):
    fig.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
             ha='center', va='center', color='black', fontsize=8)
samples_per_drug = []
for tod_drug in sorted_table.index:
    samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
# set the labels plus sample size
fig.set_yticklabels(y_labels, rotation=0, fontsize=12)
# add frame lines around heatmap
fig.axhline(y=0, color='k', linewidth=2)
fig.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
fig.axvline(x=0, color='k', linewidth=2)
fig.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
fig.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
fig.xaxis.tick_top()  # x axis on top
fig.xaxis.set_label_position('top')
fig.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
fig.set_title(f"LDA leave-one-out classification results by cell line, Cisplatin (CCLE)\nsample size in parenthesis")
for file in formats:
    plt.savefig(f'./results_revised_review/comment_5/CCLE_in_Cisplatin_out/lda/CCLE_Cisplatin_loocv_'
                f'classification_plot_lda.{file}', bbox_inches='tight', dpi=300)
plt.close()

# Plot for Circadian input
loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
columns = set()
for col in loocv_samples_circadian.columns:
    for index_object in loocv_samples_circadian[col]:
        if isinstance(index_object, pd.Index):
            columns.update(index_object)
columns = list(columns)
indeces = loocv_samples_circadian.index
this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
for drug in indeces:
    for sample in columns:
        # case of correct samples and wrong samples
        if sample in loocv_samples_circadian.loc[drug]['corr_samples'].__str__():
            this_table_circadian.loc[drug][sample] += 1
        # case of the N/A
        if sample not in loocv_samples_circadian.loc[drug]['corr_samples'].__str__() and sample not in \
                loocv_samples_circadian.loc[drug]['not_corr_samples'].__str__():
            this_table_circadian.loc[drug][sample] = np.nan
# add accuracy
this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
sorted_table = sorted_table.dropna(axis=0, how='all')
# table figure
# Create a mask for the 10th column and nans
nanmask = sorted_table.isna()
mask = np.zeros_like(sorted_table, dtype=bool)
mask[:, -1] = ~nanmask.iloc[:, -1]
combinedmask = mask | nanmask
fig = sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                  annot=False, vmin=0, center=.5, vmax=1, linewidth=.003)
# set the np.nan tiles to grey
fig.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
           extent=fig.get_xlim() + fig.get_ylim(), alpha=0.5)
for i in range(sorted_table.shape[0]):
    fig.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
             ha='center', va='center', color='black', fontsize=8)
samples_per_drug = []
for tod_drug in sorted_table.index:
    samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
# set the labels plus sample size
fig.set_yticklabels(y_labels, rotation=0, fontsize=12)
# add frame lines around heatmap
fig.axhline(y=0, color='k', linewidth=2)
fig.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
fig.axvline(x=0, color='k', linewidth=2)
fig.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
fig.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
fig.xaxis.tick_top()  # x axis on top
fig.xaxis.set_label_position('top')
fig.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
fig.set_title(f"LDA leave-one-out classification results by cell line, Cisplatin (BMAL1-PER2 mean)\n"
              f"sample size in parenthesis")
for file in formats:
    plt.savefig(f'./results_revised_review/comment_5/Circadian_in_Cisplatin_out/lda/Circadian_Cisplatin_BMAL1_PER2_'
                f'loocv_classification_plot_lda.{file}', bbox_inches='tight', dpi=300)
plt.close()

# - comment 5: LOOCV LDA on data of manuscript (CCLE drug screening)
# - comment 5: include AUC and EC50 for CCLE, select same drug as in IC50, add LOOCV

drugs = remaining_drugs.copy()
tod_scores_for_loocv = ['EC50 (uM)', 'IC50 (uM)', 'ActArea']

# limit EC50 and ActArea drugs to the ones where IC50 is at least available
only_ic50_drugs = [d for d in remaining_ccle_drugs['Compound'].unique() if len(
    remaining_ccle_drugs[remaining_ccle_drugs['Compound'] == d]['IC50 (uM)']) > 3 and
                   remaining_ccle_drugs[remaining_ccle_drugs['Compound'] == d]['IC50 (uM)'].nunique() > 2]

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)
loocv_results_grid_circadian = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)

# For each drug
for measure in tod_scores_for_loocv:
    print(f'* Analysing stat {measure}...')
    df = remaining_ccle_drugs.copy()
    df = df[['Primary Cell Line Name', 'Compound', measure]]
    for drug in [x for x in remaining_drugs if x in only_ic50_drugs]:
        print(f'*** Analysing for drug {drug}...')
        # get this drug only
        this_drug = df.loc[df['Compound'] == drug]
        this_drug_end = this_drug.dropna()
        # the number of samples must be more than the number of classes, so we always need at least 3 samples
        if len(this_drug_end) <= 3:
            # somtimes, IC50 has mostly the value 8, resulting in all classes being True / False
            text = f'only {len(this_drug_end)} cell line was' if len(this_drug_end) == 1 else \
                f'only {len(this_drug_end)} cell line were' if len(this_drug_end) > 1 else 'no cell lines were'
            print(f'For the drug {drug}, {text} measured for {measure}. This example is skipped for the LDA analysis.')
            continue  # skip if we only remain with 1 or 0 samples (most cases are EC50 (uM) not being measured)
        print(f'** Discretize drug sensitivity by median...')
        this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median()
        if len(this_drug_binarized.unique()) == 1:  # apply the -0.001 trick only if necessary!
            this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median() - 0.001
            print('## Trying "median - 0.001" for discretiztation cutoff ##')
        if len(this_drug_binarized.unique()) == 1:
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in one '
                  f'single class. This example is skipped for the LDA analysis.')
            continue
        if 1 in set(this_drug_binarized.value_counts()):  # skip LDA loocv of samples where ther is only 1 true or false
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in two classes but '
                  f'one is only represented once. This example is skipped for the LDA LOOCV analysis.')
            continue
        subtypes_here = pd.DataFrame(
            [dict_to_map_line_and_subtype[x] for x in this_drug_end['Primary Cell Line Name']],
            index=this_drug_end['Primary Cell Line Name'], columns=['subtype'])
        this_drug_binarized.index = subtypes_here.index
        # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
        print('**** LDA on CCLE expression...')
        ccle_expression = ccle_data.loc[this_drug_end['Primary Cell Line Name']]
        # LDA LOOCV machinery
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LDA_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                      target=measure)
        loocv_results_grid_ccle[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
        # •	By circadian features mean
        print(f'**** LDA on BMAL1_PER2 mean parameters...\n')
        mat = circadian_data.copy()
        matrix_circ = mat.dropna(axis=1, how='all')
        matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
        matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
        matrix_circ_norm = matrix_circ_norm.loc[subtypes_here.index]
        # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
        # sample from the drug to keep consistent data sets
        all_together = pd.concat([matrix_circ_norm, this_drug_binarized, subtypes_here], axis=1)
        all_together_final = all_together.dropna(axis=0)
        matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
        this_drug_binarized_consistent = all_together_final[measure]
        if 1 in this_drug_binarized_consistent.value_counts().values or len(
                all_together_final) <= 3:  # also here <=3 instead of 2
            print(f'/!\\ LDA on circadian parameters targeting {drug} is only possible '
                  f'for\n{len(all_together_final)} samples where parameters and {measure} are available,'
                  f'\nbut the binarization of these remaining samples is left with 1 class only using median.'
                  f'\nLDA requires at least 2 classes and 3 samples. /!\\')
            continue
        # LDA LOOCV machinery
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LDA_loocv(data=matrix_circ_norm_final, y=all_together_final, target=measure)
        loocv_results_grid_circadian[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # plot at this level for each individual measure (target), 2 per measure (CCLE and circ (including channels))
    # Plot for CCLE input
    loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_ccle.columns:
        for index_object in loocv_samples_ccle[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_ccle.index
    this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for drug in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_ccle.loc[drug]['corr_samples'].__str__():
                this_table_ccle.loc[drug][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_ccle.loc[drug]['corr_samples'].__str__() and sample not in \
                    loocv_samples_ccle.loc[drug][
                        'not_corr_samples'].__str__():
                this_table_ccle.loc[drug][sample] = np.nan
    # add accuracy
    this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
    sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig = sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                      annot=False, vmin=0, center=.5, vmax=1, linewidth=.003)
    # set the np.nan tiles to grey
    fig.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=fig.get_xlim() + fig.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        fig.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    fig.set_yticklabels(y_labels, fontsize=12)
    # add frame lines around heatmap
    fig.axhline(y=0, color='k', linewidth=2)
    fig.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    fig.axvline(x=0, color='k', linewidth=2)
    fig.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    fig.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    fig.xaxis.tick_top()  # x axis on top
    fig.xaxis.set_label_position('top')
    fig.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    txt = ', drugs limited to IC50 availability' if measure != 'IC50 (uM)' else ''
    fig.set_title(f"LDA leave-one-out classification results by cell line, {measure} (CCLE)\n"
                  f"sample size in parenthesis" + txt)
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/CCLE_in_CCLE_out/lda/CCLE_{measure}_'
                    f'loocv_classification_plot_lda.{file}', bbox_inches='tight', dpi=300)
    plt.close()
    # Plot for Circadian input
    loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_circadian.columns:
        for index_object in loocv_samples_circadian[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_circadian.index
    this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for drug in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_circadian.loc[drug]['corr_samples'].__str__():
                this_table_circadian.loc[drug][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_circadian.loc[drug]['corr_samples'].__str__() and sample not in \
                    loocv_samples_circadian.loc[drug][
                        'not_corr_samples'].__str__():
                this_table_circadian.loc[drug][sample] = np.nan
    # add accuracy
    this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
    sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig = sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                      annot=False, vmin=0, center=.5, vmax=1, linewidth=.003)
    # set the np.nan tiles to grey
    fig.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=fig.get_xlim() + fig.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        fig.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    fig.set_yticklabels(y_labels, fontsize=12)
    # add frame lines around heatmap
    fig.axhline(y=0, color='k', linewidth=2)
    fig.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    fig.axvline(x=0, color='k', linewidth=2)
    fig.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    fig.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    fig.xaxis.tick_top()  # x axis on top
    fig.xaxis.set_label_position('top')
    fig.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    txt = ', drugs limited to IC50 availability' if measure != 'IC50 (uM)' else ''
    fig.set_title(f"LDA leave-one-out classification results by cell line on {measure} (BMAL1-PER2 mean)\n"
                  f"sample size in parenthesis" + txt)
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/Circadian_in_CCLE_out/lda/circadian_{measure}_BMAL1_PER2_'
                    f'loocv_classification_plot_lda.{file}', bbox_inches='tight', dpi=300)
    plt.close()

# - comment 5: Logistic regression on binarized data, extract staitstical values, and perform LOOCV (Cisplatin own data)

drug = 'Cisplatin'
tod_scores_for_loocv = list(data_sensitivty.keys())

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                       index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_circadian = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                            index=metrics, columns=tod_scores_for_loocv)

original_hatch_width = 1.0
matplotlib.rcParams['hatch.linewidth'] = 8.0

# Without sample MCF10A
for measure in list(data_sensitivty.keys()):
    print(f'* Analysing stat {measure}...')
    df = data_sensitivty[measure].T.copy()
    df = df.drop('MCF10A', axis=0)
    this_drug = df[drug]
    this_drug_end = this_drug.dropna()
    print(f'** Discretize drug sensitivity by median...')
    this_drug_binarized = this_drug_end > this_drug_end.median()
    # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
    print('**** Logistic regression on CCLE expression...')
    ccle_expression = ccle_data.loc[this_drug_end.index]
    # Logreg on full data
    logreg = LogisticRegression(random_state=seed).fit(ccle_expression, this_drug_binarized)
    probs = logreg.predict_proba(ccle_expression)
    # Extract probabilities for the positive class
    positive_probs = probs[:, 1]
    # Extract feature importance
    importance = logreg.coef_[0]
    sorted_indices = np.argsort(np.abs(importance))[::-1]  # Sort by absolute importance
    top_importance = importance[sorted_indices]
    top_feature_names = np.array(ccle_expression.columns)[sorted_indices]
    # Plot probability distributions and feature importance
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.histplot(positive_probs[this_drug_binarized == 1], color="green", label="Positive Class", kde=True, bins='fd',
                 ax=axes[0])
    sns.histplot(positive_probs[this_drug_binarized == 0], color="red", label="Negative Class", kde=True, bins='fd',
                 ax=axes[0])
    axes[0].set_title(f"Logsitic regression training performance, Cisplatin {measure} (CCLE)")
    axes[0].set_xlabel("Predicted Probability")
    axes[0].set_ylabel("Frequency")
    axes[0].legend()
    # Feature importance
    axes[1].bar(range(0, ccle_expression.shape[1]), np.abs(top_importance), color='lavender', edgecolor='k',
                linewidth=1,  hatch=['\\' if i < 0 else None for i in top_importance])
    for i in range(0, ccle_expression.shape[1]):
        axes[1].annotate(str("{:.2f}".format(round(top_importance[i], 2))),
                         xy=(i, np.abs(top_importance[i])), ha='center', va='bottom', size=8, weight='normal')
    axes[1].set_title('Logistic regression feature weights', fontsize=10)
    axes[1].set_xlabel(None)
    axes[1].set_ylabel('Weights')
    axes[1].set_xticks(range(0, ccle_expression.shape[1]))
    axes[1].set_xticklabels(top_feature_names, rotation=90)
    axes[1].legend([Patch(facecolor='lavender', edgecolor='k'),
                    Patch(facecolor='lavender', edgecolor='k', hatch='\\')], ['Positive', 'Negative'])
    plt.tight_layout()
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/CCLE_in_Cisplatin_out/logreg/CCLE_Cisplatin_{measure}_logreg.'
                    f'{file}', bbox_inches='tight', dpi=300)
    plt.close()
    # LogReg LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogReg_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                     target=drug)
    loocv_results_grid_ccle[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # •	By circadian features mean
    print(f'**** Logistic regression on BMAL1_PER2 mean parameters...\n')
    mat = circadian_data.copy()
    matrix_circ = mat.dropna(axis=1, how='all')
    matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
    matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
    matrix_circ_norm = matrix_circ_norm.loc[this_drug_end.index]
    # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
    # sample from the drug to keep consistent data sets
    all_together = pd.concat([matrix_circ_norm, this_drug_binarized], axis=1)
    all_together_final = all_together.dropna(axis=0)
    matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
    this_drug_binarized_consistent = all_together_final[drug]
    # Logreg on full data
    logreg = LogisticRegression(random_state=seed).fit(matrix_circ_norm_final, this_drug_binarized_consistent)
    probs = logreg.predict_proba(matrix_circ_norm_final)
    # Extract probabilities for the positive class
    positive_probs = probs[:, 1]
    # Extract feature importance
    importance = logreg.coef_[0]
    sorted_indices = np.argsort(np.abs(importance))[::-1]  # Sort by absolute importance
    top_importance = importance[sorted_indices]
    top_feature_names = np.array(matrix_circ_norm_final.columns)[sorted_indices]
    # Plot probability distributions and feature importance
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.histplot(positive_probs[this_drug_binarized_consistent == 1], bins='fd',
                 color="green", label="Positive Class", kde=True, ax=axes[0])
    sns.histplot(positive_probs[this_drug_binarized_consistent == 0], bins='fd',
                 color="red", label="Negative Class", kde=True, ax=axes[0])
    axes[0].set_title(f"Logsitic regression training performance, Cisplatin {measure} (BMAL1-PER2 mean)")
    axes[0].set_xlabel("Predicted Probability")
    axes[0].set_ylabel("Frequency")
    axes[0].legend()
    # Feature importance
    axes[1].bar(range(0, matrix_circ_norm_final.shape[1]), np.abs(top_importance), color='lavender', edgecolor='k',
                linewidth=1, hatch=['\\' if i < 0 else None for i in top_importance])
    for i in range(0, matrix_circ_norm_final.shape[1]):
        axes[1].annotate(str("{:.2f}".format(round(top_importance[i], 2))),
                         xy=(i, np.abs(top_importance[i])), ha='center', va='bottom', size=8, weight='normal')
    axes[1].set_title('Logistic regression feature weights', fontsize=10)
    axes[1].set_xlabel(None)
    axes[1].set_ylabel('Weights')
    axes[1].set_xticks(range(0, matrix_circ_norm_final.shape[1]))
    axes[1].set_xticklabels(top_feature_names, rotation=90)
    axes[1].legend([Patch(facecolor='lavender', edgecolor='k'),
                    Patch(facecolor='lavender', edgecolor='k', hatch='\\')], ['Positive', 'Negative'])
    plt.tight_layout()
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/Circadian_in_Cisplatin_out/logreg/circadian_Cisplatin_'
                    f'BMAL1_PER2_{measure}_logreg.{file}', bbox_inches='tight', dpi=300)
    plt.close()
    # logreg LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogReg_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug)
    loocv_results_grid_circadian[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
matplotlib.rcParams['hatch.linewidth'] = original_hatch_width

# Plot for CCLE input
loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
columns = set()
for col in loocv_samples_ccle.columns:
    for index_object in loocv_samples_ccle[col]:
        if isinstance(index_object, pd.Index):
            columns.update(index_object)
columns = list(columns)
indeces = loocv_samples_ccle.index
this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
for drug in indeces:
    for sample in columns:
        # case of correct samples and wrong samples
        if sample in loocv_samples_ccle.loc[drug]['corr_samples'].__str__():
            this_table_ccle.loc[drug][sample] += 1
        # case of the N/A
        if sample not in loocv_samples_ccle.loc[drug]['corr_samples'].__str__() and sample not in \
                loocv_samples_ccle.loc[drug]['not_corr_samples'].__str__():
            this_table_ccle.loc[drug][sample] = np.nan
# add accuracy
this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
sorted_table = sorted_table.dropna(axis=0, how='all')
# table figure
# Create a mask for the 10th column and nans
nanmask = sorted_table.isna()
mask = np.zeros_like(sorted_table, dtype=bool)
mask[:, -1] = ~nanmask.iloc[:, -1]
combinedmask = mask | nanmask
fig, axe = plt.subplots(1, 1, figsize=(6, 6))
sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
            annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
# set the np.nan tiles to grey
axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
           extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
for i in range(sorted_table.shape[0]):
    axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
             ha='center', va='center', color='black', fontsize=8)
samples_per_drug = []
for tod_drug in sorted_table.index:
    samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
# set the labels plus sample size
axe.set_yticklabels(y_labels, rotation=0, fontsize=12)
# add frame lines around heatmap
axe.axhline(y=0, color='k', linewidth=2)
axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
axe.axvline(x=0, color='k', linewidth=2)
axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
axe.xaxis.tick_top()  # x axis on top
axe.xaxis.set_label_position('top')
axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
axe.set_title(f"LogReg leave-one-out classification results by cell line, Cisplatin (CCLE)\nsample size in parenthesis")
for file in formats:
    plt.savefig(f'./results_revised_review/comment_5/CCLE_in_Cisplatin_out/logreg/CCLE_Cisplatin_loocv_'
                f'classification_plot_LogReg.{file}', bbox_inches='tight', dpi=300)
plt.close()

# Plot for Circadian input
loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
columns = set()
for col in loocv_samples_circadian.columns:
    for index_object in loocv_samples_circadian[col]:
        if isinstance(index_object, pd.Index):
            columns.update(index_object)
columns = list(columns)
indeces = loocv_samples_circadian.index
this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
for drug in indeces:
    for sample in columns:
        # case of correct samples and wrong samples
        if sample in loocv_samples_circadian.loc[drug]['corr_samples'].__str__():
            this_table_circadian.loc[drug][sample] += 1
        # case of the N/A
        if sample not in loocv_samples_circadian.loc[drug]['corr_samples'].__str__() and sample not in \
                loocv_samples_circadian.loc[drug]['not_corr_samples'].__str__():
            this_table_circadian.loc[drug][sample] = np.nan
# add accuracy
this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
sorted_table = sorted_table.dropna(axis=0, how='all')
# table figure
# Create a mask for the 10th column and nans
nanmask = sorted_table.isna()
mask = np.zeros_like(sorted_table, dtype=bool)
mask[:, -1] = ~nanmask.iloc[:, -1]
combinedmask = mask | nanmask
fig, axe = plt.subplots(1, 1, figsize=(6, 6))
sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
            annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
# set the np.nan tiles to grey
axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
           extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
for i in range(sorted_table.shape[0]):
    axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
             ha='center', va='center', color='black', fontsize=8)
samples_per_drug = []
for tod_drug in sorted_table.index:
    samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
# set the labels plus sample size
axe.set_yticklabels(y_labels, rotation=0, fontsize=12)
# add frame lines around heatmap
axe.axhline(y=0, color='k', linewidth=2)
axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
axe.axvline(x=0, color='k', linewidth=2)
axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
axe.xaxis.tick_top()  # x axis on top
axe.xaxis.set_label_position('top')
axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
axe.set_title(f"LogReg leave-one-out classification results by cell line, Cisplatin (BMAL1-PER2 mean)\n"
              f"sample size in parenthesis")
for file in formats:
    plt.savefig(f'./results_revised_review/comment_5/Circadian_in_Cisplatin_out/logreg/Circadian_Cisplatin_'
                f'BMAL1_PER2_loocv_classification_plot_LogReg.{file}', bbox_inches='tight', dpi=300)
plt.close()

# - comment 5: Logistic regression on binarized data, extract staitstical values, and perform LOOCV (CCLE sensitivity)

drugs = remaining_drugs.copy()
tod_scores_for_loocv = ['EC50 (uM)', 'IC50 (uM)', 'ActArea']

# this time we are NOT going to limit EC50 and ActArea drugs to the ones where IC50 is at least available

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)
loocv_results_grid_circadian = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)

original_hatch_width = 1.0
matplotlib.rcParams['hatch.linewidth'] = 8.0

# For each drug
for measure in tod_scores_for_loocv:
    print(f'* Analysing stat {measure}...')
    df = remaining_ccle_drugs.copy()
    df = df[['Primary Cell Line Name', 'Compound', measure]]
    for drug in remaining_drugs:
        print(f'*** Analysing for drug {drug}...')
        # get this drug only
        this_drug = df.loc[df['Compound'] == drug]
        this_drug_end = this_drug.dropna()
        # the number of samples must be more than the number of classes, so we always need at least 3 samples
        if len(this_drug_end) <= 3:
            # somtimes, IC50 has mostly the value 8, resulting in all classes being True / False
            text = f'only {len(this_drug_end)} cell line was' if len(this_drug_end) == 1 else \
                f'only {len(this_drug_end)} cell line were' if len(this_drug_end) > 1 else 'no cell lines were'
            print(f'For the drug {drug}, {text} measured for {measure}. This example is skipped for the LDA analysis.')
            continue  # skip if we only remain with 1 or 0 samples (most cases are EC50 (uM) not being measured)
        print(f'** Discretize drug sensitivity by median...')
        this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median()
        this_drug_binarized.name = drug
        if len(this_drug_binarized.unique()) == 1:  # apply the -0.001 trick only if necessary!
            this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median() - 0.001
            print('## Trying "median - 0.001" for discretiztation cutoff ##')
        if len(this_drug_binarized.unique()) == 1:
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in one '
                  f'single class. This example is skipped for the LDA analysis.')
            continue
        if 1 in set(this_drug_binarized.value_counts()):  # skip LDA loocv of samples where ther is only 1 true or false
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in two classes but '
                  f'one is only represented once. This example is skipped for the LDA LOOCV analysis.')
            continue
        subtypes_here = pd.DataFrame(
            [dict_to_map_line_and_subtype[x] for x in this_drug_end['Primary Cell Line Name']],
            index=this_drug_end['Primary Cell Line Name'], columns=['subtype'])
        this_drug_binarized.index = subtypes_here.index
        # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
        print('**** Logistic regression on CCLE expression...')
        ccle_expression = ccle_data.loc[this_drug_end['Primary Cell Line Name']]
        # Logreg on full data
        logreg = LogisticRegression(random_state=seed).fit(ccle_expression, this_drug_binarized)
        probs = logreg.predict_proba(ccle_expression)
        # Extract probabilities for the positive class
        positive_probs = probs[:, 1]
        # Extract feature importance
        importance = logreg.coef_[0]
        sorted_indices = np.argsort(np.abs(importance))[::-1]  # Sort by absolute importance
        top_importance = importance[sorted_indices]
        top_feature_names = np.array(ccle_expression.columns)[sorted_indices]
        # Plot probability distributions and feature importance
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        sns.histplot(positive_probs[this_drug_binarized == 1], bins='fd',
                     color="green", label="Positive Class", kde=True, ax=axes[0])  # Freedman Diaconis for n_bins
        sns.histplot(positive_probs[this_drug_binarized == 0], bins='fd',
                     color="red", label="Negative Class", kde=True, ax=axes[0])
        axes[0].set_title(f"Logsitic regression training performance, {drug} {measure} (CCLE)")
        axes[0].set_xlabel("Predicted Probability")
        axes[0].set_ylabel("Frequency")
        axes[0].legend()
        # Feature importance
        axes[1].bar(range(0, ccle_expression.shape[1]), np.abs(top_importance), color='lavender', edgecolor='k',
                    linewidth=1, hatch=['\\' if i < 0 else None for i in top_importance])
        for i in range(0, ccle_expression.shape[1]):
            axes[1].annotate(str("{:.2f}".format(round(top_importance[i], 2))),
                             xy=(i, np.abs(top_importance[i])), ha='center', va='bottom', size=8, weight='normal')
        axes[1].set_title('Logistic regression feature weights', fontsize=10)
        axes[1].set_xlabel(None)
        axes[1].set_ylabel('Weights')
        axes[1].set_xticks(range(0, ccle_expression.shape[1]))
        axes[1].set_xticklabels(top_feature_names, rotation=90)
        axes[1].legend([Patch(facecolor='lavender', edgecolor='k'),
                        Patch(facecolor='lavender', edgecolor='k', hatch='\\')], ['Positive', 'Negative'])
        plt.tight_layout()
        for file in formats:
            plt.savefig(f'./results_revised_review/comment_5/CCLE_in_CCLE_out/logreg/CCLE_{drug}_{measure}_'
                        f'logreg.{file}', bbox_inches='tight', dpi=300)
        plt.close()
        # LogReg LOOCV machinery
        this_drug_binarized.name = drug  # make sure the targets have the right name
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogReg_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                         target=drug)
        loocv_results_grid_ccle[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
        # •	By circadian features mean
        print(f'**** Logistic regression on BMAL1_PER2 mean parameters...\n')
        mat = circadian_data.copy()
        matrix_circ = mat.dropna(axis=1, how='all')
        matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
        matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
        matrix_circ_norm = matrix_circ_norm.loc[subtypes_here.index]
        # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
        # sample from the drug to keep consistent data sets
        all_together = pd.concat([matrix_circ_norm, this_drug_binarized, subtypes_here], axis=1)
        all_together_final = all_together.dropna(axis=0)
        matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
        this_drug_binarized_consistent = all_together_final[drug]
        if 1 in this_drug_binarized_consistent.value_counts().values or len(
                all_together_final) <= 3:  # also here <=3 instead of 2
            print(f'/!\\ Logistic regression on circadian parameters targeting {drug} is only possible '
                  f'for\n{len(all_together_final)} samples where parameters and {measure} are available,'
                  f'\nbut the binarization of these remaining samples is left with 1 class only using median.'
                  f'\nLDA requires at least 2 classes and 3 samples. /!\\')
            continue
        # Logreg on full data
        logreg = LogisticRegression(random_state=seed).fit(matrix_circ_norm_final, this_drug_binarized_consistent)
        probs = logreg.predict_proba(matrix_circ_norm_final)
        # Extract probabilities for the positive class
        positive_probs = probs[:, 1]
        # Extract feature importance
        importance = logreg.coef_[0]
        sorted_indices = np.argsort(np.abs(importance))[::-1]  # Sort by absolute importance
        top_importance = importance[sorted_indices]
        top_feature_names = np.array(matrix_circ_norm_final.columns)[sorted_indices]
        # Plot probability distributions and feature importance
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        sns.histplot(positive_probs[this_drug_binarized_consistent == 1], color="green", label="Positive Class",
                     bins='fd', kde=True, ax=axes[0])
        sns.histplot(positive_probs[this_drug_binarized_consistent == 0], color="red", label="Negative Class",
                     bins='fd', kde=True, ax=axes[0])
        axes[0].set_title(f"Logsitic regression training performance, {drug} {measure} (BMAL1-PER2 mean)")
        axes[0].set_xlabel("Predicted Probability")
        axes[0].set_ylabel("Frequency")
        axes[0].legend()
        # Feature importance
        axes[1].bar(range(0, matrix_circ_norm_final.shape[1]), np.abs(top_importance), color='lavender', edgecolor='k',
                    linewidth=1, hatch=['\\' if i < 0 else None for i in top_importance])
        for i in range(0, matrix_circ_norm_final.shape[1]):
            axes[1].annotate(str("{:.2f}".format(round(top_importance[i], 2))),
                             xy=(i, np.abs(top_importance[i])), ha='center', va='bottom', size=8, weight='normal')
        axes[1].set_title('Logistic regression feature weights', fontsize=10)
        axes[1].set_xlabel(None)
        axes[1].set_ylabel('Weights')
        axes[1].set_xticks(range(0, matrix_circ_norm_final.shape[1]))
        axes[1].set_xticklabels(top_feature_names, rotation=90)
        axes[1].legend([Patch(facecolor='lavender', edgecolor='k'),
                        Patch(facecolor='lavender', edgecolor='k', hatch='\\')], ['Positive', 'Negative'])
        plt.tight_layout()
        for file in formats:
            plt.savefig(f'./results_revised_review/comment_5/Circadian_in_CCLE_out/logreg/circadian_{drug}_'
                        f'BMAL1_PER2_{measure}_logreg.{file}', bbox_inches='tight', dpi=300)
        plt.close()
        # logreg LOOCV machinery
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogReg_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug)
        loocv_results_grid_circadian[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # plot at this level for each individual measure (target), 2 per measure (CCLE and circ (including channels))
    # Plot for CCLE sensitivity input
    loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_ccle.columns:
        for index_object in loocv_samples_ccle[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_ccle.index
    this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for d in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_ccle.loc[d]['corr_samples'].__str__():
                this_table_ccle.loc[d][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_ccle.loc[d]['corr_samples'].__str__() and sample not in \
                    loocv_samples_ccle.loc[d][
                        'not_corr_samples'].__str__():
                this_table_ccle.loc[d][sample] = np.nan
    # add accuracy
    this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
    sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig, axe = plt.subplots(1, 1, figsize=(6, 6))
    sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
    # set the np.nan tiles to grey
    axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    axe.set_yticklabels(y_labels, fontsize=12)
    # add frame lines around heatmap
    axe.axhline(y=0, color='k', linewidth=2)
    axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    axe.axvline(x=0, color='k', linewidth=2)
    axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    axe.xaxis.tick_top()  # x axis on top
    axe.xaxis.set_label_position('top')
    axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    axe.set_title(f"LogReg leave-one-out classification results by cell line, {measure} (CCLE)\n"
                  f"sample size in parenthesis")
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/CCLE_in_CCLE_out/logreg/CCLE_{measure}_loocv_'
                    f'classification_plot_LogReg.{file}', bbox_inches='tight', dpi=300)
    plt.close()
    # Plot for Circadian input
    loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_circadian.columns:
        for index_object in loocv_samples_circadian[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_circadian.index
    this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for d in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_circadian.loc[d]['corr_samples'].__str__():
                this_table_circadian.loc[d][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_circadian.loc[d]['corr_samples'].__str__() and sample not in \
                    loocv_samples_circadian.loc[d][
                        'not_corr_samples'].__str__():
                this_table_circadian.loc[d][sample] = np.nan
    # add accuracy
    this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
    sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig, axe = plt.subplots(1, 1, figsize=(6, 6))
    sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
    # set the np.nan tiles to grey
    axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    axe.set_yticklabels(y_labels, fontsize=12)
    # add frame lines around heatmap
    axe.axhline(y=0, color='k', linewidth=2)
    axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    axe.axvline(x=0, color='k', linewidth=2)
    axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    axe.xaxis.tick_top()  # x axis on top
    axe.xaxis.set_label_position('top')
    axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    axe.set_title(f"LogReg leave-one-out classification results by cell line on {measure} (BMAL1-PER2 mean)\n"
                  f"sample size in parenthesis")
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/Circadian_in_CCLE_out/logreg/circadian_{measure}_BMAL1_PER2_'
                    f'loocv_classification_plot_LogReg.{file}', bbox_inches='tight', dpi=300)
    plt.close()
matplotlib.rcParams['hatch.linewidth'] = original_hatch_width


# - comment 5: Ridge/Lasso on continuous senstivity on CISPLATIN only (n=15), get a goodness-of-fit param and do LOOCV

# ! In case of Ridge error: sym_pos=True unknown argument, go to the sklearn '_ridge.py' script and exchange
# 'sym_pos=True' with 'assume_a='pos' !

seed = 42

drug = 'Cisplatin'
tod_scores_for_loocv = list(data_sensitivty.keys())

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['mse', 'mae', 'r2', 'preds', 'abs_diffs', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle_ridge = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                             index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_circadian_ridge = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                                  index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_ccle_lasso = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                             index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_circadian_lasso = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                                  index=metrics, columns=tod_scores_for_loocv)

original_hatch_width = 1.0
matplotlib.rcParams['hatch.linewidth'] = 8.0

correctness_cutoff = 'below ten'  # 'mean diff', 'below five', or 'below ten'

# Without sample MCF10A
for measure in list(data_sensitivty.keys()):
    print(f'* Analysing stat {measure}...')
    df = data_sensitivty[measure].T.copy()
    df = df.drop('MCF10A', axis=0)
    this_drug = df[drug]
    this_drug_end = this_drug.dropna()
    # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
    print('**** Lasso and Ridge regression on CCLE expression...')
    ccle_expression = ccle_data.loc[this_drug_end.index]
    lasso = Lasso(alpha=0.1, random_state=seed).fit(ccle_expression, this_drug_end)  # L1 regularization
    ridge = Ridge(alpha=1.0, random_state=seed).fit(ccle_expression, this_drug_end)  # L2 regularization
    # Evaluate lasso and ridge
    plot_linear_regressor(estimator=lasso, x=ccle_expression, y=this_drug_end,
                          title=f'Lasso regression performance, CCLE on Cisplatin {measure}')
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/CCLE_in_Cisplatin_out/lasso/CCLE_Cisplatin_{measure}_lasso'
                    f'.{file}', bbox_inches='tight', dpi=300)
    plt.close()
    plot_linear_regressor(estimator=ridge, x=ccle_expression, y=this_drug_end,
                          title=f'Ridge regression performance, CCLE on Cisplatin {measure}')
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/CCLE_in_Cisplatin_out/ridge/CCLE_Cisplatin_{measure}_ridge'
                    f'.{file}', bbox_inches='tight', dpi=300)
    plt.close()
    mse, mae, r2, preds, abs_diffs, corr_samples, not_corr_samples = \
        regressor_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_end], axis=1), target=drug,
                        regressor='lasso', cut_off=correctness_cutoff)
    loocv_results_grid_ccle_lasso[measure] = [mse, mae, r2, preds, abs_diffs, corr_samples, not_corr_samples]
    mse, mae, r2, preds, abs_diffs, corr_samples, not_corr_samples = \
        regressor_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_end], axis=1), target=drug,
                        regressor='ridge', cut_off=correctness_cutoff)
    loocv_results_grid_ccle_ridge[measure] = [mse, mae, r2, preds, abs_diffs, corr_samples, not_corr_samples]
    # •	By circadian features mean
    print(f'**** Lasso and Ridge on BMAL1_PER2 mean parameters...\n')
    mat = circadian_data.copy()
    matrix_circ = mat.dropna(axis=1, how='all')
    matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
    matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
    matrix_circ_norm = matrix_circ_norm.loc[this_drug_end.index]
    # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
    # sample from the drug to keep consistent data sets
    all_together = pd.concat([matrix_circ_norm, this_drug_end], axis=1)
    all_together_final = all_together.dropna(axis=0)
    matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
    this_drug_end_consistent = all_together_final[drug]
    lasso = Lasso(alpha=0.1, random_state=seed).fit(matrix_circ_norm_final, this_drug_end_consistent)  # L1
    ridge = Ridge(alpha=1.0, random_state=seed).fit(matrix_circ_norm_final, this_drug_end_consistent)  # L2
    # Evaluate lasso and ridge
    plot_linear_regressor(estimator=lasso, x=matrix_circ_norm_final, y=this_drug_end_consistent,
                          title=f'Lasso regression performance, circadian on Cisplatin {measure}')
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/Circadian_in_Cisplatin_out/lasso/circadian_Cisplatin_{measure}'
                    f'_lasso.{file}', bbox_inches='tight', dpi=300)
    plt.close()
    plot_linear_regressor(estimator=ridge, x=matrix_circ_norm_final, y=this_drug_end_consistent,
                          title=f'Ridge regression performance, circadian on Cisplatin {measure}')
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/Circadian_in_Cisplatin_out/ridge/circadian_Cisplatin_{measure}'
                    f'_ridge.{file}', bbox_inches='tight', dpi=300)
    plt.close()
    mse, mae, r2, preds, abs_diffs, corr_samples, not_corr_samples = \
        regressor_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug, regressor='lasso',
                        cut_off=correctness_cutoff)
    loocv_results_grid_circadian_lasso[measure] = [mse, mae, r2, preds, abs_diffs, corr_samples, not_corr_samples]
    mse, mae, r2, preds, abs_diffs, corr_samples, not_corr_samples = \
        regressor_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug, regressor='ridge',
                        cut_off=correctness_cutoff)
    loocv_results_grid_circadian_ridge[measure] = [mse, mae, r2, preds, abs_diffs, corr_samples, not_corr_samples]
matplotlib.rcParams['hatch.linewidth'] = original_hatch_width

# Get tables of the predictions and true values
sorted_tables_predicted = {name: pd.DataFrame(columns=['GRAOC', 'GRINF', 'GEC50', 'GR50', 'Hill']) for name in
                           ['CCLE lasso', 'CCLE ridge', 'Circadian lasso', 'Circadian ridge']}
sorted_tables_actual = {name: pd.DataFrame(columns=['GRAOC', 'GRINF', 'GEC50', 'GR50', 'Hill']) for name in
                        ['CCLE lasso', 'CCLE ridge', 'Circadian lasso', 'Circadian ridge']}

# Finally, we will plot again the LOOCV results as heatmaps
loocv_grids = [loocv_results_grid_ccle_lasso, loocv_results_grid_ccle_ridge,
               loocv_results_grid_circadian_lasso, loocv_results_grid_circadian_ridge]

fig, axes = plt.subplots(2, 2, figsize=(14, 6))
for grid, ax, name in zip(loocv_grids, axes.ravel(), ['CCLE lasso', 'CCLE ridge',
                                                      'Circadian lasso', 'Circadian ridge']):
    preds_here = grid.T[['corr_samples', 'not_corr_samples']]
    abs_diffs_here = grid.T['abs_diffs']
    pred_values_here = grid.T['preds']
    r2_here = grid.T['r2']
    columns = set()
    for col in preds_here.columns:
        for index_object in preds_here[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = preds_here.index
    this_table = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for drug in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in preds_here.loc[drug]['corr_samples'].__str__():
                this_table.loc[drug][sample] += 1
    # add accuracy and r2
    acc_here = np.mean(this_table, axis=1).astype(float)
    this_table['r2'] = r2_here.astype(float)
    this_table['accuracy'] = acc_here
    sorted_table = this_table.sort_values(by='r2', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    sorted_tables_predicted[name].index = sorted_table.columns[:-2]
    sorted_tables_actual[name].index = sorted_table.columns[:-2]
    # table figure
    # Create a mask for the 2 last columns and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -2:] = ~nanmask.iloc[:, -2:]
    combinedmask = mask | nanmask
    # After we got the mask where heatmap should not be colorized, we could change the colors to adapt for continuous
    for meas in sorted_table.index:
        trues = sorted_table.loc[meas][:-2].copy()
        for cell in sorted_table.columns[:-2]:
            sorted_table.loc[meas, cell] = abs_diffs_here.loc[meas][cell]
        mini, maxi = sorted_table.loc[meas][:-2].min(), sorted_table.loc[meas][:-2].max()
        if correctness_cutoff == 'mean diff':
            center = sorted_table.loc[meas][:-2].mean()
        elif correctness_cutoff in ('below five', 'below ten'):
            center = sorted_table.loc[meas][
                trues.index[trues == 0]].min() if (1 in trues.values and 0 in trues.values) else \
                sorted_table.loc[meas][trues.index[trues == 1]].max() if 1 in trues.values else \
                sorted_table.loc[meas][trues.index[trues == 0]].min()
        else:
            center = 0.1
        mask_here = combinedmask.copy()
        mask_here[mask_here.index != meas] = True
        sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn_r'), mask=mask_here, cbar=False,
                    annot=False, norm=TwoSlopeNorm(vmin=mini if sum(trues.values) >= 1 else 0.999 * mini,
                                                   vcenter=center,
                                                   vmax=maxi if sum(trues.values) != 0 else 1.001 * maxi),
                    linewidth=.003, ax=ax)
    for i in range(sorted_table.shape[0]):
        ax.text(sorted_table.shape[1] - 1.5, i + 0.5, f'{sorted_table.iloc[i, -2]:.2f}',
                ha='center', va='center', color='black', fontsize=8)
        ax.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                ha='center', va='center', color='black', fontsize=8)
        for t in range(sorted_table.shape[1] - 2):
            true_here = data_sensitivty[sorted_table.index[i]].loc["Cisplatin", sorted_table.columns[t]]
            predicted_here = pred_values_here.loc[sorted_table.index[i]][
                abs_diffs_here.loc[sorted_table.index[i]].index == sorted_table.columns[t]][0]
            abs_difference_here = abs_diffs_here.loc[sorted_table.index[i]][sorted_table.columns[t]]
            # update the table for later loocv regression plots
            sorted_tables_predicted[name].loc[sorted_table.columns[t], sorted_table.index[i]] = predicted_here
            sorted_tables_actual[name].loc[sorted_table.columns[t], sorted_table.index[i]] = true_here
            ax.text(t + 0.5, i + 0.5, f't: {true_here:.2f}\np: {predicted_here:.2f}\nd: {abs_difference_here:.2f}',
                    ha='center', va='center', color='black', fontsize=6)
    # add frame lines around heatmap
    ax.axhline(y=0, color='k', linewidth=2)
    ax.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    ax.axvline(x=0, color='k', linewidth=2)
    ax.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    ax.axvline(x=sorted_table.iloc[:, :-2].shape[1], color='k', linewidth=1)
    ax.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    ax.xaxis.tick_top()  # x axis on top
    ax.xaxis.set_label_position('top')
    # set the labels plus sample size
    ax.set_yticklabels(sorted_table.index, rotation=0, fontsize=10)
    if name in ('CCLE lasso', 'CCLE ridge'):
        ax.set_xticklabels(sorted_table.columns, fontsize=10, rotation=75)
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])
    ax.set_title('$\\bf{Lasso\\ regression}$' if name == 'CCLE lasso' else
                 '$\\bf{Ridge\\ regression}$' if name == 'CCLE ridge' else '', fontsize=12)
    ax.set_ylabel('$\\bf{CCLE\\ input}$' if name == 'CCLE lasso' else
                  '$\\bf{Circadian\\ input}$' if name == 'Circadian lasso' else '',
                  fontsize=12)
fig.suptitle(f"Lasso and Ridge leave-one-out Cisplatin sensitivity prediction by cell line\n"
             f"(n={len(columns)}, threshold: {correctness_cutoff})")
fig.tight_layout()
for file in formats:
    plt.savefig(f'./results_revised_review/comment_5/Lass_Ridge_cisplatin_loocv_classification_plot_'
                f'{correctness_cutoff.replace(" ", "-")}.{file}', bbox_inches='tight', dpi=300)
plt.close()

# Plot the loocv predictions versus actual (the correctness cut off does not matter here, sorted_tables are identical)
for name in ['CCLE lasso', 'CCLE ridge', 'Circadian lasso', 'Circadian ridge']:
    for meas in sorted_tables_actual[name].columns:
        preds = sorted_tables_predicted[name][meas].astype(float)
        trues = sorted_tables_actual[name][meas].astype(float)
        abs_diff = np.abs(preds - trues)
        plot_regressor_loocv(preds, trues, title=f'LOOCV {"Lasso" if "lasso" in name else "Ridge"} predicted '
                                                 f'versus actual Cisplatin sensitivity ({name.split(" ")[0]})')
        for file in formats:
            plt.savefig(f'./results_revised_review/comment_5/LOOCV_plot_Cisplatin_{name.replace(" ", "_")}_{meas}.'
                        f'{file}', bbox_inches='tight', dpi=300)
        plt.close()

####################
# ## Alternatives ##

# - comment 5: Alternative: LogReg on LDA transformed data as LOOCV
# LDA loocv with only transformation, followed by LogReg training on the all-except-one, then predict the left out

drug = 'Cisplatin'
tod_scores_for_loocv = list(data_sensitivty.keys())

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                       index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_circadian = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                            index=metrics, columns=tod_scores_for_loocv)

# Without sample MCF10A
for measure in list(data_sensitivty.keys()):
    print(f'* Analysing stat {measure}...')
    df = data_sensitivty[measure].T.copy()
    df = df.drop('MCF10A', axis=0)
    this_drug = df[drug]
    this_drug_end = this_drug.dropna()
    print(f'** Discretize drug sensitivity by median...')
    this_drug_binarized = this_drug_end > this_drug_end.median()
    # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
    print('**** Logistic regression on CCLE expression...')
    ccle_expression = ccle_data.loc[this_drug_end.index]
    # LogReg on LDA transformed LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogReg_on_LDA_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                            target=drug)
    loocv_results_grid_ccle[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # •	By circadian features mean
    print(f'**** Logistic regression on BMAL1_PER2 mean parameters...\n')
    mat = circadian_data.copy()
    matrix_circ = mat.dropna(axis=1, how='all')
    matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
    matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
    matrix_circ_norm = matrix_circ_norm.loc[this_drug_end.index]
    # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
    # sample from the drug to keep consistent data sets
    all_together = pd.concat([matrix_circ_norm, this_drug_binarized], axis=1)
    all_together_final = all_together.dropna(axis=0)
    matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
    this_drug_binarized_consistent = all_together_final[drug]
    # logreg on LDA transformed LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogReg_on_LDA_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug)
    loocv_results_grid_circadian[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]

# Plot for CCLE input
loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
columns = set()
for col in loocv_samples_ccle.columns:
    for index_object in loocv_samples_ccle[col]:
        if isinstance(index_object, pd.Index):
            columns.update(index_object)
columns = list(columns)
indeces = loocv_samples_ccle.index
this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
for drug in indeces:
    for sample in columns:
        # case of correct samples and wrong samples
        if sample in loocv_samples_ccle.loc[drug]['corr_samples'].__str__():
            this_table_ccle.loc[drug][sample] += 1
        # case of the N/A
        if sample not in loocv_samples_ccle.loc[drug]['corr_samples'].__str__() and sample not in \
                loocv_samples_ccle.loc[drug]['not_corr_samples'].__str__():
            this_table_ccle.loc[drug][sample] = np.nan
# add accuracy
this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
sorted_table = sorted_table.dropna(axis=0, how='all')
# table figure
# Create a mask for the 10th column and nans
nanmask = sorted_table.isna()
mask = np.zeros_like(sorted_table, dtype=bool)
mask[:, -1] = ~nanmask.iloc[:, -1]
combinedmask = mask | nanmask
fig, axe = plt.subplots(1, 1, figsize=(6, 6))
sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
            annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
# set the np.nan tiles to grey
axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
           extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
for i in range(sorted_table.shape[0]):
    axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
             ha='center', va='center', color='black', fontsize=8)
samples_per_drug = []
for tod_drug in sorted_table.index:
    samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
# set the labels plus sample size
axe.set_yticklabels(y_labels, rotation=0, fontsize=12)
# add frame lines around heatmap
axe.axhline(y=0, color='k', linewidth=2)
axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
axe.axvline(x=0, color='k', linewidth=2)
axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
axe.xaxis.tick_top()  # x axis on top
axe.xaxis.set_label_position('top')
axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
axe.set_title(f"LogReg on LDA leave-one-out classification results by cell line, Cisplatin (CCLE)\n"
              f"sample size in parenthesis")
for file in formats:
    plt.savefig(f'./results_revised_review/comment_5/CCLE_in_Cisplatin_out/logreg_on_lda/CCLE_Cisplatin_loocv_'
                f'classification_plot_LogReg.{file}', bbox_inches='tight', dpi=300)
plt.close()

# Plot for Circadian input
loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
columns = set()
for col in loocv_samples_circadian.columns:
    for index_object in loocv_samples_circadian[col]:
        if isinstance(index_object, pd.Index):
            columns.update(index_object)
columns = list(columns)
indeces = loocv_samples_circadian.index
this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
for drug in indeces:
    for sample in columns:
        # case of correct samples and wrong samples
        if sample in loocv_samples_circadian.loc[drug]['corr_samples'].__str__():
            this_table_circadian.loc[drug][sample] += 1
        # case of the N/A
        if sample not in loocv_samples_circadian.loc[drug]['corr_samples'].__str__() and sample not in \
                loocv_samples_circadian.loc[drug]['not_corr_samples'].__str__():
            this_table_circadian.loc[drug][sample] = np.nan
# add accuracy
this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
sorted_table = sorted_table.dropna(axis=0, how='all')
# table figure
# Create a mask for the 10th column and nans
nanmask = sorted_table.isna()
mask = np.zeros_like(sorted_table, dtype=bool)
mask[:, -1] = ~nanmask.iloc[:, -1]
combinedmask = mask | nanmask
fig, axe = plt.subplots(1, 1, figsize=(6, 6))
sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
            annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
# set the np.nan tiles to grey
axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
           extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
for i in range(sorted_table.shape[0]):
    axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
             ha='center', va='center', color='black', fontsize=8)
samples_per_drug = []
for tod_drug in sorted_table.index:
    samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
# set the labels plus sample size
axe.set_yticklabels(y_labels, rotation=0, fontsize=12)
# add frame lines around heatmap
axe.axhline(y=0, color='k', linewidth=2)
axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
axe.axvline(x=0, color='k', linewidth=2)
axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
axe.xaxis.tick_top()  # x axis on top
axe.xaxis.set_label_position('top')
axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
axe.set_title(f"LogReg on LDA leave-one-out classification results by cell line, Cisplatin (BMAL1-PER2 mean)\n"
              f"sample size in parenthesis")
for file in formats:
    plt.savefig(f'./results_revised_review/comment_5/Circadian_in_Cisplatin_out/logreg_on_lda/Circadian_Cisplatin_'
                f'BMAL1_PER2_loocv_classification_plot_LogReg.{file}', bbox_inches='tight', dpi=300)
plt.close()

# Logreg on LDA transformed data using CCLE drug sensitivities

drugs = remaining_drugs.copy()
tod_scores_for_loocv = ['EC50 (uM)', 'IC50 (uM)', 'ActArea']

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)
loocv_results_grid_circadian = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)

# For each drug
for measure in tod_scores_for_loocv:
    print(f'* Analysing stat {measure}...')
    df = remaining_ccle_drugs.copy()
    df = df[['Primary Cell Line Name', 'Compound', measure]]
    for drug in remaining_drugs:
        print(f'*** Analysing for drug {drug}...')
        # get this drug only
        this_drug = df.loc[df['Compound'] == drug]
        this_drug_end = this_drug.dropna()
        # the number of samples must be more than the number of classes, so we always need at least 3 samples
        if len(this_drug_end) <= 3:
            # somtimes, IC50 has mostly the value 8, resulting in all classes being True / False
            text = f'only {len(this_drug_end)} cell line was' if len(this_drug_end) == 1 else \
                f'only {len(this_drug_end)} cell line were' if len(this_drug_end) > 1 else 'no cell lines were'
            print(f'For the drug {drug}, {text} measured for {measure}. This example is skipped for the LDA analysis.')
            continue  # skip if we only remain with 1 or 0 samples (most cases are EC50 (uM) not being measured)
        print(f'** Discretize drug sensitivity by median...')
        this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median()
        this_drug_binarized.name = drug
        if len(this_drug_binarized.unique()) == 1:  # apply the -0.001 trick only if necessary!
            this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median() - 0.001
            print('## Trying "median - 0.001" for discretiztation cutoff ##')
        if len(this_drug_binarized.unique()) == 1:
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in one '
                  f'single class. This example is skipped for the LDA analysis.')
            continue
        if 1 in set(this_drug_binarized.value_counts()):  # skip LDA loocv of samples where ther is only 1 true or false
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in two classes but '
                  f'one is only represented once. This example is skipped for the LDA LOOCV analysis.')
            continue
        subtypes_here = pd.DataFrame(
            [dict_to_map_line_and_subtype[x] for x in this_drug_end['Primary Cell Line Name']],
            index=this_drug_end['Primary Cell Line Name'], columns=['subtype'])
        this_drug_binarized.index = subtypes_here.index
        # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
        print('**** Logistic regression on CCLE expression...')
        ccle_expression = ccle_data.loc[this_drug_end['Primary Cell Line Name']]
        # LogReg on LDA transformed LOOCV machinery
        this_drug_binarized.name = drug  # make sure the targets have the right name
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogReg_on_LDA_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                                target=drug)
        loocv_results_grid_ccle[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
        # •	By circadian features mean
        print(f'**** Logistic regression on BMAL1_PER2 mean parameters...\n')
        mat = circadian_data.copy()
        matrix_circ = mat.dropna(axis=1, how='all')
        matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
        matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
        matrix_circ_norm = matrix_circ_norm.loc[subtypes_here.index]
        # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
        # sample from the drug to keep consistent data sets
        all_together = pd.concat([matrix_circ_norm, this_drug_binarized, subtypes_here], axis=1)
        all_together_final = all_together.dropna(axis=0)
        matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
        this_drug_binarized_consistent = all_together_final[drug]
        if 1 in this_drug_binarized_consistent.value_counts().values or len(
                all_together_final) <= 3:  # also here <=3 instead of 2
            print(f'/!\\ Logistic regression on circadian parameters targeting {drug} is only possible '
                  f'for\n{len(all_together_final)} samples where parameters and {measure} are available,'
                  f'\nbut the binarization of these remaining samples is left with 1 class only using median.'
                  f'\nLDA requires at least 2 classes and 3 samples. /!\\')
            continue
        # logreg on LDA transformed LOOCV machinery
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogReg_on_LDA_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug)
        loocv_results_grid_circadian[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # plot at this level for each individual measure (target), 2 per measure (CCLE and circ (including channels))
    # Plot for CCLE sensitivity input
    loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_ccle.columns:
        for index_object in loocv_samples_ccle[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_ccle.index
    this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for d in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_ccle.loc[d]['corr_samples'].__str__():
                this_table_ccle.loc[d][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_ccle.loc[d]['corr_samples'].__str__() and sample not in \
                    loocv_samples_ccle.loc[d][
                        'not_corr_samples'].__str__():
                this_table_ccle.loc[d][sample] = np.nan
    # add accuracy
    this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
    sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig, axe = plt.subplots(1, 1, figsize=(6, 6))
    sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
    # set the np.nan tiles to grey
    axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    axe.set_yticklabels(y_labels, fontsize=12)
    # add frame lines around heatmap
    axe.axhline(y=0, color='k', linewidth=2)
    axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    axe.axvline(x=0, color='k', linewidth=2)
    axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    axe.xaxis.tick_top()  # x axis on top
    axe.xaxis.set_label_position('top')
    axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    axe.set_title(f"LogReg on LDA leave-one-out classification results by cell line, {measure} (CCLE)\n"
                  f"sample size in parenthesis")
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/CCLE_in_CCLE_out/logreg_on_LDA/CCLE_{measure}_loocv_'
                    f'classification_plot_LogReg.{file}', bbox_inches='tight', dpi=300)
    plt.close()
    # Plot for Circadian input
    loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_circadian.columns:
        for index_object in loocv_samples_circadian[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_circadian.index
    this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for d in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_circadian.loc[d]['corr_samples'].__str__():
                this_table_circadian.loc[d][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_circadian.loc[d]['corr_samples'].__str__() and sample not in \
                    loocv_samples_circadian.loc[d][
                        'not_corr_samples'].__str__():
                this_table_circadian.loc[d][sample] = np.nan
    # add accuracy
    this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
    sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig, axe = plt.subplots(1, 1, figsize=(6, 6))
    sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
    # set the np.nan tiles to grey
    axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    axe.set_yticklabels(y_labels, fontsize=12)
    # add frame lines around heatmap
    axe.axhline(y=0, color='k', linewidth=2)
    axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    axe.axvline(x=0, color='k', linewidth=2)
    axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    axe.xaxis.tick_top()  # x axis on top
    axe.xaxis.set_label_position('top')
    axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    axe.set_title(f"LogReg on LDA leave-one-out classification results by cell line on {measure} (BMAL1-PER2 mean)\n"
                  f"sample size in parenthesis")
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5/Circadian_in_CCLE_out/logreg_on_lda/circadian_{measure}_'
                    f'BMAL1_PER2_loocv_classification_plot_LogReg.{file}', bbox_inches='tight', dpi=300)
    plt.close()

# #################### FURTHER ALTERNATIVES: LDA vs LogReg on PCA-transformed 9 features and on top 3 obvious ##########

# ## First on the 9 features transformed by PCA ##

# On cisplatin sensitivity data

drug = 'Cisplatin'
tod_scores_for_loocv = list(data_sensitivty.keys())

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle_lda = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                           index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_circadian_lda = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                                index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_ccle_logreg = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                              index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_circadian_logreg = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                                   index=metrics, columns=tod_scores_for_loocv)

# Without sample MCF10A
for measure in list(data_sensitivty.keys()):
    print(f'* Analysing stat {measure}...')
    df = data_sensitivty[measure].T.copy()
    df = df.drop('MCF10A', axis=0)
    this_drug = df[drug]
    this_drug_end = this_drug.dropna()
    print(f'** Discretize drug sensitivity by median...')
    this_drug_binarized = this_drug_end > this_drug_end.median()
    # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
    print('**** Logistic regression and LDA on CCLE expression...')
    ccle_expression = ccle_data.loc[this_drug_end.index]
    # LogReg and LDA on PCA transformed LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogRegLDA_on_PCA_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                               target=drug, estimator='lda', pca_comp=3)
    loocv_results_grid_ccle_lda[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogRegLDA_on_PCA_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                               target=drug, estimator='logreg', pca_comp=3)
    loocv_results_grid_ccle_logreg[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # •	By circadian features mean
    print(f'**** Logistic regression and LDA on BMAL1_PER2 mean parameters...\n')
    mat = circadian_data.copy()
    matrix_circ = mat.dropna(axis=1, how='all')
    matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
    matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
    matrix_circ_norm = matrix_circ_norm.loc[this_drug_end.index]
    # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
    # sample from the drug to keep consistent data sets
    all_together = pd.concat([matrix_circ_norm, this_drug_binarized], axis=1)
    all_together_final = all_together.dropna(axis=0)
    matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
    this_drug_binarized_consistent = all_together_final[drug]
    # logreg and LDA on PCA transformed LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogRegLDA_on_PCA_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug, estimator='lda',
                               pca_comp=3)
    loocv_results_grid_circadian_lda[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogRegLDA_on_PCA_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug, estimator='logreg',
                               pca_comp=3)
    loocv_results_grid_circadian_logreg[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]

# Plot for CCLE input
for loocv_results_grid_ccle, name in zip([loocv_results_grid_ccle_lda, loocv_results_grid_ccle_logreg],
                                         ['lda', 'logreg']):
    loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_ccle.columns:
        for index_object in loocv_samples_ccle[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_ccle.index
    this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for drug in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_ccle.loc[drug]['corr_samples'].__str__():
                this_table_ccle.loc[drug][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_ccle.loc[drug]['corr_samples'].__str__() and sample not in \
                    loocv_samples_ccle.loc[drug]['not_corr_samples'].__str__():
                this_table_ccle.loc[drug][sample] = np.nan
    # add accuracy
    this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
    sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig, axe = plt.subplots(1, 1, figsize=(6, 6))
    sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
    # set the np.nan tiles to grey
    axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    axe.set_yticklabels(y_labels, rotation=0, fontsize=12)
    # add frame lines around heatmap
    axe.axhline(y=0, color='k', linewidth=2)
    axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    axe.axvline(x=0, color='k', linewidth=2)
    axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    axe.xaxis.tick_top()  # x axis on top
    axe.xaxis.set_label_position('top')
    axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    axe.set_title(f"{name.capitalize()} leave-one-out classification results on PCA data, Cisplatin (CCLE)\n"
                  f"sample size in parenthesis")
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5_alternatives/CCLE_on_cisplatin_PCA_transformed'
                    f'classification_plot_{name}.{file}', bbox_inches='tight', dpi=300)
    plt.close()

# Plot for Circadian input
for loocv_results_grid_circadian, name in zip([loocv_results_grid_circadian_lda, loocv_results_grid_circadian_logreg],
                                         ['lda', 'logreg']):
    loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_circadian.columns:
        for index_object in loocv_samples_circadian[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_circadian.index
    this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for drug in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_circadian.loc[drug]['corr_samples'].__str__():
                this_table_circadian.loc[drug][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_circadian.loc[drug]['corr_samples'].__str__() and sample not in \
                    loocv_samples_circadian.loc[drug]['not_corr_samples'].__str__():
                this_table_circadian.loc[drug][sample] = np.nan
    # add accuracy
    this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
    sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig, axe = plt.subplots(1, 1, figsize=(6, 6))
    sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
    # set the np.nan tiles to grey
    axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    axe.set_yticklabels(y_labels, rotation=0, fontsize=12)
    # add frame lines around heatmap
    axe.axhline(y=0, color='k', linewidth=2)
    axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    axe.axvline(x=0, color='k', linewidth=2)
    axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    axe.xaxis.tick_top()  # x axis on top
    axe.xaxis.set_label_position('top')
    axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    axe.set_title(f"{name.capitalize()} leave-one-out classification results on PCA data, Cisplatin (BMAL1-PER2 mean)\n"
                  f"sample size in parenthesis")
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5_alternatives/Circadian_Cisplatin_PCAtransformed_'
                    f'BMAL1_PER2_loocv_classification_plot_{name}.{file}', bbox_inches='tight', dpi=300)
    plt.close()

# Same on CCLE sensitivity data

drugs = remaining_drugs.copy()
tod_scores_for_loocv = ['EC50 (uM)', 'IC50 (uM)', 'ActArea']

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle_lda = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)
loocv_results_grid_circadian_lda = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)
loocv_results_grid_ccle_logreg = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)
loocv_results_grid_circadian_logreg = \
    pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)

# For each drug
for measure in tod_scores_for_loocv:
    print(f'* Analysing stat {measure}...')
    df = remaining_ccle_drugs.copy()
    df = df[['Primary Cell Line Name', 'Compound', measure]]
    for drug in remaining_drugs:
        print(f'*** Analysing for drug {drug}...')
        # get this drug only
        this_drug = df.loc[df['Compound'] == drug]
        this_drug_end = this_drug.dropna()
        # the number of samples must be more than the number of classes, so we always need at least 3 samples
        if len(this_drug_end) <= 3:
            # somtimes, IC50 has mostly the value 8, resulting in all classes being True / False
            text = f'only {len(this_drug_end)} cell line was' if len(this_drug_end) == 1 else \
                f'only {len(this_drug_end)} cell line were' if len(this_drug_end) > 1 else 'no cell lines were'
            print(f'For the drug {drug}, {text} measured for {measure}. This example is skipped for the LDA analysis.')
            continue  # skip if we only remain with 1 or 0 samples (most cases are EC50 (uM) not being measured)
        print(f'** Discretize drug sensitivity by median...')
        this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median()
        this_drug_binarized.name = drug
        if len(this_drug_binarized.unique()) == 1:  # apply the -0.001 trick only if necessary!
            this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median() - 0.001
            print('## Trying "median - 0.001" for discretiztation cutoff ##')
        if len(this_drug_binarized.unique()) == 1:
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in one '
                  f'single class. This example is skipped for the LDA analysis.')
            continue
        if 1 in set(this_drug_binarized.value_counts()):  # skip LDA loocv of samples where ther is only 1 true or false
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in two classes but '
                  f'one is only represented once. This example is skipped for the LDA LOOCV analysis.')
            continue
        subtypes_here = pd.DataFrame(
            [dict_to_map_line_and_subtype[x] for x in this_drug_end['Primary Cell Line Name']],
            index=this_drug_end['Primary Cell Line Name'], columns=['subtype'])
        this_drug_binarized.index = subtypes_here.index
        # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
        print('**** Logistic regression on CCLE expression...')
        ccle_expression = ccle_data.loc[this_drug_end['Primary Cell Line Name']]
        # LogReg on LDA transformed LOOCV machinery
        this_drug_binarized.name = drug  # make sure the targets have the right name
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogRegLDA_on_PCA_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                                target=drug, estimator='lda', pca_comp=3)
        loocv_results_grid_ccle_lda[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogRegLDA_on_PCA_loocv(data=ccle_expression, y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                                target=drug, estimator='logreg', pca_comp=3)
        loocv_results_grid_ccle_logreg[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
        # •	By circadian features mean
        print(f'**** Logistic regression on BMAL1_PER2 mean parameters...\n')
        mat = circadian_data.copy()
        matrix_circ = mat.dropna(axis=1, how='all')
        matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
        matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
        matrix_circ_norm = matrix_circ_norm.loc[subtypes_here.index]
        # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
        # sample from the drug to keep consistent data sets
        all_together = pd.concat([matrix_circ_norm, this_drug_binarized, subtypes_here], axis=1)
        all_together_final = all_together.dropna(axis=0)
        matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
        this_drug_binarized_consistent = all_together_final[drug]
        if 1 in this_drug_binarized_consistent.value_counts().values or len(
                all_together_final) <= 3:  # also here <=3 instead of 2
            print(f'/!\\ Logistic regression on circadian parameters targeting {drug} is only possible '
                  f'for\n{len(all_together_final)} samples where parameters and {measure} are available,'
                  f'\nbut the binarization of these remaining samples is left with 1 class only using median.'
                  f'\nLDA requires at least 2 classes and 3 samples. /!\\')
            continue
        # logreg on LDA transformed LOOCV machinery
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogRegLDA_on_PCA_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug, estimator='lda',
                                   pca_comp=3)
        loocv_results_grid_circadian_lda[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogRegLDA_on_PCA_loocv(data=matrix_circ_norm_final, y=all_together_final, target=drug, estimator='logreg',
                                   pca_comp=3)
        loocv_results_grid_circadian_logreg[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # plot at this level for each individual measure (target), 2 per measure (CCLE and circ (including channels))
    # Plot for CCLE sensitivity input
    for loocv_results_grid_ccle, name in zip([loocv_results_grid_ccle_lda, loocv_results_grid_ccle_logreg],
                                             ['lda', 'logreg']):
        loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
        columns = set()
        for col in loocv_samples_ccle.columns:
            for index_object in loocv_samples_ccle[col]:
                if isinstance(index_object, pd.Index):
                    columns.update(index_object)
        columns = list(columns)
        indeces = loocv_samples_ccle.index
        this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
        for d in indeces:
            for sample in columns:
                # case of correct samples and wrong samples
                if sample in loocv_samples_ccle.loc[d]['corr_samples'].__str__():
                    this_table_ccle.loc[d][sample] += 1
                # case of the N/A
                if sample not in loocv_samples_ccle.loc[d]['corr_samples'].__str__() and sample not in \
                        loocv_samples_ccle.loc[d][
                            'not_corr_samples'].__str__():
                    this_table_ccle.loc[d][sample] = np.nan
        # add accuracy
        this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
        sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
        sorted_table = sorted_table.dropna(axis=0, how='all')
        # table figure
        # Create a mask for the 10th column and nans
        nanmask = sorted_table.isna()
        mask = np.zeros_like(sorted_table, dtype=bool)
        mask[:, -1] = ~nanmask.iloc[:, -1]
        combinedmask = mask | nanmask
        fig, axe = plt.subplots(1, 1, figsize=(6, 6))
        sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                    annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
        # set the np.nan tiles to grey
        axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
                   extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
        for i in range(sorted_table.shape[0]):
            axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                     ha='center', va='center', color='black', fontsize=8)
        samples_per_drug = []
        for tod_drug in sorted_table.index:
            samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
        y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
        # set the labels plus sample size
        axe.set_yticklabels(y_labels, fontsize=12)
        # add frame lines around heatmap
        axe.axhline(y=0, color='k', linewidth=2)
        axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
        axe.axvline(x=0, color='k', linewidth=2)
        axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
        axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
        axe.xaxis.tick_top()  # x axis on top
        axe.xaxis.set_label_position('top')
        axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
        axe.set_title(f"{name.capitalize()} leave-one-out classification results on PCA data, {measure} (CCLE)\n"
                      f"sample size in parenthesis")
        for file in formats:
            plt.savefig(f'./results_revised_review/comment_5_alternatives/CCLE_{measure}_PCAtransformed_loocv_'
                        f'classification_plot_{name}.{file}', bbox_inches='tight', dpi=300)
        plt.close()
    # Plot for Circadian input
    for loocv_results_grid_circadian, name in zip([loocv_results_grid_circadian_lda,
                                                   loocv_results_grid_circadian_logreg], ['lda', 'logreg']):
        loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
        columns = set()
        for col in loocv_samples_circadian.columns:
            for index_object in loocv_samples_circadian[col]:
                if isinstance(index_object, pd.Index):
                    columns.update(index_object)
        columns = list(columns)
        indeces = loocv_samples_circadian.index
        this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
        for d in indeces:
            for sample in columns:
                # case of correct samples and wrong samples
                if sample in loocv_samples_circadian.loc[d]['corr_samples'].__str__():
                    this_table_circadian.loc[d][sample] += 1
                # case of the N/A
                if sample not in loocv_samples_circadian.loc[d]['corr_samples'].__str__() and sample not in \
                        loocv_samples_circadian.loc[d][
                            'not_corr_samples'].__str__():
                    this_table_circadian.loc[d][sample] = np.nan
        # add accuracy
        this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
        sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
        sorted_table = sorted_table.dropna(axis=0, how='all')
        # table figure
        # Create a mask for the 10th column and nans
        nanmask = sorted_table.isna()
        mask = np.zeros_like(sorted_table, dtype=bool)
        mask[:, -1] = ~nanmask.iloc[:, -1]
        combinedmask = mask | nanmask
        fig, axe = plt.subplots(1, 1, figsize=(6, 6))
        sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                    annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
        # set the np.nan tiles to grey
        axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
                   extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
        for i in range(sorted_table.shape[0]):
            axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                     ha='center', va='center', color='black', fontsize=8)
        samples_per_drug = []
        for tod_drug in sorted_table.index:
            samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
        y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
        # set the labels plus sample size
        axe.set_yticklabels(y_labels, fontsize=12)
        # add frame lines around heatmap
        axe.axhline(y=0, color='k', linewidth=2)
        axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
        axe.axvline(x=0, color='k', linewidth=2)
        axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
        axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
        axe.xaxis.tick_top()  # x axis on top
        axe.xaxis.set_label_position('top')
        axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
        axe.set_title(f"{name.capitalize()} leave-one-out classification results on PCA data, {measure} "
                      f"(BMAL1-PER2 mean)\nsample size in parenthesis")
        for file in formats:
            plt.savefig(f'./results_revised_review/comment_5_alternatives/circadian_{measure}_PCAtransformed'
                        f'BMAL1_PER2_loocv_classification_plot_{name}.{file}', bbox_inches='tight', dpi=300)
        plt.close()


# ## Then on the top 3 features we used in k means, without PCA transformation ##

# We know which one to take for the circadian parameters, but what about the gene expression? Let's see with PCA

#  - rerun PCA on standard scaled data (careful with arrow scaling) (+ stats if possible, if so also for the min-max)

width = 0.3  # for the barplots
best_components = 5  # enough to describe 95% of the variance

random.seed(seed)
np.random.seed(seed)

# ------- pca plot to determine good gene expression candidates
data_here = X_short
matrix = data_here.copy()
full_ind = []
for x in matrix.index:
    full_ind.append(df_matrix_bmal1_per2.index[[x in k for k in df_matrix_bmal1_per2.index]][0])
matrix.index = full_ind
df = matrix.dropna(axis=1, how='all')
df_norm = df.dropna(axis=0)  # need to have all possible nans removed
# Get pca of the data
pca = PCA(n_components=best_components, random_state=seed)
pca.fit(df_norm)
PCs = pca.fit_transform(df_norm)
PCdf = pd.DataFrame(data=PCs, columns=["PC" + str(i) for i in range(1, PCs.shape[1] + 1)])
PCdf.index = df_norm.index
# Match PC names to loadings
pc_loadings = dict(zip(PCdf.columns, pca.components_))
# Matrix of corr coefficients between pcs and features
loadings_df = pd.DataFrame.from_dict(pc_loadings)
loadings_df['feature_names'] = df_norm.columns
loadings_df = loadings_df.set_index('feature_names')

# plot the loading plot with the scatterplot
target_markers = {'TNBC-BL1': 'o',
                  'TNBC-BL2': '<',
                  'TNBC-M': 'd',
                  'BC-LumA': '^'}
target_colors = {'TNBC-BL1': 'brown',
                 'TNBC-BL2': 'orange',
                 'TNBC-M': 'green',
                 'BC-LumA': 'blue'}
xs = pca.components_[0]
ys = pca.components_[1]
targets = pd.DataFrame(index=PCdf.index, columns=['subtype'], data=[x.split('_')[0] for x in df_norm.index])

g = sns.lmplot(x='PC1', y='PC2', data=pd.concat([PCdf, targets], axis=1), fit_reg=False, hue='subtype',
               palette=target_colors, markers=list(target_markers.values()),
               scatter_kws={'linewidths': 1, 'edgecolor': 'k', "s": 75})
sns.move_legend(g, "upper right", bbox_to_anchor=(1.02, 0.98), frameon=True)
text_points = []
text_arrows = []
# label points
for ind in PCdf.index:
    text_points.append(plt.annotate(ind.split('_')[1], (PCdf.loc[ind]['PC1'], PCdf.loc[ind]['PC2']), fontsize=8,
                       ha='left', va='center', weight='bold'))
# label arrows
for i, varnames in enumerate(df_norm.columns):
    plt.arrow(0, 0, xs[i], ys[i], color='k', head_width=0.08, alpha=0.3, linestyle='-', linewidth=2)
    text_arrows.append(plt.text(xs[i], ys[i], varnames, fontsize=8))
plt.xlabel(f'PC1 ({"{:.1f}".format(round(pca.explained_variance_ratio_[0] * 100, 2))}%)')
plt.ylabel(f'PC2 ({"{:.1f}".format(round(pca.explained_variance_ratio_[1] * 100, 2))}%)')
plt.title(f'Circadian cell model mapping')
plt.tight_layout()
adjust_text(text_arrows + text_points)
for form in formats:
    plt.savefig(f'./results_revised_review/comment_5_alternatives/pca_on_gene_expression.{form}',
                bbox_inches='tight', dpi=dpi)
plt.close()

# PCA bars, here seed reset is not required as we continue to use transformed data from the biplot
original_hatch_width = 1.0
matplotlib.rcParams['hatch.linewidth'] = 8.0
fig, ax = plt.subplots(figsize=(8, 6))
loadings_df = loadings_df.reindex(loadings_df['PC1'].abs().sort_values(ascending=False).index)  # sort by absolute PC1
rect1 = ax.bar(np.arange(len(loadings_df.index)) - width/2, loadings_df["PC1"].abs(), width=width, label="PC1",
               color='steelblue', edgecolor="k", hatch=['\\' if i < 0 else None for i in loadings_df["PC1"]])
rect2 = ax.bar(np.arange(len(loadings_df.index)) + width/2,
               loadings_df["PC2"].abs(), width=width, label="PC2", edgecolor="k",
               color='lemonchiffon', hatch=['\\' if i < 0 else None for i in loadings_df["PC2"]])
ax.set_xticks(np.arange(len(loadings_df.index)))
ax.set_xticklabels(loadings_df.index)
ax.tick_params(axis="x", labelsize=10, labelrotation=45)
lgd = ax.legend(labels=['PC1', 'PC2'], fontsize=12, ncol=3)
handles, labs = lgd.axes.get_legend_handles_labels()
handles.append(Patch(facecolor='lightgrey', edgecolor='k', hatch='\\'))
labs.append('Negative')
lgd._legend_box = None
lgd._init_legend_box(handles, labs)
lgd._set_loc(lgd._loc)
lgd.set_title(lgd.get_title().get_text())
for num, ha in enumerate(lgd.legendHandles):
    if num < len(lgd.legendHandles) - 1:
        ha.set_hatch(None)
plt.xlabel('')
plt.ylabel('Absolute PC loading')
plt.title('Ranking of PC loadings')
plt.tight_layout()
for form in formats:
    plt.savefig(f'./results_revised_review/comment_5_alternatives/pca_gene_expression_bars.{form}',
                bbox_inches='tight', dpi=300)
plt.close()
matplotlib.rcParams['hatch.linewidth'] = original_hatch_width
# -------

# On cisplatin sensitivity data

drug = 'Cisplatin'
tod_scores_for_loocv = list(data_sensitivty.keys())

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

# obvious candidates in BMAL1_PER2
obvious_candidates_comb = ['period_median', 'phdiff_coeffvar', 'mra_circadian']
obvious_genes = ['PER1', 'CSNK1D', 'CLOCK']  # previously tested: PER1, BMAL1, CLOCK

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle_lda = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                           index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_circadian_lda = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                                index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_ccle_logreg = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                              index=metrics, columns=tod_scores_for_loocv)
loocv_results_grid_circadian_logreg = pd.DataFrame(data=np.zeros((len(metrics), len(tod_scores_for_loocv))),
                                                   index=metrics, columns=tod_scores_for_loocv)

# Without sample MCF10A
for measure in list(data_sensitivty.keys()):
    print(f'* Analysing stat {measure}...')
    df = data_sensitivty[measure].T.copy()
    df = df.drop('MCF10A', axis=0)
    this_drug = df[drug]
    this_drug_end = this_drug.dropna()
    print(f'** Discretize drug sensitivity by median...')
    this_drug_binarized = this_drug_end > this_drug_end.median()
    # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
    print('**** Logistic regression and LDA on CCLE expression...')
    ccle_expression = ccle_data.loc[this_drug_end.index]
    # LogReg and LDA on PCA transformed LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LDA_loocv(data=ccle_expression[obvious_genes], y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                  target=drug)
    loocv_results_grid_ccle_lda[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogReg_loocv(data=ccle_expression[obvious_genes], y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                     target=drug)
    loocv_results_grid_ccle_logreg[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # •	By circadian features mean
    print(f'**** Logistic regression and LDA on BMAL1_PER2 mean parameters...\n')
    mat = circadian_data.copy()
    matrix_circ = mat.dropna(axis=1, how='all')
    matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
    matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
    matrix_circ_norm = matrix_circ_norm.loc[this_drug_end.index]
    # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
    # sample from the drug to keep consistent data sets
    all_together = pd.concat([matrix_circ_norm, this_drug_binarized], axis=1)
    all_together_final = all_together.dropna(axis=0)
    matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
    this_drug_binarized_consistent = all_together_final[drug]
    # logreg and LDA on PCA transformed LOOCV machinery
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LDA_loocv(data=matrix_circ_norm_final[obvious_candidates_comb], y=all_together_final, target=drug)
    loocv_results_grid_circadian_lda[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
        LogReg_loocv(data=matrix_circ_norm_final[obvious_candidates_comb], y=all_together_final, target=drug)
    loocv_results_grid_circadian_logreg[measure] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]

# Plot for CCLE input
for loocv_results_grid_ccle, name in zip([loocv_results_grid_ccle_lda, loocv_results_grid_ccle_logreg],
                                         ['lda', 'logreg']):
    loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_ccle.columns:
        for index_object in loocv_samples_ccle[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_ccle.index
    this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for drug in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_ccle.loc[drug]['corr_samples'].__str__():
                this_table_ccle.loc[drug][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_ccle.loc[drug]['corr_samples'].__str__() and sample not in \
                    loocv_samples_ccle.loc[drug]['not_corr_samples'].__str__():
                this_table_ccle.loc[drug][sample] = np.nan
    # add accuracy
    this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
    sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig, axe = plt.subplots(1, 1, figsize=(6, 6))
    sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
    # set the np.nan tiles to grey
    axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    axe.set_yticklabels(y_labels, rotation=0, fontsize=12)
    # add frame lines around heatmap
    axe.axhline(y=0, color='k', linewidth=2)
    axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    axe.axvline(x=0, color='k', linewidth=2)
    axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    axe.xaxis.tick_top()  # x axis on top
    axe.xaxis.set_label_position('top')
    axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    axe.set_title(f"{name.capitalize()} leave-one-out classification results on obvious 3 candidates, Cisplatin "
                  f"(CCLE)\nsample size in parenthesis")
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5_alternatives/CCLE_on_cisplatin_obvious_3_'
                    f'classification_plot_{name}.{file}', bbox_inches='tight', dpi=300)
    plt.close()

# Plot for Circadian input
for loocv_results_grid_circadian, name in zip([loocv_results_grid_circadian_lda, loocv_results_grid_circadian_logreg],
                                         ['lda', 'logreg']):
    loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
    columns = set()
    for col in loocv_samples_circadian.columns:
        for index_object in loocv_samples_circadian[col]:
            if isinstance(index_object, pd.Index):
                columns.update(index_object)
    columns = list(columns)
    indeces = loocv_samples_circadian.index
    this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
    for drug in indeces:
        for sample in columns:
            # case of correct samples and wrong samples
            if sample in loocv_samples_circadian.loc[drug]['corr_samples'].__str__():
                this_table_circadian.loc[drug][sample] += 1
            # case of the N/A
            if sample not in loocv_samples_circadian.loc[drug]['corr_samples'].__str__() and sample not in \
                    loocv_samples_circadian.loc[drug]['not_corr_samples'].__str__():
                this_table_circadian.loc[drug][sample] = np.nan
    # add accuracy
    this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
    sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
    sorted_table = sorted_table.dropna(axis=0, how='all')
    # table figure
    # Create a mask for the 10th column and nans
    nanmask = sorted_table.isna()
    mask = np.zeros_like(sorted_table, dtype=bool)
    mask[:, -1] = ~nanmask.iloc[:, -1]
    combinedmask = mask | nanmask
    fig, axe = plt.subplots(1, 1, figsize=(6, 6))
    sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
    # set the np.nan tiles to grey
    axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
               extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
    for i in range(sorted_table.shape[0]):
        axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                 ha='center', va='center', color='black', fontsize=8)
    samples_per_drug = []
    for tod_drug in sorted_table.index:
        samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
    y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
    # set the labels plus sample size
    axe.set_yticklabels(y_labels, rotation=0, fontsize=12)
    # add frame lines around heatmap
    axe.axhline(y=0, color='k', linewidth=2)
    axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
    axe.axvline(x=0, color='k', linewidth=2)
    axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
    axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
    axe.xaxis.tick_top()  # x axis on top
    axe.xaxis.set_label_position('top')
    axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
    axe.set_title(f"{name.capitalize()} leave-one-out classification results on obvious 3 candidates, Cisplatin "
                  f"(BMAL1-PER2 mean)\nsample size in parenthesis")
    for file in formats:
        plt.savefig(f'./results_revised_review/comment_5_alternatives/Circadian_Cisplatin_'
                    f'BMAL1_PER2_obvious_3_loocv_classification_plot_{name}.{file}', bbox_inches='tight', dpi=300)
    plt.close()

# Same on CCLE sensitivity data

drugs = remaining_drugs.copy()
tod_scores_for_loocv = ['EC50 (uM)', 'IC50 (uM)', 'ActArea']

ccle_data = X_short.copy()
circadian_data = df_matrix_bmal1_per2.copy()

metrics = ['acc', 'bsl', 'accs', 'bsls', 'fbeta1', 'corr_samples', 'not_corr_samples']
loocv_results_grid_ccle_lda = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)
loocv_results_grid_circadian_lda = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)
loocv_results_grid_ccle_logreg = pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)
loocv_results_grid_circadian_logreg = \
    pd.DataFrame(data=np.zeros((len(metrics), len(drugs))), index=metrics, columns=drugs)

# For each drug
for measure in tod_scores_for_loocv:
    print(f'* Analysing stat {measure}...')
    df = remaining_ccle_drugs.copy()
    df = df[['Primary Cell Line Name', 'Compound', measure]]
    for drug in remaining_drugs:
        print(f'*** Analysing for drug {drug}...')
        # get this drug only
        this_drug = df.loc[df['Compound'] == drug]
        this_drug_end = this_drug.dropna()
        # the number of samples must be more than the number of classes, so we always need at least 3 samples
        if len(this_drug_end) <= 3:
            # somtimes, IC50 has mostly the value 8, resulting in all classes being True / False
            text = f'only {len(this_drug_end)} cell line was' if len(this_drug_end) == 1 else \
                f'only {len(this_drug_end)} cell line were' if len(this_drug_end) > 1 else 'no cell lines were'
            print(f'For the drug {drug}, {text} measured for {measure}. This example is skipped for the LDA analysis.')
            continue  # skip if we only remain with 1 or 0 samples (most cases are EC50 (uM) not being measured)
        print(f'** Discretize drug sensitivity by median...')
        this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median()
        this_drug_binarized.name = drug
        if len(this_drug_binarized.unique()) == 1:  # apply the -0.001 trick only if necessary!
            this_drug_binarized = this_drug_end[measure] > this_drug_end[measure].median() - 0.001
            print('## Trying "median - 0.001" for discretiztation cutoff ##')
        if len(this_drug_binarized.unique()) == 1:
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in one '
                  f'single class. This example is skipped for the LDA analysis.')
            continue
        if 1 in set(this_drug_binarized.value_counts()):  # skip LDA loocv of samples where ther is only 1 true or false
            print(f'For the drug {drug} measured for {measure}, the binarization by median resulted in two classes but '
                  f'one is only represented once. This example is skipped for the LDA LOOCV analysis.')
            continue
        subtypes_here = pd.DataFrame(
            [dict_to_map_line_and_subtype[x] for x in this_drug_end['Primary Cell Line Name']],
            index=this_drug_end['Primary Cell Line Name'], columns=['subtype'])
        this_drug_binarized.index = subtypes_here.index
        # •	By circadian clock gene expression levels (as done for ToD framework manuscript)
        print('**** Logistic regression on CCLE expression...')
        ccle_expression = ccle_data.loc[this_drug_end['Primary Cell Line Name']]
        # LogReg on LDA transformed LOOCV machinery
        this_drug_binarized.name = drug  # make sure the targets have the right name
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LDA_loocv(data=ccle_expression[obvious_genes], y=pd.concat([ccle_expression, this_drug_binarized], axis=1),
                      target=drug)
        loocv_results_grid_ccle_lda[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogReg_loocv(data=ccle_expression[obvious_genes],
                         y=pd.concat([ccle_expression, this_drug_binarized], axis=1), target=drug)
        loocv_results_grid_ccle_logreg[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
        # •	By circadian features mean
        print(f'**** Logistic regression on BMAL1_PER2 mean parameters...\n')
        mat = circadian_data.copy()
        matrix_circ = mat.dropna(axis=1, how='all')
        matrix_circ_norm = (matrix_circ - matrix_circ.min()) / (matrix_circ.max() - matrix_circ.min())
        matrix_circ_norm.index = [x.split('_')[1] for x in matrix_circ_norm.index]
        matrix_circ_norm = matrix_circ_norm.loc[subtypes_here.index]
        # Now we need to drop the samples in circadian parameters that have NaNs, and also remove that specific
        # sample from the drug to keep consistent data sets
        all_together = pd.concat([matrix_circ_norm, this_drug_binarized, subtypes_here], axis=1)
        all_together_final = all_together.dropna(axis=0)
        matrix_circ_norm_final = matrix_circ_norm.dropna(axis=0)
        this_drug_binarized_consistent = all_together_final[drug]
        if 1 in this_drug_binarized_consistent.value_counts().values or len(
                all_together_final) <= 3:  # also here <=3 instead of 2
            print(f'/!\\ Logistic regression on circadian parameters targeting {drug} is only possible '
                  f'for\n{len(all_together_final)} samples where parameters and {measure} are available,'
                  f'\nbut the binarization of these remaining samples is left with 1 class only using median.'
                  f'\nLDA requires at least 2 classes and 3 samples. /!\\')
            continue
        # logreg on LDA transformed LOOCV machinery
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LDA_loocv(data=matrix_circ_norm_final[obvious_candidates_comb], y=all_together_final, target=drug)
        loocv_results_grid_circadian_lda[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
        acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples = \
            LogReg_loocv(data=matrix_circ_norm_final[obvious_candidates_comb], y=all_together_final, target=drug)
        loocv_results_grid_circadian_logreg[drug] = [acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples]
    # plot at this level for each individual measure (target), 2 per measure (CCLE and circ (including channels))
    # Plot for CCLE sensitivity input
    for loocv_results_grid_ccle, name in zip([loocv_results_grid_ccle_lda, loocv_results_grid_ccle_logreg],
                                             ['lda', 'logreg']):
        loocv_samples_ccle = loocv_results_grid_ccle.T[['corr_samples', 'not_corr_samples']]
        columns = set()
        for col in loocv_samples_ccle.columns:
            for index_object in loocv_samples_ccle[col]:
                if isinstance(index_object, pd.Index):
                    columns.update(index_object)
        columns = list(columns)
        indeces = loocv_samples_ccle.index
        this_table_ccle = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
        for d in indeces:
            for sample in columns:
                # case of correct samples and wrong samples
                if sample in loocv_samples_ccle.loc[d]['corr_samples'].__str__():
                    this_table_ccle.loc[d][sample] += 1
                # case of the N/A
                if sample not in loocv_samples_ccle.loc[d]['corr_samples'].__str__() and sample not in \
                        loocv_samples_ccle.loc[d][
                            'not_corr_samples'].__str__():
                    this_table_ccle.loc[d][sample] = np.nan
        # add accuracy
        this_table_ccle['accuracy'] = np.mean(this_table_ccle, axis=1)
        sorted_table = this_table_ccle.sort_values(by='accuracy', ascending=False)
        sorted_table = sorted_table.dropna(axis=0, how='all')
        # table figure
        # Create a mask for the 10th column and nans
        nanmask = sorted_table.isna()
        mask = np.zeros_like(sorted_table, dtype=bool)
        mask[:, -1] = ~nanmask.iloc[:, -1]
        combinedmask = mask | nanmask
        fig, axe = plt.subplots(1, 1, figsize=(6, 6))
        sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                    annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
        # set the np.nan tiles to grey
        axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
                   extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
        for i in range(sorted_table.shape[0]):
            axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                     ha='center', va='center', color='black', fontsize=8)
        samples_per_drug = []
        for tod_drug in sorted_table.index:
            samples_per_drug.append(len(this_table_ccle.loc[tod_drug][:-1].dropna()))
        y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
        # set the labels plus sample size
        axe.set_yticklabels(y_labels, fontsize=12)
        # add frame lines around heatmap
        axe.axhline(y=0, color='k', linewidth=2)
        axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
        axe.axvline(x=0, color='k', linewidth=2)
        axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
        axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
        axe.xaxis.tick_top()  # x axis on top
        axe.xaxis.set_label_position('top')
        axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
        axe.set_title(f"{name.capitalize()} leave-one-out classification results on obvious 3 candidates, {measure} "
                      f"(CCLE)\nsample size in parenthesis")
        for file in formats:
            plt.savefig(f'./results_revised_review/comment_5_alternatives/CCLE_{measure}_obvious_3_loocv_'
                        f'classification_plot_{name}.{file}', bbox_inches='tight', dpi=300)
        plt.close()
    # Plot for Circadian input
    for loocv_results_grid_circadian, name in zip([loocv_results_grid_circadian_lda,
                                                   loocv_results_grid_circadian_logreg], ['lda', 'logreg']):
        loocv_samples_circadian = loocv_results_grid_circadian.T[['corr_samples', 'not_corr_samples']]
        columns = set()
        for col in loocv_samples_circadian.columns:
            for index_object in loocv_samples_circadian[col]:
                if isinstance(index_object, pd.Index):
                    columns.update(index_object)
        columns = list(columns)
        indeces = loocv_samples_circadian.index
        this_table_circadian = pd.DataFrame(data=np.zeros((len(indeces), len(columns))), index=indeces, columns=columns)
        for d in indeces:
            for sample in columns:
                # case of correct samples and wrong samples
                if sample in loocv_samples_circadian.loc[d]['corr_samples'].__str__():
                    this_table_circadian.loc[d][sample] += 1
                # case of the N/A
                if sample not in loocv_samples_circadian.loc[d]['corr_samples'].__str__() and sample not in \
                        loocv_samples_circadian.loc[d][
                            'not_corr_samples'].__str__():
                    this_table_circadian.loc[d][sample] = np.nan
        # add accuracy
        this_table_circadian['accuracy'] = np.mean(this_table_circadian, axis=1)
        sorted_table = this_table_circadian.sort_values(by='accuracy', ascending=False)
        sorted_table = sorted_table.dropna(axis=0, how='all')
        # table figure
        # Create a mask for the 10th column and nans
        nanmask = sorted_table.isna()
        mask = np.zeros_like(sorted_table, dtype=bool)
        mask[:, -1] = ~nanmask.iloc[:, -1]
        combinedmask = mask | nanmask
        fig, axe = plt.subplots(1, 1, figsize=(6, 6))
        sns.heatmap(sorted_table, cmap=plt.get_cmap('RdYlGn'), mask=combinedmask, cbar=False,
                    annot=False, vmin=0, center=.5, vmax=1, linewidth=.003, ax=axe)
        # set the np.nan tiles to grey
        axe.imshow(~nanmask, cmap='gray', vmin=0, vmax=1, aspect='auto',
                   extent=axe.get_xlim() + axe.get_ylim(), alpha=0.5)
        for i in range(sorted_table.shape[0]):
            axe.text(sorted_table.shape[1] - 0.5, i + 0.5, f'{sorted_table.iloc[i, -1]:.2f}',
                     ha='center', va='center', color='black', fontsize=8)
        samples_per_drug = []
        for tod_drug in sorted_table.index:
            samples_per_drug.append(len(this_table_circadian.loc[tod_drug][:-1].dropna()))
        y_labels = [x + f' ({s})' for x, s in zip(sorted_table.index, samples_per_drug)]
        # set the labels plus sample size
        axe.set_yticklabels(y_labels, fontsize=12)
        # add frame lines around heatmap
        axe.axhline(y=0, color='k', linewidth=2)
        axe.axhline(y=sorted_table.shape[0], color='k', linewidth=2)
        axe.axvline(x=0, color='k', linewidth=2)
        axe.axvline(x=sorted_table.iloc[:, :-1].shape[1], color='k', linewidth=1)
        axe.axvline(x=sorted_table.shape[1], color='k', linewidth=2)
        axe.xaxis.tick_top()  # x axis on top
        axe.xaxis.set_label_position('top')
        axe.set_xticklabels(sorted_table.columns, fontsize=12, rotation=75)
        axe.set_title(f"{name.capitalize()} leave-one-out classification results on obvious 3 candidates, {measure} "
                      f"(BMAL1-PER2 mean)\nsample size in parenthesis")
        for file in formats:
            plt.savefig(f'./results_revised_review/comment_5_alternatives/circadian_{measure}_obvious_3_'
                        f'BMAL1_PER2_loocv_classification_plot_{name}.{file}', bbox_inches='tight', dpi=300)
        plt.close()

########################################################################################################################
# END OF SCRIPT FOR ANALYSIS PIPELINE FOR THE SECOND PAPER COLLAB WITH CHARITE-BERLIN ##################################
########################################################################################################################
