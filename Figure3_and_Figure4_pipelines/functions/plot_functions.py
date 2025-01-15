# Functions for the analyses related to TNBC and circadian rythmicity
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import umap
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from adjustText import adjust_text
from sklearn.metrics import adjusted_rand_score


def grad_color_map(cluster_df, targets):
    """
    Function to generate gradient color map for continuous target vectors.

    Parameters
    ----------
    cluster_df : pandas.core.frame.DataFrame
        Feature matrix
    targets : list
        list of target values to be mapped to a fix color

    Returns
    -------
    colors : list
        list of mapped colors
    """
    # create gradient color map if continuous targets
    cmap = {cluster_df.index[m]: targets[m] for m in range(len(cluster_df))}  # mapping correct color to correct point
    sm = ScalarMappable(norm=Normalize(vmin=min(list(cmap.values())), vmax=max(list(cmap.values()))),
                        cmap=sns.cubehelix_palette(as_cmap=True))
    colors = [sm.to_rgba(cmap[obsv_id]) for obsv_id in cluster_df.index]
    return colors


def plot_PCA(X=None, y=None, label=None, comp=5, title='my_plot', seed=42):
    """
    Function to calculate and plot PCA and to colorize by multiple target labels.

    Parameters
    ----------
    X : pandas.core.frame.DataFrame
        Feature matrix
    y : pandas.core.frame.DataFrame
        Target matrix
    label : string
        Column name of target matrix to colorize points
    comp : int
        Number of PCA components, default 5
    title : string
        Figure title
   seed : int
        Random number generator seed for reproducibility
    """
    if X is not None and y is not None and label is not None:
        # fit PCA
        pca = PCA(n_components=comp, random_state=seed)  # we now we only need 14 here, the 15th is very close to zero
        pca.fit(X)
        PCs = pca.fit_transform(X)
        PCdf = pd.DataFrame(data=PCs, columns=["PC"+str(i) for i in range(1, PCs.shape[1]+1)])

        targets = list(y[label])
        all_colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']

        if len(np.unique(targets)) <= 7:
            colors = all_colors[:len(np.unique(targets))]
        else:
            # create gradient color map if continuous targets
            colors = grad_color_map(PCdf, targets)

        # draw the 4 PCA plots (PC1 vs PC2, PC2 vs PC3, PC1 vs PC3, % variance)
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
        for target, color in zip(np.unique(targets) if len(np.unique(targets)) <= 7 else targets, colors):
            idx = y[label] == target
            # pca plots
            ax1.scatter(PCdf.loc[idx.tolist(), 'PC1'], PCdf.loc[idx.tolist(), 'PC2'], color=color, s=50)
            ax1.set_title('PC1 v PC2', fontsize=10)
            ax1.set_xlabel(f'PC1 ({"{:.2f}".format(round(pca.explained_variance_ratio_[0] * 100, 2))}%)')
            ax1.set_ylabel(f'PC2 ({"{:.2f}".format(round(pca.explained_variance_ratio_[1] * 100, 2))}%)')
            ax2.scatter(PCdf.loc[idx.tolist(), 'PC2'], PCdf.loc[idx.tolist(), 'PC3'], color=color, s=50)
            ax2.set_title('PC2 v PC3', fontsize=10)
            ax2.set_xlabel(f'PC2 ({"{:.2f}".format(round(pca.explained_variance_ratio_[1] * 100, 2))}%)')
            ax2.set_ylabel(f'PC3 ({"{:.2f}".format(round(pca.explained_variance_ratio_[2] * 100, 2))}%)')
            ax3.scatter(PCdf.loc[idx.tolist(), 'PC1'], PCdf.loc[idx.tolist(), 'PC3'], color=color, s=50)
            ax3.set_title('PC1 v PC3', fontsize=10)
            ax3.set_xlabel(f'PC1 ({"{:.2f}".format(round(pca.explained_variance_ratio_[0] * 100, 2))}%)')
            ax3.set_ylabel(f'PC3 ({"{:.2f}".format(round(pca.explained_variance_ratio_[2] * 100, 2))}%)')
            # variance bar plot
            ax4.bar(range(1, PCs.shape[1] + 1), pca.explained_variance_ratio_ * 100, color='skyblue')
            for i in range(PCs.shape[1]):
                ax4.annotate(str("{:.2f}".format(round(pca.explained_variance_ratio_[i] * 100, 2))),
                             xy=(i + 1, pca.explained_variance_ratio_[i] * 100), ha='center', va='bottom',
                             size=8, weight='normal')
            ax4.set_title('Explained variance by principal components', fontsize=10)
            ax4.set_xlabel('Principal components')
            ax4.set_ylabel('Variance [%]')
            # control ticks of axis 4
            plt.sca(ax4)
            plt.xticks(range(1, PCs.shape[1] + 1))
            plt.suptitle(title, fontsize=14)
            if len(np.unique(targets)) <= 7:
                f.legend(np.unique(targets), loc='upper right', ncol=2, fontsize=7)  # other font size to not overlap with titles
            else:
                leg = f.legend([f'min: {"{:.2f}".format(round(min(targets), 2))}',
                                f'max: {"{:.2f}".format(round(max(targets), 2))}'],
                               labelcolor=[min(zip(targets, colors))[1], max(zip(targets, colors))[1]],
                               loc='upper right', ncol=1, fontsize=9)
                leg.legendHandles[0].set_color(min(zip(targets, colors))[1])
                leg.legendHandles[1].set_color(max(zip(targets, colors))[1])
            f.tight_layout()
    else:
        raise ValueError('Either X or y or label is not set.')


def plot_UMAP(X=None, y=None, label=None, n_neighbors=15, seed=42, title='my_plot', custom_color_marker=False, **kwargs):
    """
    Function to calculate and plot UMAP and to colorize by multiple target labels.

    Parameters
    ----------
    X : pandas.core.frame.DataFrame
        Feature matrix
    y : pandas.core.frame.DataFrame
        Target matrix
    label : string
        Column name of target matrix to colorize points
    n_neighbors : int
        Control how UMAP balances local versus global structure in the data
    seed : int
        Random number generator seed for reproducibility
    title : string
        Figure title
    """
    if X is not None and y is not None and label is not None:
        # fit UMAP
        reducer = umap.UMAP(n_neighbors=n_neighbors, random_state=seed)
        UMAPs = reducer.fit_transform(X)
        UMAPdf = pd.DataFrame(data=UMAPs, columns=["UMAP1", "UMAP2"])

        targets = list(y[label])
        # all_colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
        all_colors = ['r', 'm', 'b', 'orange', 'g', 'y', 'k']

        if len(np.unique(targets)) <= 7:
            colors = all_colors[:len(np.unique(targets))]
        else:
            # create gradient color map if continuous targets
            colors = grad_color_map(UMAPdf, targets)

        # draw the UMAP plot
        plt.figure(figsize=(6, 6))
        texts = []
        if custom_color_marker:
            for target in X.index:
                marker = '^' if 'BC-LumA' in target.split('_')[0] else 's' if 'Epithelial' in target.split('_')[0] else \
                    'o' if 'TNBC-BL1' in target.split('_')[0] else '<' if 'TNBC-BL2' in target.split('_')[0] else 'd'
                color = 'lightseagreen' if target.split('_')[-1] in ('MCF7', 'MCF10A', 'MDAMB468', 'T47D', 'HCC1937', 'HCC1143') else \
                    'purple' if target.split('_')[-1] in ('HCC1806', 'MDAMB231', 'CAL51') else \
                    'violet' if target.split('_')[-1] in ('HCC1428', 'MDAMB436', 'BT549', 'HCC70') else \
                    'yellow'
                idx = y[label].index == target
                plt.scatter(UMAPdf.loc[idx.tolist(), 'UMAP1'], UMAPdf.loc[idx.tolist(), 'UMAP2'], color=color, s=75,
                            marker=marker, edgecolor='k')
                for xx, yy, ss in zip(UMAPdf.loc[idx.tolist(), 'UMAP1'], UMAPdf.loc[idx.tolist(), 'UMAP2'],
                                      X.loc[idx.tolist()].index):
                    texts.append(plt.text(xx, yy, ss.split('_')[1],
                                          fontsize=8, fontweight='bold', ha='center', va='bottom', color='k'))
            # marker legends
            handles_marker = []
            handles_colors = []
            for m, clas in zip(['o', '<', 'd', '^', 's'], ['TNBC-BL1', 'TNBC-BL2', 'TNBC-M', 'LumA', 'Epithelial']):
                handles_marker.append(Line2D([0], [0], label=clas, marker=m,
                                             markersize=10, markeredgecolor='k', markerfacecolor='lightgrey', linestyle=''))
            for c, lab in zip(['lightseagreen', 'purple', 'violet', 'yellow'],
                              ['Functional', 'Unstable', 'Weak', 'Dysfunctional']):
                handles_colors.append(Patch(color=c, label=lab))
            plt.legend(handles=handles_marker + handles_colors, loc='best', ncol=1, fontsize=9)
        else:
            for target, color in zip(np.unique(targets) if len(np.unique(targets)) <= 7 else targets, colors):
                idx = y[label] == target
                marker = 'o'
                color = color
                plt.scatter(UMAPdf.loc[idx.tolist(), 'UMAP1'], UMAPdf.loc[idx.tolist(), 'UMAP2'], color=color, s=75,
                            marker=marker, edgecolor='k')
                for xx, yy, ss in zip(UMAPdf.loc[idx.tolist(), 'UMAP1'], UMAPdf.loc[idx.tolist(), 'UMAP2'],
                                      X.loc[idx.tolist()].index):
                    texts.append(plt.text(xx, yy, ss.split('_')[1],
                                          fontsize=8, fontweight='bold', ha='center', va='bottom', color='k'))
            if len(np.unique(targets)) <= 7:
                plt.legend(np.unique(targets), loc='best', ncol=2, fontsize=9)
            else:
                # in this case, we weirdly have to create the legend handles ourselves, else it will only yield 1 handle
                mini = Line2D([], [], color=min(zip(targets, colors))[1], marker='o', ls='',
                              label=f'min: {"{:.2f}".format(round(min(targets), 2))}')
                maxi = Line2D([], [], color=max(zip(targets, colors))[1], marker='o', ls='',
                              label=f'min: {"{:.2f}".format(round(max(targets), 2))}')
                plt.legend(handles=[mini, maxi],
                           labelcolor=[min(zip(targets, colors))[1], max(zip(targets, colors))[1]],
                           loc='best', fontsize=9)
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.xticks([])
        plt.yticks([])
        plt.title(title, fontsize=14)
        plt.tight_layout()
        adjust_text(texts, **kwargs)
    else:
        raise ValueError('Either X or y or label is not set.')


def plot_tSNE(X=None, y=None, label=None, seed=42, title='my_plot', **kwargs):
    """
    Function to calculate and plot t-SNE and to colorize by multiple target labels.

    Parameters
    ----------
    X : pandas.core.frame.DataFrame
        Feature matrix
    y : pandas.core.frame.DataFrame
        Target matrix
    label : string
        Column name of target matrix to colorize points
    seed : int
        Random number generator seed for reproducibility
    title : string
        Figure title
    """
    if X is not None and y is not None and label is not None:
        # fit t-SNE
        tsne = TSNE(random_state=seed, **kwargs)  # default 2 components
        tsne_res = tsne.fit_transform(X)
        tsne_df = pd.DataFrame(data=tsne_res, columns=["t-SNE1", "t-SNE2"])

        targets = list(y[label])
        # all_colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
        all_colors = ['r', 'm', 'b', 'orange', 'g', 'y', 'k']

        if len(np.unique(targets)) <= 7:
            colors = all_colors[:len(np.unique(targets))]
        else:
            # create gradient color map if continuous targets
            colors = grad_color_map(tsne_df, targets)
        # draw the tSNE plot
        plt.figure(figsize=(6, 6))
        texts = []
        for target, color in zip(np.unique(targets) if len(np.unique(targets)) <= 7 else targets, colors):
            idx = y[label] == target
            # tSNE plot
            plt.scatter(tsne_df.loc[idx.tolist(), 't-SNE1'], tsne_df.loc[idx.tolist(), 't-SNE2'], color=color, s=50)
            for xx, yy, ss in zip(tsne_df.loc[idx.tolist(), 't-SNE1'], tsne_df.loc[idx.tolist(), 't-SNE2'],
                                  X.loc[idx.tolist()].index):
                texts.append(plt.text(xx, yy, ss.split('_')[1],
                                      fontsize=8, fontweight='bold', ha='center', va='bottom', color='k'))
            plt.xlabel('t-SNE 1')
            plt.ylabel('t-SNE 2')
            plt.title(title, fontsize=14)
            if len(np.unique(targets)) <= 7:
                plt.legend(np.unique(targets), loc='best', ncol=2, fontsize=9)
            else:
                # in this case, we weirdly have to create the legend handles ourselves, else it will only yield 1 handle
                mini = Line2D([], [], color=min(zip(targets, colors))[1], marker='o', ls='',
                              label=f'min: {"{:.2f}".format(round(min(targets), 2))}')
                maxi = Line2D([], [], color=max(zip(targets, colors))[1], marker='o', ls='',
                              label=f'min: {"{:.2f}".format(round(max(targets), 2))}')
                plt.legend(handles=[mini, maxi],
                           labelcolor=[min(zip(targets, colors))[1], max(zip(targets, colors))[1]],
                           loc='best', fontsize=9)
            plt.tight_layout()
            adjust_text(texts)
    else:
        raise ValueError('Either X or y or label is not set.')


def plot_LDA(X=None, y=None, label=None, title='my_plot', lda_output_target=None, seed=42, annot_subtype=False,
             annot_cellline=False, data_info=None, return_loading=False, custom_color_marker=False, leg_labels=None,
             **kwargs):
    """
    Function to calculate and plot LDA and to colorize by multiple target labels.

    Parameters
    ----------
    X : pandas.core.frame.DataFrame
        Feature matrix
    y : pandas.core.frame.DataFrame
        Target matrix
    label : string
        Column name of target matrix to colorize points
    title : string
        Figure title
    lda_output_target : string
        LDA-specific output target to control on what feature the supervised LDA is fitted, as it would not work for
        continuous features
    seed : int
        Random number generator seed for reproducibility
    annot_subtype : bool
        Indication whether or not labels should be added to a 1-D LD plot
    annot_cellline : bool
        Indication whether or not cell line should be added beneath the subtype, currently only works if annot_subtype
        is True and requires data_info with the cell line name mapped to the DepMap ID.
    data_info : pandas.core.frame.DataFrame
        Dataframe containing the cell line names and DepMap IDs of each cell line.
    return_loading : bool
        If normalized loadings should be returned
    **kwargs : kwargs
        Further arguments for adjusted_text()
    """
    # got quite some fixes for manuscript2
    if X is not None and y is not None and label is not None:
        # fit LDA as supervised or unsupervised classifier!
        comp = len(y[lda_output_target].unique()) - 1
        LDA = LinearDiscriminantAnalysis(n_components=comp)
        # lda components = number of classes - 1 (only for subtypes)
        LDAs = LDA.fit_transform(X, y[lda_output_target])  # Make sure LDA is only fit to categorical subtypes
        LDAdf = pd.DataFrame(data=LDAs, columns=[f'LD{number + 1}' for number in np.arange(comp)])
        # If LDA only results in 1 linear discriminant dimension, (which is the case in this study as we are looking for
        # binarized targets, therefore we only get 1 linear discriminant dimension; except if the target is Subtype, in
        # that case we will analyse differenly) we need to make sure that not only the first feature
        # is reported, but also those that may have similar loading scores (c.f. PCA summed weights). Therefore we will
        # first normalize the lda loadings in case of 1 component (summed weight if multiple components).
        if comp == 1:
            LDAloadings = LDA.coef_
            normalized_lda_loading = abs(LDAloadings[0]) / sum(abs(LDAloadings[0]))
            sorted_normalized_idx = np.argsort(normalized_lda_loading)[::-1]
        else:
            LDAloadings, normalized_lda_loading, sorted_normalized_idx = [None] * 3

        targets = list(y[label])
        all_colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
        marker = 'o'
        if custom_color_marker:
            all_colors = ['teal', 'orange']
        if len(np.unique(targets)) <= 7:
            colors = all_colors[:len(np.unique(targets))]
        else:
            # create gradient color map if continuous targets
            colors = grad_color_map(LDAdf, targets)

        if comp >= 3:
            # draw the 4 LDA plots (LD1 vs LD2, LD2 vs LD3, LD1 vs LD3, % variance)
            f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
            for target, color in zip(np.unique(targets) if len(np.unique(targets)) <= 7 else targets, colors):
                idx = y[label] == target
                # pca plots
                ax1.scatter(LDAdf.loc[idx.tolist(), 'LD1'], LDAdf.loc[idx.tolist(), 'LD2'], color=color, s=75,
                            marker=marker, edgecolor='k')
                ax1.set_title('LD1 v LD2', fontsize=10)
                ax1.set_xlabel(f'LD1 ({"{:.2f}".format(round(LDA.explained_variance_ratio_[0] * 100, 2))}%)')
                ax1.set_ylabel(f'LD2 ({"{:.2f}".format(round(LDA.explained_variance_ratio_[1] * 100, 2))}%)')
                ax2.scatter(LDAdf.loc[idx.tolist(), 'LD2'], LDAdf.loc[idx.tolist(), 'LD3'], color=color, s=75,
                            marker=marker, edgecolor='k')
                ax2.set_title('LD2 v LD3', fontsize=10)
                ax2.set_xlabel(f'LD2 ({"{:.2f}".format(round(LDA.explained_variance_ratio_[1] * 100, 2))}%)')
                ax2.set_ylabel(f'LD3 ({"{:.2f}".format(round(LDA.explained_variance_ratio_[2] * 100, 2))}%)')
                ax3.scatter(LDAdf.loc[idx.tolist(), 'LD1'], LDAdf.loc[idx.tolist(), 'LD3'], color=color, s=75,
                            marker=marker, edgecolor='k')
                ax3.set_title('LD1 v LD3', fontsize=10)
                ax3.set_xlabel(f'LD1 ({"{:.2f}".format(round(LDA.explained_variance_ratio_[0] * 100, 2))}%)')
                ax3.set_ylabel(f'LD3 ({"{:.2f}".format(round(LDA.explained_variance_ratio_[2] * 100, 2))}%)')
            # variance bar plot
            ax4.bar(range(1, LDAs.shape[1] + 1), LDA.explained_variance_ratio_ * 100, color='skyblue')
            for i in range(LDAs.shape[1]):
                ax4.annotate(str("{:.2f}".format(round(LDA.explained_variance_ratio_[i] * 100, 2))),
                             xy=(i + 1, LDA.explained_variance_ratio_[i] * 100), ha='center', va='bottom',
                             size=8, weight='normal')
            ax4.set_title('Explained variance by linear discriminant components', fontsize=10)
            ax4.set_xlabel('Discriminant components')
            ax4.set_ylabel('Variance [%]')
            plt.sca(ax2)
            plt.suptitle(title, fontsize=14)
            if len(np.unique(targets)) <= 7:
                f.legend(np.unique(targets), loc='upper right' if len(np.unique(targets)) > 2 else 'upper left',
                         ncol=1 if len(np.unique(targets)) == 2 else 2, fontsize=7 if len(np.unique(targets)) > 2 else 9)
            else:
                leg = f.legend([f'min: {"{:.2f}".format(round(min(targets), 2))}',
                                f'max: {"{:.2f}".format(round(max(targets), 2))}'],
                               labelcolor=[min(zip(targets, colors))[1], max(zip(targets, colors))[1]],
                               loc='upper left', ncol=1, fontsize=9)
                leg.legendHandles[0].set_color(min(zip(targets, colors))[1])
                leg.legendHandles[1].set_color(max(zip(targets, colors))[1])
            f.tight_layout()
        else:
            # draw only 1 LDA in case of binarized feature
            np.random.seed(seed)
            jittered_y = pd.DataFrame(np.zeros(len(LDAdf)) + 0.1 * np.random.rand(len(LDAdf)) - 0.05,
                                      columns=[f'LD{number + 1}' for number in np.arange(comp)])
            f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
            for target, color in zip(np.unique(targets) if len(np.unique(targets)) <= 7 else targets, colors):
                idx = y[label] == target
                if custom_color_marker:
                    for sample in LDAdf.loc[idx.tolist(), 'LD1'].index:
                        marker = '*' if X.index[sample] in ('HCC1937', 'MCF7', 'MDAMB468', 'HCC1143', 'T47D') else \
                            '>' if X.index[sample] in ('MDAMB231', 'HCC1806', 'CAL51') else \
                            'v' if X.index[sample] in ('MDAMB436', 'HCC70', 'HCC1428', 'BT549') else \
                            'x'
                        s = 80 if marker in ('v', '>', 'x') else 120
                        ax1.scatter(LDAdf.loc[sample, 'LD1'], jittered_y.loc[sample, 'LD1'],  color=color,
                                    s=s, alpha=0.75, marker=marker, edgecolor='k' if marker != 'x' else None)
                else:
                    ax1.scatter(LDAdf.loc[idx.tolist(), 'LD1'], jittered_y.loc[idx.tolist(), 'LD1'], color=color,
                                s=75, alpha=0.75, marker=marker, edgecolor='k')
                ax1.set_title('First LD component', fontsize=10)
                ax1.set_xlabel(f'LD1')
                ax1.set_ylabel('')
                ax1.set_ylim([-0.1, 0.1])
                ax1.set_yticks([])
            if annot_subtype:
                texts = []
                # add to each target its belonging subtype
                for num, txt in enumerate(y['subtype']):
                    if annot_cellline:
                        txt = f"{data_info['cell_line_name'][data_info['stripped_cell_line_name'] == X.index[num]][0]}"
                    texts.append(ax1.text(LDAdf['LD1'][num], jittered_y['LD1'][num], txt, fontsize=8,
                                          fontweight='bold', ha='center', va='center'))
            # lda loading plot instead of variance plot
            ax2.bar(range(LDAloadings.shape[1])[:20], normalized_lda_loading[sorted_normalized_idx][:20] * 100,
                    color='lavender', edgecolor='k', linewidth=1)
            for i in range(LDAloadings.shape[1])[:20]:
                ax2.annotate(str("{:.2f}".format(round(normalized_lda_loading[sorted_normalized_idx][i] * 100, 3))),
                             xy=(i, normalized_lda_loading[sorted_normalized_idx][i] * 100), ha='center',
                             va='bottom', size=8, weight='normal')
                ax2.set_title(f'Top {len(range(LDAloadings.shape[1])[:20])} feature loadings by the linear discriminant '
                              f'component', fontsize=10)
                ax2.set_xlabel(None)
                ax2.set_ylabel('Distinctive contribution (%)')
            # control ticks of axis 2
            ax2.set_xticks(range(LDAloadings.shape[1])[:20])
            ax2.set_xticklabels(X.columns[sorted_normalized_idx][:20], rotation=90)
            plt.sca(ax2)
        plt.suptitle(title, fontsize=14)
        if len(np.unique(targets)) <= 7:
            if custom_color_marker:
                legend_elements = [plt.plot([], marker="", ls="")[0]]
                for c, typ in zip(colors, ['cytotoxic', 'cytostatic']):
                    legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=8,
                                                  label=typ))
                l = ax1.legend(handles=legend_elements, labels=leg_labels,
                               loc='upper left', ncol=3, fontsize=9, frameon=True, columnspacing=0.3, handletextpad=0.1)
                for vpack in l._legend_handle_box.get_children()[:1]:
                    for hpack in vpack.get_children():
                        hpack.get_children()[0].set_width(0)
                legend_elements = []
                for m, typ in zip(['*', '>', 'v', 'x'], ['Functional', 'Unstable', 'Weak', 'Dysfunctional']):
                    legend_elements.append(Line2D([0], [0], marker=m, color='w', markerfacecolor='lightgrey',
                                                  markersize=12 if m == '*' else 8, label=typ,
                                                  markeredgecolor='k' if marker != 'x' else None))
                ax1.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.07, -0.2),
                           title='$\\bf{Clock\\ phenotype}$', ncol=4, fontsize=9, frameon=False, columnspacing=0.3,
                           handletextpad=0.1)
                ax1.add_artist(l)
            else:
                f.legend(np.unique(targets), loc='upper right' if len(np.unique(targets)) > 2 else 'upper left',
                         ncol=1 if len(np.unique(targets)) == 2 else 2, fontsize=7 if len(np.unique(targets)) > 2 else 9)
        else:
            leg = f.legend([f'min: {"{:.2f}".format(round(min(targets), 2))}',
                            f'max: {"{:.2f}".format(round(max(targets), 2))}'],
                           labelcolor=[min(zip(targets, colors))[1], max(zip(targets, colors))[1]],
                           loc='upper left', ncol=1, fontsize=9)
            leg.legendHandles[0].set_color(min(zip(targets, colors))[1])
            leg.legendHandles[1].set_color(max(zip(targets, colors))[1])
        f.tight_layout()
        if annot_subtype:
            adjust_text(texts, ax=ax1, **kwargs, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
            # adjust_text must be used at the very end, else you will spend hours to figure this out...
        if return_loading:
            normalized_lda_loading_to_return = pd.DataFrame(normalized_lda_loading, index=X.columns,
                                                            columns=[lda_output_target])
            return normalized_lda_loading_to_return.T
    else:
        raise ValueError('Either X or y or label is not set.')


def cluster_elbow_curve(data, k, transformer, truth, seed):
    """
    Function to calculate and plot elbow curve, within cluster sum of squares and silhouette score.

    Parameters
    ----------
    data : pandas.core.frame.DataFrame
        Transformed data matrix
    k : int
        number of clusters to check for
    transformer : str
        name of transformer used (PCA, LDA, ...)
    truth : pandas.core.frame.DataFrame
        dataframe of ground truth, must be factorized
    """
    wcss = []
    silhouette = [0]
    adj_rand_score = [0]
    for i in range(1, k):
        kmeans = KMeans(n_clusters=i, init='k-means++', random_state=seed)
        kmeans.fit(data)
        wcss.append(kmeans.inertia_)
        preds = kmeans.fit_predict(data)
        adj_rand_score.append(adjusted_rand_score(truth, preds))
        if i > 1:
            silhouette.append(silhouette_score(data, preds))
    plt.figure(figsize=(10, 5))
    plt.plot(range(1, k), wcss)
    for i, z in enumerate(wcss):
        if i > 1:
            plt.text(i, wcss[i-1], "Silh.: %.2f\nARI: %.2f" % (silhouette[i], adj_rand_score[i]), fontsize=8)
    plt.title(f'{transformer}\nElbow Curve, silhouette scores and adjusted Rand Index for k > 1')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Within Cluster Sum of Squares')
    plt.xticks(range(1, k), range(1, k))
    plt.show()


def plot_3D(X, y=None, tar=None, transformer=None, df_set_name=None,n_components=5, seed=42):
    if transformer == 'lda':
        # fit LDA as supervised or unsupervised classifier!
        comp = len(y[tar].unique()) - 1
        trans = LinearDiscriminantAnalysis(n_components=comp)
        # lda components = number of classes - 1 (only for subtypes)
        LDAs = trans.fit_transform(X, y[tar])  # Make sure LDA is only fit to categorical subtypes
        trans_df = pd.DataFrame(data=LDAs, columns=[f'LD{number + 1}' for number in np.arange(comp)])
    elif transformer == 'pca':
        # fit PCA
        trans = PCA(n_components=n_components, random_state=seed)  # we now we only need 14 here, the 15th is very close to zero
        trans.fit(X)
        PCs = trans.fit_transform(X)
        trans_df = pd.DataFrame(data=PCs, columns=["PC"+str(i) for i in range(1, PCs.shape[1]+1)])

    targets = list(y[tar])
    # all_colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
    all_colors = ['r', 'm', 'b', 'orange', 'g', 'y', 'k']
    if len(np.unique(targets)) <= 7:
        colors = all_colors[:len(np.unique(targets))]
    else:
        # create gradient color map if continuous targets
        colors = grad_color_map(trans_df, targets)

    # draw the 3D plpot and % variances
    f = plt.figure(figsize=(12, 6))
    ax1 = f.add_subplot(1, 2, 1, projection='3d')
    ax2 = f.add_subplot(1, 2, 2)
    list_of_legs = []
    for target, color in zip(np.unique(targets) if len(np.unique(targets)) <= 7 else targets, colors):
        idx = y[tar] == target
        # lda plots
        xx, yy, zz = trans_df.loc[idx.tolist(), trans_df.columns[0]],\
                     trans_df.loc[idx.tolist(), trans_df.columns[1]],\
                     trans_df.loc[idx.tolist(), trans_df.columns[2]]
        g = ax1.scatter(xx, yy, zz, color=color, s=75)
        list_of_legs.append(g)
        # lines
        z2 = np.ones(shape=xx.shape) * min(trans_df[trans_df.columns[2]])
        for i, j, k, h, ind in zip(xx, yy, zz, z2, X.loc[idx.tolist(), X.columns[0]].index):
            ax1.plot([i, i], [j, j], [k, h], color='grey', linestyle='dashed', alpha=.75)
        for num, (i, j, k, ind) in enumerate(zip(xx, yy, zz, X.loc[idx.tolist(), X.columns[0]].index)):
            txt = ind.split('_')[1] if '_' in ind else ind
            ax1.text(i, j * 1.05, k * 1.075 if num // 2 == 0 else k - (0.075 * k),
                       txt, fontsize=7, fontweight='bold', ha='center', va='top', color='k', zorder=10)
    ax1.set_title(f'3D {"Linear Discriminant Analysis" if transformer == "lda" else "Principal Component Analysis Plot"}', fontsize=14)
    ax1.set_xlabel(f'{trans_df.columns[0]} ({"{:.2f}".format(round(trans.explained_variance_ratio_[0] * 100, 2))}%)')
    ax1.set_ylabel(f'{trans_df.columns[1]} ({"{:.2f}".format(round(trans.explained_variance_ratio_[1] * 100, 2))}%)')
    ax1.set_zlabel(f'{trans_df.columns[2]} ({"{:.2f}".format(round(trans.explained_variance_ratio_[2] * 100, 2))}%)')
    ax1.legend(list_of_legs, np.unique(targets), loc='best', ncol=3)
    # variance bar plot
    ax2.bar(range(1, LDAs.shape[1] + 1 if transformer == 'lda' else PCs.shape[1] + 1),
            trans.explained_variance_ratio_ * 100, color='skyblue')
    for i in range(LDAs.shape[1] if transformer == 'lda' else PCs.shape[1]):
        ax2.annotate(str("{:.2f}".format(round(trans.explained_variance_ratio_[i] * 100, 2))),
                     xy=(i + 1, trans.explained_variance_ratio_[i] * 100), ha='center', va='bottom',
                     size=8, weight='normal')
    ax2.set_title(f'Explained variance by {transformer.upper()} components', fontsize=14)
    ax2.set_xlabel(f'{"Principal" if transformer == "pca" else "Discriminant"} components')
    ax2.set_ylabel('Variance [%]')
    f.suptitle(f"{transformer.upper()} on {df_set_name}, colored by {tar}", fontsize=16)
    plt.show()
