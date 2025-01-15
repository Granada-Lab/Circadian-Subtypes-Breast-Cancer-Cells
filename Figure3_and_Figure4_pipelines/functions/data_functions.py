import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import scipy.stats as ss
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import make_scorer, f1_score


def reformat_drugs(dataframe):
    """reshapes a CCLE pandas dataframe from 'one line per datapoint' to a more convenient
    'one line per sample' format, meaning the response of a given cell line to different drugs
    will be placed on the same line in different columns."""

    drug_names = dataframe["Compound"].unique()
    dataframe["Compound"].value_counts()
    # concatenate the drug info with one line per cell line
    merged = pd.DataFrame()
    for thisDrug in drug_names:
        dataframe_spec = dataframe.loc[dataframe["Compound"] == thisDrug]
        dataframe_spec_clean = dataframe_spec.drop(
            columns=[
                "Primary Cell Line Name",
                "Compound",
                "Target",
                "Activity SD",
                "Num Data",
                "FitType",
            ]
        )
        dataframe_spec_clean.columns = [
            "CCLE Cell Line Name",
            thisDrug + "_dr_doses",
            thisDrug + "_dr_responses",
            thisDrug + "_EC50",
            thisDrug + "_IC50",
            thisDrug + "_Amax",
            thisDrug + "_ActArea",
        ]

        if merged.empty:
            merged = dataframe_spec_clean.copy()
        else:
            merged = pd.merge(
                merged,
                dataframe_spec_clean,
                how="left",
                on="CCLE Cell Line Name",
                sort=False,
                suffixes=("_x", "_y"),
                copy=True,
            )
    merged_dataframe = merged.set_index("CCLE Cell Line Name")

    actarea_cols = [x for x in merged_dataframe.columns if 'ActArea' in x]

    return merged_dataframe[actarea_cols]


def remove_columns_with_zeros(df, max_zeros):
    """
    Removes columns that contain more than a certain number of zeros
    Made by ChatGPT
    :param df: some DataFrame
    :param max_zeros: an arbitrary threshold
    :return: a much cleaner DataFrame
    """
    # Calculating the number of zeros in each column
    num_zeros = (df == 0).sum(axis=0)
    # Selecting columns with less than or equal to max_zeros zeros
    df = df.loc[:, num_zeros <= max_zeros]

    return df


def binarize_with_kmeans(df, n_clusters=2, keep_order=True, random_state=None):
    """
    Binarize a DataFrame using k-means clustering independently for each column.
    df : pandas DataFrame
    n_clusters : int, optional
        The number of clusters to use for k-means clustering. Default is 2.

    Returns:
    --------
    pandas DataFrame
        The binarized DataFrame with the same shape as the input DataFrame.
    """
    df_binary = pd.DataFrame(index=df.index)

    for column in df.columns:
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
        labels = kmeans.fit_predict(df[column][df[column].notnull()].values.reshape(-1, 1))
        centroids = kmeans.cluster_centers_
        is_inverted = centroids[0] > centroids[1]
        if is_inverted and keep_order:
            labels = np.logical_not(labels).astype(int)
        if df[column].isnull().values.any():
            # suppress SettingWithCopyWarning warning
            pd.options.mode.chained_assignment = None
            nan_mask = df[column].isna()
            df_binary[column] = nan_mask
            # put nan back where they were
            df_binary[column][nan_mask] = np.nan
            # on inverted nan location, add the labels
            df_binary[column][~nan_mask] = labels
        else:
            df_binary[column] = labels
    return df_binary


# Corrected Cramer's V correlation between categorical features
def cramers_corrected_stat(x, y):
    """
    Function to calculate corrected Cramers V statistic for categorical-categorical association. Uses correction
    from Bergsma and Wicher, Journal of the Korean Statistical Society 42 (2013): 323-328.

    Parameters
    ----------
    x : np.array
        array of first vector or column to analyze Cramers V correlation with the second
    y : np.array
        array of second vector or column to analyze Cramers V correlation with the first

    Returns
    -------
    result : float
        float value of the corrected Cramers V correlation coefficient between x and y
    """
    result = -1
    if len(np.unique(x)) == 1:
        print("First variable is constant")
    elif len(np.unique(y)) == 1:
        print("Target variable is constant")
    else:
        conf_matrix = pd.crosstab(x, y)
        if conf_matrix.shape[0] == 2:
            correct = False
        else:
            correct = True
        chi_2 = ss.chi2_contingency(conf_matrix, correction=correct)[0]
        n = sum(conf_matrix.sum())
        phi2 = chi_2 / n
        r, k = conf_matrix.shape
        phi2corr = max(0, phi2 - ((k - 1) * (r - 1)) / (n - 1))
        r_corr = r - ((r - 1) ** 2) / (n - 1)
        k_corr = k - ((k - 1) ** 2) / (n - 1)
        result = np.sqrt(phi2corr / min((k_corr - 1), (r_corr - 1)))
    return round(result, 6)


def LDA_loocv(data=None, y=None, target=None):
    """
    Function to perform leave one out cross validation using LDA prediction.

    Parameters
    ----------
    data : pandas.core.frame.DataFrame
        Feature matrix
    y : pandas.core.frame.DataFrame
        Target matrix
    target : string
        LDA-specific output target to control on what feature the supervised LDA is fitted, as it would not work for
        continuous features

    Returns
    -------
    acc : float
        float value rounded to 3 decimals, accuracy of the leave one out LDA cross-validation
    bsl : float
        brier score loss
    accs : np.array
        accuracy scores
    bsls : np.array
        brier score loss scores
    fbeta1 : float
        F-beta-1 score
    corr_samples : list
        list of correctly predicted loocv samples
    not_corr_samples : list
        list of not correctly predicted loocv samples
    """
    # define cross validation method
    cv = LeaveOneOut()
    # initialize LDA model
    comp = len(y[target].astype(int).unique()) - 1
    LDA = LinearDiscriminantAnalysis(n_components=comp)
    # Use loocv to evaluate model
    # f beta 1
    scoring = {'f1_score': make_scorer(f1_score, average='weighted')}

    fbeta1 = cross_val_score(estimator=LDA, X=data, y=y[target].astype(int),
                             scoring=scoring['f1_score'], cv=cv, n_jobs=-1)
    # return the accuracy score
    accs = cross_val_score(estimator=LDA, X=data, y=y[target].astype(int),
                             scoring='accuracy', cv=cv, n_jobs=-1)
    acc = np.mean(accs)
    # return brier score loss
    bsls = cross_val_score(estimator=LDA, X=data, y=y[target].astype(int),
                             scoring='neg_brier_score', cv=cv, n_jobs=-1)
    bsl = np.mean(bsls)
    # predictions
    preds = cross_val_predict(LDA, data, y[target].astype(int), cv=cv, n_jobs=-1)
    # compare the preds to the true labels
    corr_preds_idx = (preds == y[target].astype(int))
    corr_samples = y.index[np.where(corr_preds_idx)[0]]
    not_corr_samples = y.index[np.where(~corr_preds_idx)[0]]

    return acc, bsl, accs, bsls, fbeta1, corr_samples, not_corr_samples


def make_mut_table(cell_lines, file):
    # Import the CSV file
    sample_info = pd.read_csv("../data/sample_info.csv").set_index('CCLE_Name')
    genelist = pd.read_excel("../data/circadian_clock_gene_lists_v2.xlsx", "genelist_long")
    genes = genelist['Gene']
    alt_genes = genelist['Alternative gene name']
    dep_map_ids = sample_info[sample_info['cell_line_name'].isin(cell_lines) | sample_info['stripped_cell_line_name'].isin(cell_lines)]['DepMap_ID']
    selected_columns = ['HugoSymbol', 'VariantInfo', 'ModelID']

    ccle_mutations = pd.read_csv(file)
    ccle_mutations = ccle_mutations[selected_columns]

    # Filter the dataframe based on the DepMap_IDs
    filtered_data = ccle_mutations[ccle_mutations['ModelID'].isin(dep_map_ids)]
    filtered_data = filtered_data[filtered_data['HugoSymbol'].isin(genes) | filtered_data['HugoSymbol'].isin(alt_genes)]

    # put back the cell line names
    mapping_series = sample_info.set_index('DepMap_ID')['cell_line_name']
    filtered_data['ModelID'] = filtered_data['ModelID'].map(mapping_series)


    # Pivot the table
    pivoted_data = filtered_data.pivot_table(index='HugoSymbol', columns='ModelID', values='VariantInfo', aggfunc=lambda x: ' + '.join(x.astype(str)))

    id_to_name_map = dict(zip(sample_info['DepMap_ID'], sample_info['stripped_cell_line_name']))
    pivoted_data.columns = [id_to_name_map.get(col, col) for col in pivoted_data.columns]

    return pivoted_data


def get_LDA_metrics(classes, LDs):
    """
    function to retrieve the distance between clusters, within clusters, and the ratio of the two
    in binary one-dimension LDA analysis.
    classes: numpy array with True/False
        example: classes = np.array([True, True, True, False, False, False])
    LDs: numpy array with the transformed data
        example: LDs = np.array([-1.5, -1, -0.5, 0.5, 1, 1.5])
    returns: bcd, wcd, ratio

    """
    classes = np.asarray(classes)
    LDs = np.asarray(LDs)
    if not len(classes) == len(LDs):
        raise ValueError("The 'classes' and 'LDs' arrays must have the same length.")

    # Calculate centroids
    mean_class_true = LDs[classes].mean()
    mean_class_false = LDs[~classes].mean()

    # Between-cluster distance
    bcd = np.abs(mean_class_true - mean_class_false)

    # Within-cluster distance
    distances_true = np.abs(LDs[classes] - mean_class_true)
    distances_false = np.abs(LDs[~classes] - mean_class_false)
    wcd = (distances_true.sum() + distances_false.sum()) / len(LDs)
    ratio = bcd / wcd

    return bcd, wcd, ratio


