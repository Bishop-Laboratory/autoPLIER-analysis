import scipy.spatial
import pandas as pd
import numpy as np
import general_ontology_tools as got
from numpy.linalg import inv
from sklearn.preprocessing import MultiLabelBinarizer

# calculate embedding


def compute_embedding(Y_df, Z_df, mat4):

    Z_df = Z_df
    Z_df = Z_df[Z_df.index.isin(Y_df.columns)]
    Z_mat = Z_df.values

    # sorting columns according to Z
    Y_df = Y_df[Z_df.index]
    mat4 = mat4[Z_df.index]

    # computing normalization values on train dataset
    per_gene_mean = mat4.mean().values
    per_gene_std = mat4.std().values

    y_vectors = (Y_df.values - per_gene_mean) / per_gene_std

    # Using regression to find B
    X = Z_mat

    # linear least squares
    b = inv(X.T.dot(X)).dot(X.T).dot(y_vectors.T)

    return b, Y_df, mat4


def get_pearson_dists(vectors_a, vectors_b):

    vectors_pairs = zip(vectors_a, vectors_b)
    return [scipy.spatial.distance.correlation(a, b) for a, b in vectors_pairs]


def get_pearson_dists_mat(vectors_a, vectors_b):

    return np.array([
            # repeat the A vector for each B vector
            get_pearson_dists(np.array(vector_a[np.newaxis, :]).repeat(len(vectors_b), axis=0), vectors_b)
            for vector_a in vectors_a
    ])


# cell type classification
def classify(B, b, Y_samples, mat4_samples, cell_types_ids, normalize):

    B_mat = B.values
    pearson_dists = get_pearson_dists_mat(B_mat, b.T)

    dists_df = pd.DataFrame(pearson_dists)
    dists_df.columns = Y_samples
    dists_df.index = mat4_samples

    sample2types = {
        sample: list(map(got.get_term_name, types_ids))
        for sample, types_ids in cell_types_ids.items()
    }

    b_T = pd.DataFrame(b.T)
    mlb = MultiLabelBinarizer()
    types_per_b_samples = B.index.map(sample2types).values

    train_dummies = pd.DataFrame(mlb.fit_transform(types_per_b_samples), columns=mlb.classes_, index=B.index)
    types_dummies = pd.concat([train_dummies, b_T])[train_dummies.columns]

    test_sample2closest_train = {
        # get the column for current test sample and find the index of raw the closest value
        sample_id: dists_df[sample_id].sort_values().index[:10]
        for sample_id in Y_samples
    }

    cell_types_y_predicted = pd.DataFrame([
        types_dummies.loc[test_sample2closest_train[test_id]].sum()
        for test_id in Y_samples
    ])

    if normalize:
        # idf normalization

        count_per_ctype = train_dummies.sum()
        idf = len(train_dummies) / count_per_ctype
        idf[:] = np.log(idf.values)

        labels = cell_types_y_predicted.sum() / cell_types_y_predicted.sum().sum()

        # normalize values with TF-IDF

        labels = labels * idf[labels.index]
        labels = labels[labels.iloc[:] > 0]
        labels = labels.sort_values(ascending=False)
        return labels

    else:

        # return sorted cell types

        sum_per_ctype = cell_types_y_predicted.sum()
        labels = sum_per_ctype[(sum_per_ctype > 1)].sort_values(ascending=False)
        return labels


# get largest LVs per Celltype
def getTopLVs(cell_types_ids, cell_label, B, n_LVs):

    sample2types = {
        sample: list(map(got.get_term_name, types_ids))
        for sample, types_ids in cell_types_ids.items()
    }
    types_per_b_samples = B.index.map(sample2types).values

    mlb = MultiLabelBinarizer()
    train_dummies = pd.DataFrame(mlb.fit_transform(types_per_b_samples), columns=mlb.classes_, index=B.index)
    samples = train_dummies[cell_label].index
    b = B.transpose()
    LVcount = pd.DataFrame(columns=("LV", "count"))

    for sample in samples:
        topLVs = b[sample].sort_values(ascending=False)[:n_LVs]
        for LV in topLVs.index:
            if "," in LV:
                LV = LV.split(",")[1]
            if LV not in LVcount["LV"].unique():
                newrow = {"LV": LV, "count": 1}
                LVcount = LVcount.append(newrow, ignore_index=True)
            else:
                LVcount.loc[LVcount['LV'] == LV, 'count'] += 1
    LVcount["Fraction"] = LVcount["count"]/len(samples)
    return LVcount


# compare to cell type
def cellcompare(cell_types_ids, cell_label, B, n_LVs, b):

    Result = pd.DataFrame(columns=("LV", "Score"))
    LVcountsamples = getTopLVs(cell_types_ids, cell_label, B, n_LVs)
    b_T = pd.DataFrame(b.T, columns=B.columns)
    b_TT = b_T.transpose()
    for sample in b_TT:

        topLvs = b_TT[sample].sort_values(ascending=False)[:n_LVs]

        for topLv in topLvs.index:
            if "," in topLv:
                topLv = topLv.split(",")[1]
            try:
                score = LVcountsamples.loc[LVcountsamples['LV'] == topLv]["Fraction"].item()
            except:
                score = -1

            if topLv not in Result["LV"].unique():
                newrow = {"LV": topLv, "Score": score}
                Result = Result.append(newrow, ignore_index=True)
    Result = Result.sort_values(by=['Score'], ascending=False).reset_index()

    Result_finalized = Result["Score"]
    Result_finalized.index = Result["LV"]
    return Result_finalized
