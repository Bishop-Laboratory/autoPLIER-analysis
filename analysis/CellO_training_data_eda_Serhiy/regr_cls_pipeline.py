import numpy as np
from tqdm.auto import tqdm
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MultiLabelBinarizer
from onto_lib import general_ontology_tools as got
from pathlib import Path
import pandas as pd
import pickle
import json
from numpy.linalg import inv
from sklearn.metrics import precision_score, recall_score


def pickle_load(filename):
    return pickle.load(open(str(filename), 'rb'))


def pickle_dump(filename, data):
    return pickle.dump(data, open(str(filename), 'wb'))


data_dir = Path('../../data')
cello_dir = data_dir / 'CellO_data/bulk_RNA_seq_training_set'
split_dir = cello_dir / 'pretraining_validation_split'


mat4_assets_path = data_dir / 'mat4_assets.pkl'

if not mat4_assets_path.is_file():
    print(f'{mat4_assets_path} is missing. reading from full mat4')
    # transposing, so we have columns as features and rows as sample vectors
    full_Y_df = pd.read_csv(data_dir / 'mat4.csv').transpose()

    with open(split_dir / 'validation_bulk_experiments.json', 'r') as f:
        validation_egs = json.load(f)

    test_Y_df = full_Y_df[full_Y_df.index.isin(validation_egs)]
    train_Y_df = full_Y_df[~full_Y_df.index.isin(validation_egs)]

    per_gene_mean = train_Y_df.mean()
    per_gene_std = train_Y_df.std()

    pickle_dump(mat4_assets_path, (test_Y_df, per_gene_mean, per_gene_std))

test_Y_df, per_gene_mean, per_gene_std = pickle_load(mat4_assets_path)


def fscore(p, r):
    denom = p + r or 1

    return 2 * (p * r) / denom


def LR_classify(solver, penalty, tol, C, fit_intercept, intercept_scaling, max_iter, class_weight, all_types,
                type2samples, train_Y_df, test_Y_df, oversample):
    #     warnings.filterwarnings('ignore')  # gets rid of sklearn convergence warning

    successful_celltypes = []
    unsuccessful_celltypes = []
    test_precision = []
    test_recall = []
    fscores = []

    for cell_type in tqdm(all_types):
        train_target, train_test, train_data, test_data = set_target(
            cell_type, type2samples, train_Y_df, test_Y_df, oversample=oversample,
            max_neg_pos_ratio=2
        )

        if (1 in train_target and 1 in train_test and 0 in train_target and 0 in train_test):
            # lasso penalty
            clf = LogisticRegression(solver=solver, penalty=penalty, random_state=111, tol=tol,
                                     C=C, fit_intercept=fit_intercept, intercept_scaling=intercept_scaling,
                                     max_iter=max_iter, class_weight=class_weight)

            clf.fit(train_data, train_target)
            target_pred = clf.predict(test_data)
            test_precision += [precision_score(train_test, target_pred)]
            test_recall += [recall_score(train_test, target_pred)]
            fscores += [fscore(test_precision[-1], test_recall[-1])]
            successful_celltypes += [cell_type]
        else:
            unsuccessful_celltypes += [cell_type]

    p = np.mean(test_precision)
    r = np.mean(test_recall)
    f = np.mean(fscores)

    return (p, r, f)


def set_target(celltype, type2samples, train_Y_df, test_Y_df, oversample, max_neg_pos_ratio):
    samplelist = type2samples[celltype]
    samplelist_train = [x for x in samplelist if x in train_Y_df.index.values]
    samplelist_test = [x for x in samplelist if x in test_Y_df.index.values]

    if oversample == True and len(samplelist_train) > 0 and len(samplelist_test) > 0:
        len_negative_train = len(train_Y_df) - len(samplelist_train)
        len_positive_train = len(samplelist_train)

        neg_pos_ratio = len_negative_train / len_positive_train

        # limit the max oversampling ratio
        neg_pos_ratio = min(neg_pos_ratio, max_neg_pos_ratio)

        if neg_pos_ratio > 1:
            df_train = train_Y_df.loc[samplelist_train].sample(
                n=int((neg_pos_ratio - 1) * len_positive_train),
                replace=True, random_state=111
            )
            train_Y_df = train_Y_df.append(df_train)

    target_train = pd.Series(0, index=train_Y_df.index)
    target_train.loc[samplelist_train] = 1

    target_test = pd.Series(0, index=test_Y_df.index)
    target_test.loc[samplelist_test] = 1

    return target_train.values, target_test.values, train_Y_df, test_Y_df


def get_all_types():
    with open(cello_dir / 'bulk_labels.json', 'r') as f:
        sample2cell_types_ids = json.load(f)

    # create list of samples by celltypes
    sample2types = {
        sample: list(map(got.get_term_name, types_ids))
        for sample, types_ids in sample2cell_types_ids.items()
    }

    all_samples = pd.Series(sample2cell_types_ids.keys())
    types_per_samples = all_samples.map(sample2types).values

    mlb = MultiLabelBinarizer()
    samples_dummies = pd.DataFrame(mlb.fit_transform(types_per_samples), columns=mlb.classes_, index=all_samples)

    all_types = samples_dummies.columns.tolist()

    type2samples = {
        type_: samples_dummies.index[samples_dummies[type_] == 1].tolist()
        for type_ in all_types
    }

    return all_types, type2samples


def get_classification_scores(
        B_train_df, b_test_df, solver='liblinear', penalty='l1', tol=0.0001, C=1.0,
        fit_intercept=True, intercept_scaling=1.0, class_weight=None, max_iter=100,
        oversample=True):

    return LR_classify(solver, penalty, tol, C, fit_intercept, intercept_scaling, max_iter, class_weight,
                       all_types, type2samples, B_train_df, b_test_df, oversample)


all_types, type2samples = get_all_types()

plier_train_folders = [
    'plierResult-cello_train',
    'plierResult-cello_train-256_dims'
]

# regression and classification for each PLIER train
for plier_assets in plier_train_folders:
    Z_df = pd.read_csv(data_dir / plier_assets / 'Z.csv')
    Z_mat = Z_df.values
    Z_genes = Z_df.index

    test_Y_df = test_Y_df[Z_genes]
    print(f'test data shape {test_Y_df.shape}')

    # transposing, so we have columns as features and rows as sample vectors
    B_train_df = pd.read_csv(data_dir / plier_assets / 'B.csv').transpose()
    print(f'train embeddings data shape {B_train_df.shape}')

    # making PLIER train set is not overlapping with our test set
    B_samples = B_train_df.index.values
    test_samples = test_Y_df.index.values
    assert not set(B_samples) & set(test_samples), 'train set is overlapping with test set'

    test_y_vectors = (test_Y_df.values - per_gene_mean[Z_genes].values) / per_gene_std[Z_genes].values

    # linear least squares: https://machinelearningmastery.com/solve-linear-regression-using-linear-algebra/
    b_test = inv(Z_mat.T.dot(Z_mat)).dot(Z_mat.T).dot(test_y_vectors.T)

    b_test_df = pd.DataFrame(b_test.T, index=test_Y_df.index)

    results = get_classification_scores(B_train_df, b_test_df)

    # TODO: log results

    print(results)
