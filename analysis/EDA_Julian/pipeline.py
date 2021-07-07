import cello_multiplier as cm
import pandas as pd
import json
from sklearn.preprocessing import MultiLabelBinarizer
from onto_lib import general_ontology_tools as got
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.metrics import precision_score, recall_score
import numpy as np
import warnings
from tqdm.auto import tqdm
from itertools import repeat
from onto_lib.load_ontology import load

# fscore metric
def fscore(p, r):
    denom = p + r or 1

    return 2 * (p * r) / denom


# create target list for given cell type with 1 being that cell type and 0 being any other cell type
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

    train_Y_transformed = train_Y_df
    test_Y_transformed = test_Y_df

    return target_train.values, target_test.values, train_Y_transformed, test_Y_transformed


def LR_classify(solver, penalty, tol, C, fit_intercept, intercept_scaling, max_iter, class_weight, all_types,
                type2samples, train_Y_df, test_Y_df, oversample):
    warnings.filterwarnings('ignore')  # gets rid of sklearn convergence warning

    successful_celltypes = []
    unsuccessful_celltypes = []
    test_prs = []
    test_precision = []
    test_recall = []
    fscores = []

    for cell_type in tqdm(all_types):
        train_target, train_test, train_data, test_data = set_target(
            cell_type, type2samples, train_Y_df, test_Y_df, oversample= oversample,
            max_neg_pos_ratio=2
        )

        if (1 in train_target and 1 in train_test and 0 in train_target and 0 in train_test):
            # lasso penalty
            clf = LogisticRegression(solver=solver, penalty=penalty, random_state=111, tol=tol,
                                     C=C, fit_intercept=fit_intercept, intercept_scaling=intercept_scaling,
                                     max_iter=max_iter, class_weight=class_weight)

            clf.fit(train_data, train_target)
            target_pred = clf.predict(test_data)
            test_pr = metrics.average_precision_score(train_test, target_pred)
            test_precision += [precision_score(train_test, target_pred)]
            test_recall += [recall_score(train_test, target_pred)]
            fscores += [fscore(test_precision[-1], test_recall[-1])]
            successful_celltypes += [cell_type]
            # plot precision recall curve for celltype
            # if(cell_type in curves):
            # disp = plot_precision_recall_curve(clf, test_data, train_test)
            # disp.ax_.set_title('2-class Precision-Recall curve '+str(cell_type)+': AP={0:0.2f}'.format(test_pr))
        else:
            unsuccessful_celltypes += [cell_type]

    p = np.mean(test_precision)
    r = np.mean(test_recall)
    f = np.mean(fscores)
    f_micro = fscore(p, r)

    report = pd.DataFrame(list(zip(successful_celltypes, test_precision, test_recall, fscores)),
                          columns=["celltype", "precision score", "recall score", "f score"])
    # print(f'precision: {p:.4f}, recall: {r:.4f}, f1: {f:.4f}, f1 micro avg: {f_micro:.4f}')#

    return (p, r, f, f_micro)


def classification_pipeline(B_df, split_dir, solver = "liblinear", penalty = "l1", tol = 0.0001, C = 1.0,
                            fit_intercept = True, intercept_scaling = 1.0, class_weight = None, max_iter = 100,
                            oversample = True):
    # load CellO Data
    B_df1, Z_df, labels, per_gene_mean, per_gene_std, train_dummies_specific, train_dummies_full, classifiers = \
        cm.get_default_mats()


    # split
    with open(split_dir / 'validation_bulk_experiments.json', 'r') as f:
        validation_egs = json.load(f)

    with open(split_dir/ 'pre_training_bulk_experiments.json', 'r') as f:
        train_egs = json.load(f)

    train_Y_df = B_df[B_df.index.isin(train_egs)]
    test_Y_df = B_df[B_df.index.isin(validation_egs)]

    ont_name_to_ont_id = {"EFO_CL_DOID_UBERON_CVCL": "17"}
    ont_id_to_og = {
        x: load(x)[0]
        for x in list(ont_name_to_ont_id.values())
    }

    # create list of samples by celltypes
    sample2types = {
        sample: list(map(got.get_term_name, types_ids, repeat(ont_id_to_og)))
        for sample, types_ids in labels.items()
    }

    types_per_b_samples = B_df.index.map(sample2types).values

    mlb = MultiLabelBinarizer()
    types_per_b_samples = B_df.index.map(sample2types).values

    samples_dummies = pd.DataFrame(mlb.fit_transform(types_per_b_samples), columns=mlb.classes_, index=B_df.index)

    celltypes = samples_dummies.columns.tolist()


    type2samples = {
        type_: samples_dummies.index[samples_dummies[type_] == 1].tolist()
        for type_ in celltypes
    }

    types_sizes = samples_dummies.sum()
    types_with_data = types_sizes.index.values
    all_types = samples_dummies.columns

    p, r, f, f_micro = LR_classify(solver, penalty, tol, C, fit_intercept, intercept_scaling, max_iter, class_weight,
                                   all_types, type2samples, train_Y_df, test_Y_df, oversample)

    return p, r, f, f_micro



