from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(self, n_estimators=10, max_depth=None, max_features=None, random_state=None):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit_each_tree(self, args) -> tuple:
        """
        Function for fit each tree in Random Forest
        param: args: tuple (features, labels, number of estimators)
        return: tuple (tree, random features for this tree)
        """
        X, y, n_estimator = args
        n_features = X.shape[1]
        np.random.seed(self.random_state + n_estimator)
        random_features = np.random.choice(np.arange(n_features), self.max_features, replace=False)
        pseudo_index_random = np.random.choice(np.arange(X.shape[0]), X.shape[0], replace=True)
        pseudo_x_samples = X[pseudo_index_random][:, random_features]
        pseudo_y_samples = y[pseudo_index_random]
        dt_model = DecisionTreeClassifier(max_depth=self.max_depth, max_features=self.max_features,
                                          random_state=self.random_state + n_estimator)
        dt_model.fit(pseudo_x_samples, pseudo_y_samples)
        return dt_model, random_features

    def fit(self, X, y, n_jobs: int = 1, parallel_process: bool = True) -> BaseEstimator:
        """Fit Random Forest
        param: X : array-like of shape (n_samples, n_features)
        param: y: array-like of shape (n_samples,)
        param: n_jobs : int, default=1
        param: parallel_process: (Default: True) : parallelization on processes (True) or on threads (False).
        return: fitted RandomForestClassifierCustom object
        """
        self.classes_ = sorted(np.unique(y))
        args_for_fit_tree = zip([X] * self.n_estimators, [y] * self.n_estimators, range(self.n_estimators))

        if parallel_process:
            func_parallel = ProcessPoolExecutor
        else:
            func_parallel = ThreadPoolExecutor

        with func_parallel(n_jobs) as pool:
            results = pool.map(self.fit_each_tree, args_for_fit_tree)

        for dt_model, random_features in results:
            self.trees.append(dt_model)
            self.feat_ids_by_tree.append(random_features)

        return self

    def predict_proba_tree(self, args, n_jobs=1) -> np.array:
        """Predict probabilities for each tree estimator
        :param args: tuple(features, tree, indexs of features for tree)
        :param n_jobs: int, default=1
        :return: np.array of predicted probabilities for tree:
        """
        X, tree, index_features = args
        prediction = tree.predict_proba(X[:, index_features])
        return prediction

    def predict_proba(self, X, n_jobs=1, parallel_process=True) -> np.array:
        """Predict probabilities for RandomForest
        :param: X : array-like of shape (n_samples, n_features)
        :param n_jobs: int, default=1
        :param parallel_process: bool, default=True : parallelization on processes (True) or on threads (False)
        """
        args_for_predict_proba = zip([X] * self.n_estimators, self.trees, self.feat_ids_by_tree)

        if parallel_process:
            func_parallel = ProcessPoolExecutor
        else:
            func_parallel = ThreadPoolExecutor

        with func_parallel(n_jobs) as pool:
            predict_proba_results = pool.map(self.predict_proba_tree, args_for_predict_proba)
        return np.mean(list(predict_proba_results), axis=0)

    def predict(self, X, n_jobs=1, parallel_process=True) -> np.ndarray:
        """Predict labels based on predicted probabilities
        :param: X : array-like of shape (n_samples, n_features)
        :param n_jobs: int, default=1
        :param parallel_process: bool, default=True : parallelization on processes(True) or on threads (False)
        """
        probas = self.predict_proba(X, n_jobs=n_jobs, parallel_process=parallel_process)
        predictions = np.argmax(probas, axis=1)
        return predictions
