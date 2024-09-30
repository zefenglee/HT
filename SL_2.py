from sklearn.preprocessing import StandardScaler,MinMaxScaler,RobustScaler,LabelEncoder
from sklearn.metrics import accuracy_score,auc,roc_curve,precision_score,recall_score,f1_score
import warnings
warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
import numpy as np
import pandas as pd
#LR,DTC,RF,GBDT,XGB,CAT,ETC,LGB,SVC,KNN,MLP,NB
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import ExtraTreesClassifier ,StackingClassifier
from sklearn.svm import SVC
from catboost import CatBoostClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import StackingClassifier
from xgboost.sklearn import XGBClassifier
from sklearn.naive_bayes import GaussianNB
import lightgbm as lgb
import random
import hyperopt
from hyperopt import hp, fmin, tpe, Trials, partial
from hyperopt.early_stop import no_progress_loss
from sklearn.model_selection import KFold,cross_validate,LeaveOneOut,GridSearchCV,train_test_split,RandomizedSearchCV
import time
from sklearn.preprocessing  import StandardScaler
from sklearn.metrics import roc_auc_score
import joblib
from joblib import Parallel,delayed
import sys
sys.path.append('F:\\my projects\\02.论文\\论文四_桥本甲状腺炎\\Project\\Codes\\202312_Hashimotos_thyroiditis\\05 modeling')
from Modeling_functions import *
import os

class SL_2():
    def __init__(self):
        self.estimators = self.model_init()

    def model_init(self):
        path = 'F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/07.Modeling/Layer1/'
        estimators = [joblib.load(path+i) for i in os.listdir('F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Results/07.Modeling/Layer1')]
        return estimators

    def lay1(self, X, y,random_state = 42):
        cv = KFold(n_splits=10, shuffle=True, random_state=random_state)
        index_split = cv.split(X, y)
        lay1_res_list = []
        idx_valid = []
        for estimator, (idx_tra, idx_val) in zip(self.estimators, index_split):
            estimator.fit(X[idx_tra, :], y[idx_tra])
            lay1_res_list.append(estimator.predict_proba(X[idx_val, :])[:, 1])
            idx_valid.extend(idx_val.tolist())
        lay1_output = np.concatenate(lay1_res_list, axis=0)
        lay1_output = np.array(pd.Series(lay1_output, index=idx_valid).sort_index()).reshape(-1, 1)
        return lay1_output

    def lay2(self, X, y):
        lay1_output = self.lay1(X, y)
        lay2_input = np.concatenate([lay1_output, X], axis=1)
        return lay2_input

    def fit(self, X, y, random_state=42):
        X = self.lay2(X, y)
        final_estimator = LogisticRegression(random_state=random_state,C=0.005,class_weight='balanced')
        base_models = self.model_init()
        base_model = [('lr', base_models[0]), ('rf', base_models[1]), ('gbdt', base_models[2]), ('etc', base_models[3]),
                      ('svc', base_models[4])
            , ('xgb', base_models[5]), ('mlp', base_models[6]), ('knn', base_models[7]), ('nb', base_models[8]),
                      ('cat', base_models[9])]
        self.clf_stack = StackingClassifier(estimators=base_model
                                            , final_estimator=final_estimator
                                            , stack_method="auto"
                                            )
        self.clf_stack.fit(X, y)

    def predict_proba(self, X):
        lay1_prob = np.mean([estimator.predict_proba(X)[:, 1] for estimator in self.estimators], 0).reshape(-1, 1)
        X = np.concatenate([lay1_prob, X], axis=1)
        return self.clf_stack.predict_proba(X)

    def predict(self, X):
        prob = self.predict_proba(X)
        return np.argmax(prob, axis=1)

    def score(self, X, y):
        pred = self.predict(X)
        return np.mean(pred == y)