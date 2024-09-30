import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler,MinMaxScaler,RobustScaler,LabelEncoder
from sklearn.metrics import accuracy_score,auc,roc_curve,precision_score,recall_score,f1_score,roc_auc_score
# from sklearn.linear_model import LogisticRegression
import hyperopt
from hyperopt import hp, fmin, tpe, Trials, partial
from hyperopt.early_stop import no_progress_loss
from sklearn.model_selection import KFold,cross_validate,LeaveOneOut,GridSearchCV,train_test_split
def load_data(path='F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/subdata/Intersect_2p/subdata_fpkm.csv'
              # ,path2='F:/my projects/02.论文/论文四_桥本甲状腺炎/Project/Datas/Datas/subdata/Intersect_2p/subdata_fpkm_derived.csv'
              ):
    data = pd.read_csv(path)
    # data2 = pd.read_csv(path2)
    target = data.Group
    lbe = LabelEncoder()
    target = lbe.fit_transform(target)
    feature = data.drop(columns='Group')
    # feature2 = data2
    X_train,X_test,y_train,y_test = train_test_split(feature,target,test_size=20,stratify=target,random_state=42)
    # X_train2,X_test2,y_train2,y_test2 = train_test_split(feature2,target,test_size=20,stratify=target,random_state=42)
    scaler = RobustScaler()
    X_train_scale = scaler.fit_transform(X_train)
    X_test_scale = scaler.transform(X_test)
    # X_train2_scale = scaler.fit_transform(X_train2)
    # X_test2_scale = scaler.transform(X_test2)
    return (X_train_scale,y_train
            ,X_test_scale,y_test
            )


# def param_hyperopt(max_evals=3000):
#     # 保存迭代过程
#     trials = Trials()
#     # 设置提前停止
#     early_stop_fn = no_progress_loss(500)
#     # 定义代理模型
#     params_best = fmin(hyperopt_objective
#                        , space=param_grid_simple
#                        , algo=tpe.suggest
#                        , max_evals=max_evals
#                        , verbose=True
#                        , trials=trials
#                        , early_stop_fn=early_stop_fn
#                        )
#     # 打印最优参数，fmin会自动打印最佳分数
#     print("\n", "\n", "best params: ", params_best,
#           "\n")
#     return params_best, trials


def cv_res(X_train, y_train, clf):
    cv = KFold(n_splits=k, shuffle=True, random_state=random_state)
    res = cross_validate(clf
                         , X=X_train
                         , y=y_train
                         , cv=cv
                         , n_jobs=16
                         , scoring='accuracy')

    return np.mean(res['test_score'])


def test_res(X_train, y_train, X_test, y_test, clf):
    clf.fit(X_train, y_train)
    train_acc = clf.score(X_train, y_train)
    train_auc = roc_auc_score(y_train, clf.predict_proba(X_train)[:, 1])
    test_acc = clf.score(X_test, y_test)
    test_auc = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    return train_acc, train_auc, test_acc, test_auc

def clf_res(clf,X_train,y_train,X_test,y_test,k,random_state):
    cv = KFold(n_splits=k,shuffle=True,random_state=random_state)
    res = cross_validate(clf
                        ,X=X_train
                        ,y=y_train
                        ,cv=cv
                        ,scoring='accuracy')
    cv_res = np.mean(res['test_score'])
    print('cv score: ',cv_res)
    clf.fit(X_train,y_train)
    clf_train_acc = clf.score(X_train, y_train)
    clf_train_auc = roc_auc_score(y_train, clf.predict_proba(X_train)[:, 1])
    clf_test_acc = clf.score(X_test,y_test)
    clf_test_auc = roc_auc_score(y_test,clf.predict_proba(X_test)[:,1])
    print('train acc: ', clf_train_acc)
    print('train auc: ', clf_train_auc)
    print('test acc: ',clf_test_acc)
    print('test auc: ',clf_test_auc)