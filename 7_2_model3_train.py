# model3_train.py
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import label_binarize
from sklearn.metrics import classification_report, confusion_matrix, roc_curve, roc_auc_score, accuracy_score, precision_score
import joblib

input_path = './model3_fourtypes_region.csv'
prefix = './model3'
model_out = './va_nn_model3.pkl'

df = pd.read_csv(input_path)
X = df.loc[:, df.columns.str.startswith('PC')].values
y = df['region']

X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.2, random_state=42)

model = MLPClassifier(hidden_layer_sizes=(10, 5), activation='relu', solver='adam', max_iter=500, random_state=42)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)
classes = model.classes_

cm = confusion_matrix(y_test, y_pred, labels=classes)
pd.DataFrame(cm, index=[f'True_{c}' for c in classes], columns=[f'Pred_{c}' for c in classes]).to_csv(f'{prefix}_confusion_matrix.csv')

report_dict = classification_report(y_test, y_pred, output_dict=True, zero_division=0)
pd.DataFrame(report_dict).to_csv(f'{prefix}_classification_report.csv')

y_test_bin = label_binarize(y_test, classes=classes)
roc_rows, metrics_rows = [], []
for i, cls in enumerate(classes):
    fpr, tpr, _ = roc_curve(y_test_bin[:, i], y_prob[:, i])
    auc = roc_auc_score(y_test_bin[:, i], y_prob[:, i])
    precision = precision_score(y_test == cls, y_pred == cls)
    accuracy = accuracy_score(y_test, y_pred)
    roc_rows.extend([{'class': cls, 'fpr': float(f), 'tpr': float(t)} for f, t in zip(fpr, tpr)])
    metrics_rows.append({'class': cls, 'auc': float(auc), 'precision': float(precision), 'accuracy': float(accuracy)})

pd.DataFrame(roc_rows).to_csv(f'{prefix}_roc_curve.csv', index=False)
pd.DataFrame(metrics_rows).to_csv(f'{prefix}_roc_metrics.csv', index=False)

joblib.dump(model, model_out)
