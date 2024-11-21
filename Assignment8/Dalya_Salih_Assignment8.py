import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.metrics import roc_auc_score

# Load the datasets
train_file_path = 'train.csv'
test_file_path = 'test.csv'

train_data = pd.read_csv(train_file_path)
test_data = pd.read_csv(test_file_path)

# Separate features and target variable
X = train_data.drop(columns=['breast_cancer', 'id'])
y = train_data['breast_cancer']

# Identify categorical and numerical columns
categorical_cols = ['predicted_ancestry']
numerical_cols = [col for col in X.columns if col not in categorical_cols]

# Preprocess the data
preprocessor = ColumnTransformer(
    transformers=[
        ('num', StandardScaler(), numerical_cols),
        ('cat', OneHotEncoder(), categorical_cols)
    ]
)

X_preprocessed = preprocessor.fit_transform(X)

# Reduce dimensionality using Variance Threshold
selector = VarianceThreshold(threshold=0.05)  # Stricter feature reduction
X_reduced = selector.fit_transform(X_preprocessed)

# Split into train and validation sets
X_train, X_val, y_train, y_val = train_test_split(X_reduced, y, test_size=0.2, random_state=42)

# Take a subset of the training data for faster experimentation
X_train_sample, _, y_train_sample, _ = train_test_split(
    X_train, y_train, train_size=0.5, random_state=42
)

# Define logistic regression model
log_reg = LogisticRegression(penalty='elasticnet', solver='liblinear', max_iter=1000, random_state=42)

# Set up GridSearchCV with a reduced parameter grid
param_grid = {
    'C': [0.1, 1],           # Smaller range
    'l1_ratio': [0.5]        # Single ElasticNet mixing parameter
}
grid_search = GridSearchCV(
    estimator=log_reg,
    param_grid=param_grid,
    scoring='roc_auc',  # Optimize AUC-ROC
    cv=3,               # Reduced cross-validation folds
    verbose=1,
    n_jobs=-1           # Parallelize across all CPU cores
)

# Run grid search
grid_search.fit(X_train_sample, y_train_sample)

# Get the best parameters and model
best_log_reg = grid_search.best_estimator_
print(f"Best Parameters: {grid_search.best_params_}")

# Evaluate on validation set
y_val_pred_prob = best_log_reg.predict_proba(X_val)[:, 1]
auc = roc_auc_score(y_val, y_val_pred_prob)
print(f"Optimized Logistic Regression Validation AUC: {auc:.4f}")

# Preprocess the test set and predict probabilities
X_test = test_data.drop(columns=['id'])
X_test_preprocessed = preprocessor.transform(X_test)
X_test_reduced = selector.transform(X_test_preprocessed)  # Apply the same selector
test_pred_prob = best_log_reg.predict_proba(X_test_reduced)[:, 1]

# Prepare submission file
submission = pd.DataFrame({'id': test_data['id'], 'breast_cancer': test_pred_prob})
submission.to_csv('submission.csv', index=False)

print("Submission file 'submission.csv' has been saved.")
