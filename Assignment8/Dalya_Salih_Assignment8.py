import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer

# Load the datasets
train_file_path = 'train.csv'
test_file_path = 'test.csv'

# Load data
train_data = pd.read_csv(train_file_path)
test_data = pd.read_csv(test_file_path)

# Separate features and target variable from training data
X = train_data.drop(columns=['breast_cancer', 'id'])  # Drop target and ID
y = train_data['breast_cancer']                      # Target variable

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

# Split training data into train and validation sets
X_train, X_val, y_train, y_val = train_test_split(X_preprocessed, y, test_size=0.2, random_state=42)

# Train logistic regression as a baseline model
log_reg = LogisticRegression(max_iter=1000, random_state=42)
log_reg.fit(X_train, y_train)

# Predict probabilities on validation set
y_val_pred_prob = log_reg.predict_proba(X_val)[:, 1]

# Calculate AUC-ROC for validation
auc = roc_auc_score(y_val, y_val_pred_prob)
print(f"Validation AUC: {auc:.4f}")

# Preprocess the test set
X_test = test_data.drop(columns=['id'])  # Drop ID column
X_test_preprocessed = preprocessor.transform(X_test)

# Predict probabilities on the test set
test_pred_prob = log_reg.predict_proba(X_test_preprocessed)[:, 1]

# Prepare submission file
submission = pd.DataFrame({'id': test_data['id'], 'breast_cancer': test_pred_prob})
submission.to_csv('submission.csv', index=False)

print("Submission file 'submission.csv' has been saved.")
