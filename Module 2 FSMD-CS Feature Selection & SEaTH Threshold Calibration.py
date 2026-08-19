import numpy as np
import pandas as pd
from scipy.optimize import fsolve


# 1. Mahalanobis Distance & FSMD-CS Selection
def mahalanobis_distance(x, y):
    """Calculates 1D Mahalanobis Distance between two sample distributions."""
    mean_x, mean_y = np.mean(x), np.mean(y)
    var_x, var_y = np.var(x, ddof=1), np.var(y, ddof=1)
    pooled_var = (var_x + var_y) / 2.0
    if pooled_var == 0:
        return 0.0
    return np.abs(mean_x - mean_y) / np.sqrt(pooled_var)


def fsmd_cs_feature_selection(df, class_col, features):
    """Executes FSMD-CS greedy search to determine feature hierarchy."""
    classes = list(df[class_col].unique())
    selected_features = {}
    remaining_classes = classes.copy()

    while len(remaining_classes) > 1:
        best_j_avg = -1.0
        best_combo = None  # Tuple: (target_class, best_feature)

        for c in remaining_classes:
            other_classes = [k for k in classes if k != c]
            for f in features:
                distances = []
                c_data = df[df[class_col] == c][f]
                for o in other_classes:
                    o_data = df[df[class_col] == o][f]
                    distances.append(mahalanobis_distance(c_data, o_data))

                j_avg = np.mean(distances)
                if j_avg > best_j_avg:
                    best_j_avg = j_avg
                    best_combo = (c, f)

        target_class, best_feature = best_combo
        selected_features[target_class] = {
            'feature': best_feature,
            'j_avg': best_j_avg
        }
        # Exclude the successfully separated class from the candidate pool
        remaining_classes.remove(target_class)

    # Assign residual status to the final class
    selected_features[remaining_classes[0]] = {'feature': 'Residual', 'j_avg': 0.0}
    return selected_features


# 2. SEaTH Threshold Optimization & Calibration
def calculate_seath_threshold(m1, s1, m2, s2):
    """Solves for the intersection point of two Gaussian PDFs (SEaTH theoretical threshold)."""

    def pdf_diff(x):
        pdf1 = (1.0 / (s1 * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - m1) / s1) ** 2)
        pdf2 = (1.0 / (s2 * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - m2) / s2) ** 2)
        return pdf1 - pdf2

    init_guess = (m1 + m2) / 2.0
    t_star = fsolve(pdf_diff, init_guess)[0]
    return t_star


def get_calibrated_threshold(df, target_class, feature, class_col, j_score):
    """Calibrates SEaTH threshold based on class separability (J-score)."""
    target_data = df[df[class_col] == target_class][feature]
    other_data = df[df[class_col] != target_class][feature]

    m1, s1 = target_data.mean(), target_data.std()
    m2, s2 = other_data.mean(), other_data.std()

    t_star = calculate_seath_threshold(m1, s1, m2, s2)

    # 3-Tier calibration rule
    if j_score >= 1.9:
        t_calibrated = t_star
    elif 1.0 <= j_score < 1.9:
        t_calibrated = (t_star + m1) / 2.0
    else:
        t_calibrated = m1

    direction = 'lte' if m1 < m2 else 'gte'
    return t_calibrated, direction


# 3. Execution Pipeline
if __name__ == "__main__":
    # Load exported CSV dataset from GEE
    df = pd.read_csv("MPS_Sample_Features.csv")
    feature_cols = [c for c in df.columns if c not in ['class_code', 'system:index', '.geo']]

    # Step 1: Run FSMD-CS feature selection
    selection_results = fsmd_cs_feature_selection(df, 'class_code', feature_cols)

    # Step 2: Compute calibrated thresholds
    rules = []
    for cls, info in selection_results.items():
        if info['feature'] == 'Residual':
            rules.append({'class': cls, 'type': 'residual'})
            print(f"Class {cls}: Assigned as Residual (Final Layer)")
            continue

        thresh, direction = get_calibrated_threshold(
            df, cls, info['feature'], 'class_code', info['j_avg']
        )
        rules.append({
            'class': cls,
            'feature': info['feature'],
            'threshold': round(thresh, 4),
            'direction': direction
        })
        print(f"Class {cls}: Feature = {info['feature']}, Threshold = {thresh:.4f} ({direction})")
