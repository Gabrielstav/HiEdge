# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from scipy.interpolate import PchipInterpolator, UnivariateSpline
from sklearn.isotonic import IsotonicRegression
import numpy as np


class SplineDataPreparer:

    def __init__(self, binned_data):
        self.binned_data = binned_data

    def prepare_data_for_spline(self):
        x = self.binned_data["average_distance"]
        y = self.binned_data["average_contact_probability"]
        yerr = self.calculate_standard_error()

        return x, y, yerr

    def calculate_standard_error(self):
        standard_error = (self.binned_data["interaction_std"]
                          / np.sqrt(self.binned_data["interaction_count"]))
        return standard_error


class SplineFitter:
    def __init__(self, x, y, yerr, config):
        self.x = x
        self.y = y
        self.yerr = yerr
        self.config = config

    def fit_spline(self):
        if not np.all(np.diff(self.x) < 0):
            raise ValueError("x values must be in monotonically decreasing order")
        spline_error = min(self.y) ** 2
        return UnivariateSpline(self.x, self.y, s=spline_error)

    def fit_with_isotonic_regression(self, spline):
        original_spline_y = spline(self.x)
        ir = IsotonicRegression(increasing=False)
        corrected_spline_y = ir.fit_transform(self.x, original_spline_y)
        return PchipInterpolator(self.x, corrected_spline_y)

    def multiple_passes_fit(self):
        spline = self.fit_spline()
        for _ in range(self.config.pipeline_settings.spline_passes - 1):
            spline = self.fit_with_isotonic_regression(spline)
        return spline

class SplineAnalysis:
    def __init__(self, x, y, spline):
        self.x = x
        self.y = y
        self.spline = spline

    def calculate_residuals(self):
        spline_y = self.spline(self.x)
        return self.y - spline_y

    def calculate_mse(self):
        residuals = self.calculate_residuals()
        return np.mean(residuals ** 2)

    def calculate_r_squared(self):
        total_variance = np.var(self.y)
        residuals = self.calculate_residuals()
        residual_variance = np.var(residuals)
        return 1 - (residual_variance / total_variance)
