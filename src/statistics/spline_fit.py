# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.
from dask import compute
# Import modules
from scipy.interpolate import UnivariateSpline
from sklearn.isotonic import IsotonicRegression
import numpy as np


class SplineFitter:
    def __init__(self, binned_data, config):
        self.binned_data = binned_data
        self.config = config
        self.x, self.y = self._sort_data(binned_data)

    def run(self):
        spline = self._fit_spline()
        spline = self._apply_isotonic_regression(spline)
        return spline

    def _fit_spline(self):
        spline_error = min(self.y) ** 2
        return UnivariateSpline(self.x, self.y, s=spline_error)

    def _apply_isotonic_regression(self, spline):
        spline_predictions = spline(self.x)
        ir = IsotonicRegression(increasing=False)
        fitted_spline = ir.fit_transform(self.x, spline_predictions)
        return fitted_spline

    @staticmethod
    def _sort_data(binned_data):
        sorted_data = binned_data.sort_values("average_distance", ascending=True)
        x = compute(sorted_data["average_distance"])[0]
        y = compute(sorted_data["average_contact_probability"])[0]
        return x, y

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
