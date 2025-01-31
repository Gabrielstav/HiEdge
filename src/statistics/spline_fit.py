# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from scipy.interpolate import UnivariateSpline
from sklearn.isotonic import IsotonicRegression
import dask.dataframe as dd
import numpy as np

class SplineFitter:

    def __init__(self, data):
        self.data = data
        self.spline_error = None
        self.x = None
        self.y = None
        self.mono_y = None
        self._initial_spline = None
        self.mono_spline = None

    def fit_spline(self):
        self._sort_input_data()
        self._get_spline_error()
        self._fit_initial_spline()
        self._enforce_monotonicity()
        return self

    def predict(self, x_new) -> float:
        return self.mono_spline(x_new)

    def _fit_initial_spline(self):
        initial_spline = UnivariateSpline(self.x, self.y, s=self.spline_error)
        print(f"INITIAL X: {self.x}")
        print(f"INITIAL Y: {self.y}")
        self._initial_spline = initial_spline

    # TODO: Fix NaN when cycles are enabled in config (input is 0 for self-interacting bins)
    def _enforce_monotonicity(self):
        y_mono = IsotonicRegression(increasing=False, y_min=min(self.y), y_max=max(self.y)).fit_transform(self.x, self._initial_spline(self.x))
        print(f"MONO Y: {y_mono}")
        adjusted_spline = UnivariateSpline(self.x, y_mono, s=self.spline_error)
        self.mono_spline = adjusted_spline

    def _sort_input_data(self):
        if isinstance(self.data, dd.DataFrame):
            sorted_data = self.data.compute().sort_values("average_genomic_distance")
        else:
            sorted_data = self.data.sort_values("average_genomic_distance")

        self.x = sorted_data["average_genomic_distance"].values
        self.y = sorted_data["average_contact_probability"].values

    def _get_spline_error(self):
        spline_error = min(self.y) ** 2
        self.spline_error = spline_error

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
