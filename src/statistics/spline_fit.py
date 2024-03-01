# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from scipy.interpolate import UnivariateSpline
from sklearn.isotonic import IsotonicRegression
import dask.dataframe as dd
from dask import compute
import numpy as np
from collections import namedtuple

SortedData = namedtuple("SortedData", ["x", "y"])

class SplineFitter:

    def __init__(self, binned_data: dd.DataFrame, config):
        sorted_data = self._get_sorted_data(binned_data)
        self.binned_data = binned_data
        self.config = config
        self.x = sorted_data.x
        self.y = sorted_data.y
        self.spline_error = min(self.y) ** 2

    def run(self) -> UnivariateSpline:
        initial_spline = self._fit_initial_spline()
        adjusted_y = self._apply_isotonic_regression(initial_spline)
        final_spline = self._fit_final_spline(adjusted_y)
        return final_spline

    def _fit_initial_spline(self) -> UnivariateSpline:
        # check x and y for null values
        if np.isnan(self.x).any() or np.isnan(self.y).any():
            raise ValueError("X or Y contains null values")

        print(f"Y values before spline fit: {self.y}")
        print(f"X values before spline fit: {self.x}")
        return UnivariateSpline(self.x, self.y, s=self.spline_error)

    def _apply_isotonic_regression(self, initial_spline) -> UnivariateSpline:
        # print y vals from initial spline fit
        print(f"Initial spline y: {initial_spline(self.y)}")
        y_predictions = initial_spline(self.x)
        ir = IsotonicRegression(increasing=False)
        print(f"Y predictions: {y_predictions}")
        adjusted_y = ir.fit_transform(self.x, y_predictions)
        return adjusted_y

    def _fit_final_spline(self, adjusted_y) -> UnivariateSpline:
        # Fit a new spline to the isotonic regression adjusted y-values
        return UnivariateSpline(self.x, adjusted_y, s=self.spline_error)

    @staticmethod
    def _get_sorted_data(binned_data: dd.DataFrame) -> SortedData:
        """Sort data and return a namedtuple containing x (distances) and y (probabilities)."""
        sorted_data = binned_data.sort_values("average_distance", ascending=True)
        y = compute(sorted_data["average_distance"])[0]
        x = compute(sorted_data["average_contact_probability"])[0]
        return SortedData(x=x, y=y)  # Return as SortedData namedtuple

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
