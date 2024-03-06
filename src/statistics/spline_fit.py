# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from scipy.interpolate import UnivariateSpline, interp1d
from sklearn.isotonic import IsotonicRegression
import dask.dataframe as dd
from dask import compute
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt

SortedData = namedtuple("SortedData", ["x", "y"])

class SplineFitter:
    def __init__(self, binned_data: dd.DataFrame, config):
        self.adjusted_spline = None
        self.adjusted_y = None
        self.binned_data = binned_data
        self.config = config
        self.x, self.y = self._get_sorted_data(binned_data)
        self.spline_error = min(self.y) ** 2

    def run(self):

        # TODO: DEBUGGING
        plt.figure(figsize=(12, 6))

        plt.subplot(1, 2, 1)
        plt.hist(self.x, bins=50)
        plt.title('Distribution of Genomic Distances (x)')
        plt.xlabel('Genomic Distance')
        plt.ylabel('Frequency')

        plt.subplot(1, 2, 2)
        plt.hist(self.y, bins=50)
        plt.title('Distribution of Probabilities (y)')
        plt.xlabel('Probability')
        plt.ylabel('Frequency')

        plt.tight_layout()
        plt.show()

        print(f'Min X: {min(self.x)}, Max X: {max(self.x)}')
        print(f'Min Y: {min(self.y)}, Max Y: {max(self.y)}')

        plt.scatter(self.x, self.y, alpha=0.5)
        plt.title('Genomic Distance vs Probability')
        plt.xlabel('Genomic Distance')
        plt.ylabel('Probability')
        plt.show()

        initial_spline = self._fit_initial_spline()
        self.adjusted_y = self._apply_isotonic_regression(initial_spline)
        # Create an interpolation function based on adjusted y-values
        self.adjusted_spline = interp1d(self.x, self.adjusted_y, kind="linear", fill_value="extrapolate")

    def get_adjusted_spline(self):
        # Ensure run has been called
        if hasattr(self, "adjusted_spline"):
            return self.adjusted_spline
        else:
            raise ValueError("Spline has not been fitted yet.")

    def _fit_initial_spline(self):
        return UnivariateSpline(self.x, self.y, s=self.spline_error)

    def _apply_isotonic_regression(self, initial_spline) -> np.ndarray:
        y_predictions = initial_spline(self.x)
        ir = IsotonicRegression(increasing=False)
        adjusted_y = ir.fit_transform(self.x, y_predictions)
        return adjusted_y

    @staticmethod
    def _get_sorted_data(binned_data: dd.DataFrame) -> SortedData:
        sorted_data = binned_data.sort_values("average_distance", ascending=True)
        x = compute(sorted_data["average_distance"])[0]
        y = compute(sorted_data["average_contact_probability"])[0]
        return SortedData(x=x, y=y)

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
