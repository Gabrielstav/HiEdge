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
        initial_spline = self._fit_initial_spline()
        self.adjusted_y = self._apply_isotonic_regression(initial_spline)
        initial_spline.get_knots()

        x_spline = initial_spline.get_knots()
        y_spline = initial_spline(x_spline)
        print(f"Initial spline knots: {x_spline}")
        print(f"Initial spline values: {y_spline}")

        # plot spline
        # self._plot_spline()
        # plot initial data
        # self._plot_initial_data()

        # TODO: Just adjust y-values of spline instead of creating new spline
        self.adjusted_spline = UnivariateSpline(self.x, self.adjusted_y, s=self.spline_error)

        y_adusted = interp1d(self.x, self.adjusted_y, kind="cubic")



    def _fit_initial_spline(self):
        x_interpolated = np.linspace(min(self.x), max(self.x), 1000)
        y_cubic = interp1d(self.x, self.y, kind="cubic")
        initial_spline = UnivariateSpline(x_interpolated, y_cubic(x_interpolated), s=self.spline_error)
        self._plot_spline()
        # initial_spline = UnivariateSpline(self.x, self.y, s=self.spline_error)
        return initial_spline

    def _apply_isotonic_regression(self, initial_spline) -> np.ndarray:
        y_predictions = initial_spline(self.x)
        ir = IsotonicRegression(increasing=False)
        adjusted_y = ir.fit_transform(self.x, y_predictions)

        print(f"Initial y: {initial_spline(self.x)}")
        print(f"Adjusted y: {adjusted_y}")

        return adjusted_y

    @staticmethod
    def _get_sorted_data(binned_data: dd.DataFrame) -> SortedData:
        sorted_data = binned_data.sort_values("average_distance", ascending=True)
        x = compute(sorted_data["average_distance"])[0]
        y = compute(sorted_data["average_contact_probability"])[0]
        return SortedData(x=x, y=y)

    def _plot_spline(self):
        ax.plot(self.x, self.y, label="Original Data", alpha=0.5)
        pass

    def _plot_initial_data(self):
        # TODO: DEBUGGING

        print(f"Genomic distances: {self.x}")
        print(f"Probabilities before spline: {self.y}")
        print(f"Length of x: {len(self.x)}")
        print(f"Length of y: {len(self.y)}")

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
