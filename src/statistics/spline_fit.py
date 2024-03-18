# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from scipy.interpolate import UnivariateSpline, interp1d
from sklearn.isotonic import IsotonicRegression
import dask.dataframe as dd
from dask import compute
import numpy as np
from src.setup.config_loader import Config
from collections import namedtuple
import matplotlib.pyplot as plt

SortedData = namedtuple("SortedData", ["x", "y"])

# TODO: Figure out why the old implementation didn't work and make the new one cleaner
#   add plots in plots/spline_fitting
#   pre-sort data (or at least sort it once)
#   make spline, error, x and y instance variables?

class SplineFitter:
    def __init__(self, data, config):
        self.data = data
        self.config = config
        self.spline = None
        self.spline_error = None

    def fit_spline(self):
        # Ensure the DataFrame is computed if it's a Dask DataFrame to work with in-memory arrays
        if isinstance(self.data, dd.DataFrame):
            sorted_data = self.data.compute().sort_values("average_genomic_distance")
        else:
            sorted_data = self.data.sort_values("average_genomic_distance")

        x = sorted_data["average_genomic_distance"].values
        y = sorted_data["average_contact_probability"].values

        print(f"X values: {x}")
        print(f"Y values: {y}")

        spline_error = min(y) ** 2
        self.spline_error = spline_error

        # rescale y-values maybe??
        # y_rescaled = (y - y.min()) / (y.max() - y.min()) * (1 - 1e-6) + 1e-6

        self.spline = UnivariateSpline(x, y, s=spline_error)
        return self

    def enforce_monotonicity(self):
        # again, sort the data
        if isinstance(self.data, dd.DataFrame):
            sorted_data = self.data.compute().sort_values("average_genomic_distance")
        else:
            sorted_data = self.data.sort_values("average_genomic_distance")
        x = sorted_data["average_genomic_distance"].values

        # print the x values
        print(f"X values: {x}")

        # print the original y values
        print(f"Y values: {self.spline(x)}")

        y = self.spline(x)  # Evaluate spline on the original x values

        iso_reg = IsotonicRegression(increasing=True, y_min=min(y), y_max=max(y))
        y_mono = iso_reg.fit_transform(x, y)

        # print the monotonic y values
        print(f"Monotonic Y values: {y_mono}")

        self.spline = UnivariateSpline(x, y_mono, s=self.spline_error)
        return self

    def predict(self, x_new):
        return self.spline(x_new)


    

# class SplineFitter:
#     def __init__(self, binned_data: dd.DataFrame, config):
#         self.adjusted_spline = None
#         self.adjusted_y = None
#         self.binned_data = binned_data
#         self.config = config
#         self.x, self.y = self._get_sorted_data(binned_data)
#         self.spline_error = min(self.y) ** 2
#
#
#
#     def run(self):
#         initial_spline = self._fit_initial_spline()
#         self.adjusted_y = self._apply_isotonic_regression(initial_spline)
#         initial_spline.get_knots()
#
#         x_spline = initial_spline.get_knots()
#         y_spline = initial_spline(x_spline)
#         print(f"Initial spline knots: {x_spline}")
#         print(f"Initial spline values: {y_spline}")
#
#         # plot spline
#         # self._plot_spline()
#         # plot initial data
#         # self._plot_initial_data()
#
#         # TODO: Just adjust y-values of spline instead of creating new spline
#         self.adjusted_spline = UnivariateSpline(self.x, self.adjusted_y, s=self.spline_error)
#
#         y_adusted = interp1d(self.x, self.adjusted_y, kind="cubic")
#
#
#
#     def _fit_initial_spline(self):
#         x_interpolated = np.linspace(min(self.x), max(self.x), 1000)
#         y_cubic = interp1d(self.x, self.y, kind="cubic")
#         initial_spline = UnivariateSpline(x_interpolated, y_cubic(x_interpolated), s=self.spline_error)
#         return initial_spline
#
#     def _apply_isotonic_regression(self, initial_spline) -> np.ndarray:
#         y_predictions = initial_spline(self.x)
#         ir = IsotonicRegression(increasing=False)
#         adjusted_y = ir.fit_transform(self.x, y_predictions)
#
#         print(f"Initial y: {initial_spline(self.x)}")
#         print(f"Adjusted y: {adjusted_y}")
#
#         return adjusted_y
#
#     @staticmethod
#     def _get_sorted_data(binned_data: dd.DataFrame) -> SortedData:
#         sorted_data = binned_data.sort_values("average_distance", ascending=True)
#         x = compute(sorted_data["average_genomic_distance"])[0]
#         y = compute(sorted_data["average_contact_probability"])[0]
#         return SortedData(x=x, y=y)


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
