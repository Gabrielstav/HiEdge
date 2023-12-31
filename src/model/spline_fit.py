# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass

from scipy.interpolate import PchipInterpolator, UnivariateSpline
from sklearn.isotonic import IsotonicRegression
from src.setup.config_loader import Config
from src.setup.data_structures import SplineInput
import dask.dataframe as dd
import numpy as np


class SplineDataPreparer:

    def __init__(self, binned_data):
        self.binned_data = binned_data

    def prepare_data_for_spline(self):
        x = self.binned_data['average_distance']
        y = self.binned_data['average_contact_probability']
        yerr = self.calculate_standard_error()

        return x, y, yerr

    def calculate_standard_error(self):
        standard_error = (self.binned_data["interaction_std"]
                          / np.sqrt(self.binned_data["interaction_count"]))
        return standard_error


class SplineFitter:

    def __init__(self, x, y, yerr, config: Config):
        self.x = x
        self.y = y
        self.yerr = yerr
        self.config = config

    def fit_spline(self):
        # Check if x values are monotonically decreasing
        if not np.all(np.diff(self.x) < 0):
            raise ValueError("x values must be in monotonically decreasing order")

        # Setting the spline error - it could be set as min(y)^2 or another suitable value
        spline_error = min(self.y) ** 2

        # Fitting the spline
        spline = UnivariateSpline(self.x, self.y, s=spline_error)

        return spline

    def evaluate_spline(self, spline):
        # Evaluate the spline over the range of x values
        spline_y = spline(self.x)

        return spline_y

    def fit_with_isotonic_regression(self, spline):
        # Evaluate the original spline
        original_spline_y = spline(self.x)

        # Apply isotonic regression to ensure monotonically decreasing curve
        ir = IsotonicRegression(increasing=False)
        corrected_spline_y = ir.fit_transform(self.x, original_spline_y)

        # Create a new spline based on the isotonic regression output
        corrected_spline = PchipInterpolator(self.x, corrected_spline_y)

        return corrected_spline

    def multiple_passes_fit(self, num_passes=Config.pipeline_settings.spline_passes):
        spline = self.fit_spline()
        for _ in range(num_passes - 1):
            spline = self.fit_with_isotonic_regression(spline)
        return spline

    def calculate_residuals(self, spline):
        spline_y = spline(self.x)
        residuals = self.y - spline_y
        return residuals

    def calculate_mse(self, residuals):
        mse = np.mean(residuals ** 2)
        return mse

    def calculate_r_squared(self):
        # Assuming self.y has no zero variance
        total_variance = np.var(self.y)
        residuals = self.calculate_residuals(self.fit_spline())
        residual_variance = np.var(residuals)
        r_squared = 1 - (residual_variance / total_variance)
        return r_squared