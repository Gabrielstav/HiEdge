# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from src.setup.config_loader import Config


# TODO:
#   This is the "main" class in the model part of the pipeline.
#   We call this after filtering is complete, and it handles the inter- and intra data accoringly.
#   So if inter, we need statistics and then calculate the p/q values and confidence intervals.
#   If intra, we do the statistics and the input to spline, then spline, and then p/q value and confidence intervals.


class CalculateSignificance:

    def __init__(self, config: Config, data):
        self.config = config
        self.data = data

    def run(self):
        pass

    def _calculate_p_value_intra(self):
        pass

    def _calculate_p_value_inter(self):
        pass
