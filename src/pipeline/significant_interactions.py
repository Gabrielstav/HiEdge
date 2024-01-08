# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from src.setup.config_loader import Config


# TODO:
#   This is the "main" class in the model part of the pipeline.
#   We call this after filtering is complete, and it handles the inter- and intra data accoringly.
#   So if inter, we need statistics and then calculate the p/q values and confidence intervals.
#   If intra, we do the statistics and the input to spline, then spline, and then p/q value and confidence intervals.

# TODO:
#   So we need to:
#    Make bias reader class
#    Make bias applier class
#    Make spline input class (DONE)
#    Make spline fit class (DONE)
#    Make p/q value and confidence interval class
#    Make the same for inter

# Either we handle inter- and intra data in the same file, where the p/q value and confidence interval classes handle both cases, or we make separate files ocntaining the inter- and intra classes.

# if intra and binom:
#   spline and p/q vals
# if inter and binom:
#   inter stats, then binom (?) and p/q vals
# then same checks for NCHG if we implement that, need different classes for that, or we can include the calculations in spline input


