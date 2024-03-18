# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass

# TODO:
#   Make plots of initial distance/frequency data
#   Overay initial spline fit
#   Overlay spline fit after bias correction and with error estimates, confidence intervals etc


# plt.figure(figsize=(12, 6))
#
# plt.subplot(1, 2, 1)
# plt.hist(self.x, bins=50)
# plt.title('Distribution of Genomic Distances (x)')
# plt.xlabel('Genomic Distance')
# plt.ylabel('Frequency')
#
# plt.subplot(1, 2, 2)
# plt.hist(self.y, bins=50)
# plt.title('Distribution of Probabilities (y)')
# plt.xlabel('Probability')
# plt.ylabel('Frequency')
#
# plt.tight_layout()
# plt.show()
#
# print(f'Min X: {min(self.x)}, Max X: {max(self.x)}')
# print(f'Min Y: {min(self.y)}, Max Y: {max(self.y)}')
#
# plt.scatter(self.x, self.y, alpha=0.5)
# plt.title('Genomic Distance vs Probability')
# plt.xlabel('Genomic Distance')
# plt.ylabel('Probability')
# plt.show()