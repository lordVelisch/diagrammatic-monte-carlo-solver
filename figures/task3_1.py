import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def p(x_inp):
    return np.exp(-x_inp)/(1-np.exp(-5))

df = pd.read_csv("../data/task3/cdf_inverse_data.csv")

x_start = 0
x_end = 5

x = np.linspace(x_start, x_end, 1000)

fig, ax = plt.subplots()

ax.hist(df["value"], bins=50, color="tab:blue", alpha=0.6, label="samples", density=True)
ax.plot(x, p(x), color="tab:red", label="exact (normalized) $Q(\\tau)$")

ax.set_xlabel(r"$\tau$")
ax.set_ylabel("probability density")
ax.legend()

plt.tight_layout()
plt.savefig("../report/figures/task3/cdf_histogram.pdf")
plt.show()
