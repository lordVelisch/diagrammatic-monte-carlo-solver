import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("../data/task2/pi_estimation_distribution.csv")

area = df[df["method"] == "area"]
integration = df[df["method"] == "integration"]

def plot_histogram(df_to_plot):
    mean = df_to_plot["pi_estimate"].mean()
    variance = df_to_plot["pi_estimate"].var()

    plt.hist(df_to_plot["pi_estimate"], bins=50, label=f"method: {df_to_plot['method'].iloc[0]}")

    print(f"Mean value is {mean}")
    print(f"Variance is {variance}")

    plt.legend()

plt.axvline(np.pi, label="pi", color="red")

plt.title("N=1_000_000")

plot_histogram(area)
plot_histogram(integration)

plt.tight_layout()
plt.savefig("../report/figures/task2/pi_distribution.pdf")
plt.show()
