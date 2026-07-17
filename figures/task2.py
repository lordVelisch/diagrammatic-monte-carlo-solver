import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("../data/task2/pi_estimation_convergence.csv")

area = df[df["method"] == "area"]
integration = df[df["method"] == "integration"]

plt.loglog(area["N"], area["pi_estimate"], label="area method", marker="x", linestyle="None", markersize=1)

plt.xlabel("N")
plt.ylabel("Pi estimate")

plt.axhline(np.pi, label="pi", color="red")
plt.legend()
plt.tight_layout()
plt.show()


