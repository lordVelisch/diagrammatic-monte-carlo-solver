import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("../data/task2/pi_estimation_convergence.csv")

area = df[df["method"] == "area"]
integration = df[df["method"] == "integration"]

plt.loglog(area["N"][100:], area["pi_estimate"][100:], label="area method", marker="x", linestyle="-", linewidth=1, markersize=1)
plt.loglog(integration["N"][100:], integration["pi_estimate"][100:], label="integration method", marker="x", linestyle="-", linewidth=1, markersize=1)

plt.xlabel("N")
plt.ylabel("Pi estimate")

plt.axhline(np.pi, label="pi", color="red")
plt.legend()
plt.tight_layout()
plt.savefig("../report/figures/task2/pi_convergence.pdf")
plt.show()


