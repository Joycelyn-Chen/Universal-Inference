import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

# python csv2boxplot.py --csv_file "results.csv"
parser = argparse.ArgumentParser(description="Creating boxplots from CSV simulation results.")
parser.add_argument("--csv_file", type=str, help='Path to the CSV simulation results file.')
args = parser.parse_args()

# 1) load data
df = pd.read_csv(args.csv_file)

# 2) label each beta combination
df['beta_combo'] = df.apply(
    lambda r: f"β1={r.beta1},β2={r.beta2},β3={r.beta3}",
    axis=1
)

# 3) compute mean and IQR per (n, k, beta_combo)
summary = (
    df
    .groupby(['n','k','beta_combo'])['LC']
    .agg(
        mean   = 'mean',
        q25    = lambda x: np.percentile(x, 25),
        q75    = lambda x: np.percentile(x, 75)
    )
    .reset_index()
)

# 4) compute low/high errors
summary['err_low']  = summary['mean'] - summary['q25']
summary['err_high'] = summary['q75']  - summary['mean']

# 4a) clamp to zero so no negatives
summary['err_low']  = np.maximum(summary['err_low'],  0)
summary['err_high'] = np.maximum(summary['err_high'], 0)

# 5) debugging print
for _, row in summary.iterrows():
    print(
        f"n={int(row.n)}, k={int(row.k)}, {row.beta_combo} | "
        f"mean={row.mean}, "
        f"err_low={row.err_low}, "
        f"err_high={row.err_high}"
    )

# 6) fix ordering of the 9 combos
betas = [0.1, 1.0, 10.0]
combo_list   = [(b1,b2,b3) for b1 in betas for b2 in betas for b3 in betas]
combo_labels = [f"β1={b1},β2={b2},β3={b3}" for b1,b2,b3 in combo_list]

# number of bars per β₁ group
chunk = len(betas)**2  # 3×3 = 9

# 7) plot grid of bar+errorbar (mean ± IQR), y in log scale
ks = [10, 50, 100, 500, 2000, 5000]
ns = [100, 1000, 10000]
n_rows, n_cols = len(ks), len(ns)

fig, axes = plt.subplots(n_rows, n_cols,
                         figsize=(4*n_cols, 3*n_rows),
                         sharey=True)

for i, k_val in enumerate(ks):
    for j, n_val in enumerate(ns):
        ax = axes[i, j]
        sub = summary[(summary.k == k_val) & (summary.n == n_val)]
        
        # ensure the 9 combos in the right order
        sub = sub.set_index('beta_combo').reindex(combo_labels).reset_index()
        
        x = np.arange(len(combo_labels))
        means    = sub['mean'].values
        err_low  = sub['err_low'].values
        err_high = sub['err_high'].values
        
        # draw bars with errorbars
        ax.bar(x, means, color='lightgray', edgecolor='k')
        ax.errorbar(x, means, yerr=[err_low, err_high],
                    fmt='none', ecolor='black', capsize=3)
        
        # **set log scale on y-axis**
        ax.set_yscale('log')

        # draw two vertical dashed lines at the β₁ boundaries
        for boundary in [chunk - 0.5, 2*chunk - 0.5]:
            ax.axvline(boundary, color='black', linestyle='--')
        
        # formatting
        ax.set_title(f"k={k_val}, n={n_val}")
        if i == n_rows - 1:
            ax.set_xticks(x)
            ax.set_xticklabels(combo_labels, rotation=90, fontsize=6)
        else:
            ax.set_xticks([])
        if j == 0:
            ax.set_ylabel("log LC")         # mean ± IQR
        
# 8) global legend
fig.legend(combo_labels, title="beta combinations",
           loc="upper right", bbox_to_anchor=(1.15,0.95))

fig.tight_layout(rect=[0, 0, 0.88, 1])
plt.savefig("LC_summary_bargrid_logscale.png", dpi=150)
plt.show()



