import os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))

pairs = [
    (
        os.path.join(BASE, "CHNOSZ_Pourbaix_Fe.png"),
        os.path.join(BASE, "Reaktoro_Pourbaix_Fe_CHNOSZmatched.png"),
        "Fe Pourbaix (Matched Setup)",
    ),
    (
        os.path.join(BASE, "CHNOSZ_Mosaic_Fe.png"),
        os.path.join(BASE, "Reaktoro_Mosaic_Fe_CHNOSZmatched.png"),
        "Fe Mosaic (Matched Setup)",
    ),
]


def load_img(path):
    img = plt.imread(path)
    if img.ndim == 2:
        img = np.stack([img, img, img], axis=-1)
    if img.shape[-1] == 4:
        img = img[..., :3]
    return img.astype(np.float32)


def resize_nn(img, h, w):
    y = np.linspace(0, img.shape[0] - 1, h).astype(int)
    x = np.linspace(0, img.shape[1] - 1, w).astype(int)
    return img[np.ix_(y, x)]


rows = len(pairs)
fig, axs = plt.subplots(rows, 2, figsize=(13, 5.4 * rows), constrained_layout=True)
if rows == 1:
    axs = np.array([axs])

report_lines = []
for i, (ch_path, rk_path, title) in enumerate(pairs):
    ch = load_img(ch_path)
    rk = load_img(rk_path)

    h = min(ch.shape[0], rk.shape[0])
    w = min(ch.shape[1], rk.shape[1])
    ch_r = resize_nn(ch, h, w)
    rk_r = resize_nn(rk, h, w)

    mse = float(np.mean((ch_r - rk_r) ** 2))
    mae = float(np.mean(np.abs(ch_r - rk_r)))

    axs[i, 0].imshow(ch)
    axs[i, 0].set_title(f"CHNOSZ {title}")
    axs[i, 0].axis("off")

    axs[i, 1].imshow(rk)
    axs[i, 1].set_title(f"Reaktoro {title}")
    axs[i, 1].axis("off")

    report_lines.append(f"{title}: MAE={mae:.4f}, MSE={mse:.4f}")

out_panel = os.path.join(BASE, "CHNOSZ_vs_Reaktoro_matched_panels.png")
fig.suptitle("CHNOSZ vs Reaktoro (Matched Setup)", fontsize=16)
fig.savefig(out_panel, dpi=180)
plt.close(fig)

out_txt = os.path.join(BASE, "CHNOSZ_vs_Reaktoro_matched_metrics.txt")
with open(out_txt, "w", encoding="utf-8") as f:
    f.write("\n".join(report_lines) + "\n")

print("Wrote:", out_panel)
print("Wrote:", out_txt)
for line in report_lines:
    print(line)
