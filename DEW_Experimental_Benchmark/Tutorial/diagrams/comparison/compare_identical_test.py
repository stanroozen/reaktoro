import os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))

ch_path = os.path.join(BASE, "CHNOSZ_Pourbaix_Fe_identical_test.png")
rk_path = os.path.join(BASE, "Reaktoro_Pourbaix_Fe_identical_test.png")

# Pixel-level MAE threshold for "passing" the comparison.
# Diagrams sharing the same colour palette and boundaries should have
# MAE < 0.05 (images are normalised to [0,1]).  The old closed-system
# Reaktoro result had MAE ≈ 0.31 — any value below the threshold confirms
# the open-system fix is working.
MAE_PASS_THRESHOLD = 0.05


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


ch = load_img(ch_path)
rk = load_img(rk_path)

h = min(ch.shape[0], rk.shape[0])
w = min(ch.shape[1], rk.shape[1])
ch_r = resize_nn(ch, h, w)
rk_r = resize_nn(rk, h, w)

mse = float(np.mean((ch_r - rk_r) ** 2))
mae = float(np.mean(np.abs(ch_r - rk_r)))
passed = mae < MAE_PASS_THRESHOLD
status = "PASS" if passed else "FAIL"

fig, axs = plt.subplots(1, 2, figsize=(13, 5.4), constrained_layout=True)
axs[0].imshow(ch)
axs[0].set_title("CHNOSZ — open system\nlog10(a(Fe²⁺)) = 0, T=25°C, P=1 bar")
axs[0].axis("off")
axs[1].imshow(rk)
axs[1].set_title("Reaktoro — open system\nspecs.lgActivity('Fe+2'), log10(a) = 0")
axs[1].axis("off")

out_panel = os.path.join(BASE, "CHNOSZ_vs_Reaktoro_identical_test_panel.png")
fig.suptitle(
    f"CHNOSZ vs Reaktoro — Fe Pourbaix (open system, shared species)\n"
    f"MAE={mae:.4f}  MSE={mse:.4f}  [{status}  threshold={MAE_PASS_THRESHOLD}]",
    fontsize=13,
)
fig.savefig(out_panel, dpi=180)
plt.close(fig)

out_txt = os.path.join(BASE, "CHNOSZ_vs_Reaktoro_identical_test_metrics.txt")
with open(out_txt, "w", encoding="utf-8") as f:
    f.write(f"Identical test (Fe Pourbaix): MAE={mae:.4f}, MSE={mse:.4f}, {status}\n")
    f.write(f"Pass threshold (MAE): {MAE_PASS_THRESHOLD}\n")
    f.write(
        "Setup: open system, log10(a(Fe+2))=0, T=25C, P=1 bar, pH [-2,16], Eh [-2,2] V\n"
    )
    f.write("Bug fixes applied:\n")
    f.write(
        "  FIX 1: specs.lgActivity('Fe+2') — system open to Fe (implicit titrant)\n"
    )
    f.write(
        "  FIX 2: predominance via mineral amount > threshold, else max aqueous log10(a)\n"
    )

print("Wrote:", out_panel)
print("Wrote:", out_txt)
print(f"Identical test (Fe Pourbaix): MAE={mae:.4f}  MSE={mse:.4f}  [{status}]")
