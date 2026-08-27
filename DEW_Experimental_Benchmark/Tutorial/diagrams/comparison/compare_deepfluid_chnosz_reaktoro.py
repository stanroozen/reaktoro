import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

BASE = os.path.abspath(os.path.dirname(__file__))

img_rkt = os.path.join(BASE, "Reaktoro_DeepFluid_LogfO2_pH_Fe.png")
img_chz = os.path.join(BASE, "CHNOSZ_DeepFluid_LogfO2_pH_Fe.png")
out_panel = os.path.join(BASE, "CHNOSZ_vs_Reaktoro_DeepFluid_LogfO2_pH_panel.png")
out_txt = os.path.join(BASE, "CHNOSZ_vs_Reaktoro_DeepFluid_LogfO2_pH_notes.txt")

if not os.path.isfile(img_rkt):
    raise FileNotFoundError(f"Missing Reaktoro image: {img_rkt}")
if not os.path.isfile(img_chz):
    raise FileNotFoundError(f"Missing CHNOSZ image: {img_chz}")

im_rkt = Image.open(img_rkt).convert("RGB")
im_chz = Image.open(img_chz).convert("RGB")

fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.0), dpi=180)
axes[0].imshow(im_chz)
axes[0].set_title("CHNOSZ deep-fluid logfO2-pH")
axes[0].axis("off")

axes[1].imshow(im_rkt)
axes[1].set_title("Reaktoro deep-fluid logfO2-pH")
axes[1].axis("off")

fig.suptitle("Deep-fluid style Fe diagram comparison (log10(fO2) vs pH)")
plt.tight_layout()
fig.savefig(out_panel, dpi=180)
plt.close(fig)

with open(out_txt, "w", encoding="utf-8") as f:
    f.write("CHNOSZ vs Reaktoro deep-fluid diagram comparison\n")
    f.write("- Space: log10(fO2) vs pH\n")
    f.write("- Conditions target: T=350 C, P=2000 bar\n")
    f.write(
        "- This panel is visual-only (no MAE metric due palette/label differences).\n"
    )
    f.write(f"- CHNOSZ image: {img_chz}\n")
    f.write(f"- Reaktoro image: {img_rkt}\n")
    f.write(f"- Panel: {out_panel}\n")

print("Wrote:", out_panel)
print("Wrote:", out_txt)
