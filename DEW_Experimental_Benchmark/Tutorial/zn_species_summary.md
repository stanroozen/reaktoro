# DEW17HP622_Zn species scan summary

Database scanned: embedded/databases/perplex/DEW17HP622_Zn_2025-reaktoro.json

## Aqueous Zn species (15 total)

| Name | Formula | Charge | Notes |
|------|---------|--------|-------|
| Zn2+ | Zn | +2 | Free zinc ion |
| ZnOH+ | HOZn | +1 | Hydroxide complex |
| ZnO | OZn | 0 | Neutral hydroxide |
| ZnO2-2 | O2Zn | -2 | Zincate |
| HZnO2- | HO2Zn | -1 | Alkaline Zn complex |
| ZnCl+ | ClZn | +1 | Chloride complex |
| ZnCl2 | Cl2Zn | 0 | Neutral chloride |
| ZnCl3- | Cl3Zn | -1 | Chloride complex |
| ZnCl4-2 | Cl4Zn | -2 | Chloride complex |
| ZnF+ | FZn | +1 | Fluoride complex |
| ZnHCO3+ | CHO3Zn | +1 | Carbonate complex |
| Zn(HS)2 | H2S2Zn | 0 | Neutral bisulfide complex |
| Zn(HS)2OH- | H2.5OS2Zn | -1 | Bisulfide-hydroxide |
| Zn(HS)3- | H1.5S1.5Zn | -1 | Tris-bisulfide |
| Zn(HS)4-2 | HS2Zn | -2 | Tetrakis-bisulfide |

## Zn-bearing minerals present (13 total)

| Abbreviation | Formula | Name |
|-------------|---------|------|
| Frk | Fe2O4Zn | Franklinite |
| Ghn | Al2O4Zn | Gahnite |
| HZnc | H6C2O12Zn5 | Hydrozincite |
| Hrds | Si2Ca2O7Zn | Hardystonite |
| Smth | CO3Zn | Smithsonite |
| Sph | SZn | Sphalerite |
| Wlm | Si2O4Zn | Willemite |
| Wrt | Fe0.5SZn0.5 | Wurtzite |
| Zn | Zn | Native zinc |
| Zn-St | HAl9Si4O24Zn2 | Zn-stilpnomelane |
| ZnSp | TiO4Zn2 | Zn-spinel |
| Znc | OZn | Zincite |
| Znks | SO4Zn | Zinkosite |

## Database fix note
The original DEW17HP622_Zn_2025.dat was missing an `end` line after the `HZnO2-` entry, causing the converter to silently absorb `Zn(HS)2` into the preceding block. This was fixed and the database regenerated. The JSON now contains 398 species matching all 398 DAT entries.
