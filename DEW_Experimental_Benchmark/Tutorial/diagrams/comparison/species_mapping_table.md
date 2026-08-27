# Species Mapping: Reaktoro vs CHNOSZ

## Fe_Pourbaix
- exact: 11, approx: 1, missing_in_chnosz: 2, missing_in_reaktoro: 0

| Reaktoro | CHNOSZ | Status | Reason |
|---|---|---|---|
| Fe+2 | Fe+2 | exact | case-insensitive exact name match |
| Fe+3 | Fe+3 | exact | case-insensitive exact name match |
| FeO(aq) |  | missing_in_chnosz | no CHNOSZ species with equivalent normalized key |
| FeO+ | FeO+ | exact | case-insensitive exact name match |
| FeO2- | FeO2- | exact | case-insensitive exact name match |
| FeOH+ | FeOH+ | exact | case-insensitive exact name match |
| FeOH+2 | FeOH+2 | exact | case-insensitive exact name match |
| HFeO2(aq) | HFeO2 | approx | normalized alias (aq polymorph/comma/punctuation removed) |
| HFeO2- | HFeO2- | exact | case-insensitive exact name match |
| Ferropericlase |  | missing_in_chnosz | no CHNOSZ species with equivalent normalized key |
| Goethite | goethite | exact | case-insensitive exact name match |
| Hematite | hematite | exact | case-insensitive exact name match |
| Iron | iron | exact | case-insensitive exact name match |
| Magnetite | magnetite | exact | case-insensitive exact name match |

## Fe_Mosaic
- exact: 7, approx: 1, missing_in_chnosz: 0, missing_in_reaktoro: 0

| Reaktoro | CHNOSZ | Status | Reason |
|---|---|---|---|
| Fe+2 | Fe+2 | exact | case-insensitive exact name match |
| Fe+3 | Fe+3 | exact | case-insensitive exact name match |
| HFeO2- | HFeO2- | exact | case-insensitive exact name match |
| Pyrite | pyrite | exact | case-insensitive exact name match |
| Pyrrhotite,trot | pyrrhotite | approx | normalized alias (aq polymorph/comma/punctuation removed) |
| Siderite | siderite | exact | case-insensitive exact name match |
| Hematite | hematite | exact | case-insensitive exact name match |
| Magnetite | magnetite | exact | case-insensitive exact name match |
