## Download Flow Data

1. Visit [morganrj071/oscillating-hydrofoil-wake-Re11000](https://huggingface.co/datasets/morganrj071/oscillating-hydrofoil-wake-Re11000).
2. Download `oscfoil_data.mat` and `foilcoords.mat`.
3. Create `oscillating_foil_wake/data/`.
4. Place both files in that folder.

The expected layout is:

```text oscillating_foil_wake/README.md
oscillating_foil_wake/
├── data/
│   ├── oscfoil_data.mat
│   └── foilcoords.mat
├── functions/
│   ├── FTLEsb.m
│   ├── CSEsb.m
│   └── optdmdsrc/
├── customcolormap/
├── OptDMD.m
├── calcFTLE.m
├── calcMTU.m
└── README.md
```