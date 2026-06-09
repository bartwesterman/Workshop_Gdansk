# Workshop_Gdansk_2026

## Introduction

In this workshop you will use **Jupyter Notebook** to analyze public cancer dependency data from the **DepMap Project**.

The objective is to identify transcription factors that are important for the survival of specific cancer types and investigate their role in cancer biology.

No prior programming experience is required.

---

# Step 1: Install Anaconda

Download and install Anaconda:

https://www.anaconda.com/download

After installation, open the **Anaconda Prompt**.

Create a new Python environment:

```bash
conda create -n depmap python=3.12
```

Activate the environment:

```bash
conda activate depmap
```

Install the required packages:

```bash
pip install jupyter pandas matplotlib
```

Start Jupyter Notebook:

```bash
jupyter notebook
```

Create a new notebook:

```text
New → Python 3
```

Run notebook cells using:

```text
Shift + Enter
```

---

# Step 2: Download the DepMap Data

For this workshop you need two files:

```text
CRISPRGeneEffect.csv
Model.csv
```

These files are available from the instructor.

Copy the following files:

```text
CRISPRGeneEffect.csv
Model.csv
```

Create the following folder structure:

```text
Workshop_Gdansk/
├── data/
├── results/
└── notebook.ipynb
```

Place the downloaded files inside the `data` folder:

```text
Workshop_Gdansk/
├── data/
│   ├── CRISPRGeneEffect.csv
│   └── Model.csv
├── results/
└── notebook.ipynb
```

---

# Step 3: Verify the Installation

Create a new notebook cell and run:

```python
from pathlib import Path

print(Path("data/CRISPRGeneEffect.csv").exists())
print(Path("data/Model.csv").exists())
```

Expected output:

```text
True
True
```

If both values are `True`, you are ready to start the analysis.

---

# Step 4: Load the DepMap Data

Run:

```python
import pandas as pd

gene_effect = pd.read_csv("data/CRISPRGeneEffect.csv")
metadata = pd.read_csv("data/Model.csv")

print("Gene Effect data:", gene_effect.shape)
print("Model metadata:", metadata.shape)
```

This may take a few seconds because the gene-effect file is large.

---

# Step 5: Display Available Tumor Types

Run:

```python
tumor_types = sorted(
    metadata["OncotreeLineage"]
    .dropna()
    .unique()
)

for i, tumor in enumerate(tumor_types):
    print(f"{i}: {tumor}")
```

You will see a list similar to:

```text
0: Acute Myeloid Leukemia
1: Bladder Cancer
2: Breast Cancer
3: Colorectal Cancer
4: Glioma
...
```

---

# Step 6: Select a Tumor Type

Choose one tumor type:

```python
selection = int(input("Enter tumor type number: "))

selected_tumor = tumor_types[selection]

print("Selected tumor type:", selected_tumor)
```

---

# Step 7: Enter Up to 10 Transcription Factors

Enter up to ten transcription factors separated by commas:

```python
tf_input = input(
    "Enter transcription factors separated by commas: "
)

transcription_factors = [
    x.strip().upper()
    for x in tf_input.split(",")
][:10]

print(transcription_factors)
```

Example:

```text
SOX2,OLIG2,MYC,MYCN,STAT3
```

---

## Step 8: Analyze Average and Model-Level Dependency Scores

This step calculates:

- The average dependency score per transcription factor
- The dependency score for each individual DepMap model
- Tumor subtype information for each model

This is useful because broad groups such as **CNS**, **Glioma**, or **Lung Cancer** often contain multiple biological subtypes.

Run the following code:

```python
selected_models = metadata[
    metadata["OncotreeLineage"] == selected_tumor
].copy()

model_ids = set(selected_models["ModelID"])

tumor_gene_effect = gene_effect[
    gene_effect["Unnamed: 0"].isin(model_ids)
].copy()

summary_results = []
model_level_results = []

for tf in transcription_factors:

    matching_columns = [
        c for c in tumor_gene_effect.columns
        if c.startswith(tf + " ")
    ]

    if len(matching_columns) == 0:

        summary_results.append({
            "TranscriptionFactor": tf,
            "Status": "Not Found",
            "MeanGeneEffect": None,
            "MedianGeneEffect": None,
            "NumberOfCellLines": 0
        })

    else:

        col = matching_columns[0]

        scores = tumor_gene_effect[
            ["Unnamed: 0", col]
        ].copy()

        scores = scores.rename(
            columns={
                "Unnamed: 0": "ModelID",
                col: "GeneEffect"
            }
        )

        scores["TranscriptionFactor"] = tf

        scores = scores.merge(
            selected_models[
                [
                    "ModelID",
                    "CellLineName",
                    "OncotreeLineage",
                    "OncotreePrimaryDisease",
                    "OncotreeSubtype"
                ]
            ],
            on="ModelID",
            how="left"
        )

        model_level_results.append(scores)

        summary_results.append({
            "TranscriptionFactor": tf,
            "Status": "Found",
            "DepMapColumn": col,
            "MeanGeneEffect": scores["GeneEffect"].mean(),
            "MedianGeneEffect": scores["GeneEffect"].median(),
            "NumberOfCellLines": scores["GeneEffect"].notna().sum()
        })

summary_results = pd.DataFrame(summary_results)

summary_results = summary_results.sort_values(
    "MeanGeneEffect",
    ascending=True
)

model_level_results = pd.concat(
    model_level_results,
    ignore_index=True
)

model_level_results = model_level_results.sort_values(
    ["TranscriptionFactor", "GeneEffect"],
    ascending=[True, True]
)

print("Summary Results")
display(summary_results)

print("Model-Level Results")
display(model_level_results)
```

## Step 9: Save the Results

After running Step 8, save both the summary results and the model-level results.

Run the following code in a new Jupyter Notebook cell:

```python
from pathlib import Path

# Create results folder if it does not exist
Path("results").mkdir(exist_ok=True)

# Create a safe file name
safe_name = (
    selected_tumor
    .replace(" ", "_")
    .replace("/", "_")
)

# Output file names
summary_output_file = (
    f"results/{safe_name}_TF_dependency_summary.csv"
)

model_output_file = (
    f"results/{safe_name}_TF_dependency_model_level.csv"
)

# Save files
summary_results.to_csv(
    summary_output_file,
    index=False
)

model_level_results.to_csv(
    model_output_file,
    index=False
)

print("Summary results saved to:")
print(summary_output_file)

print()

print("Model-level results saved to:")
print(model_output_file)
```

### Expected Output

For example, if you selected **Glioma**, the script will create:

```text
results/
├── Glioma_TF_dependency_summary.csv
└── Glioma_TF_dependency_model_level.csv
```

### Summary File

The summary file contains one row per transcription factor.

Example:

| TranscriptionFactor | MeanGeneEffect | MedianGeneEffect |
| ------------------- | -------------- | ---------------- |
| SOX2                | -1.34          | -1.29            |
| MYC                 | -0.92          | -0.88            |
| OLIG2               | -0.71          | -0.69            |

### Model-Level File

The model-level file contains one row per cell line.

Example:

| TranscriptionFactor | CellLineName | OncotreePrimaryDisease | GeneEffect |
| ------------------- | ------------ | ---------------------- | ---------- |
| SOX2                | LN229        | Glioblastoma           | -1.75      |
| SOX2                | U87MG        | Glioblastoma           | -1.61      |
| SOX2                | SF268        | Astrocytoma            | -0.52      |

### Why Save Both Files?

The summary file provides an overview of the average dependency of each transcription factor.

The model-level file allows you to investigate:

* Differences between cell lines
* Differences between tumor subtypes
* Outlier models
* Subtype-specific transcription factor dependencies

This information is often more informative than the average score alone.


# Interpreting the Results

DepMap Gene Effect scores indicate how important a gene is for cell survival.

| Score  | Interpretation         |
| ------ | ---------------------- |
| 0      | No effect              |
| -0.5   | Moderate dependency    |
| -1.0   | Strong dependency      |
| < -1.0 | Very strong dependency |

More negative values indicate a stronger lethal effect when the gene is knocked out.

---

# Assignment

Select:

* One tumor type
* Up to 10 transcription factors

Answer the following questions:

1. Which transcription factor shows the strongest dependency?
2. Is this transcription factor already known to play a role in this cancer type?
3. Is it involved in the normal development of the tissue from which this cancer originates?
4. Could it represent a potential therapeutic target?

Prepare a short presentation (5 minutes) summarizing your findings.

---

# Useful Resources

* DepMap: https://depmap.org
* PubMed: https://pubmed.ncbi.nlm.nih.gov
* GeneCards: https://www.genecards.org
* Human Protein Atlas: https://www.proteinatlas.org

---

