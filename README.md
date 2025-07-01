---
editor_options: 
  markdown: 
    wrap: 100
---

# 2025-topic-03-group-04

# 🧬 Proteome-wide Screen for RNA-dependent Proteins

Welcome to the **Proteome-wide Screen for RNA-dependent Proteins** project! This repository will
serve as the central place for exploring, analyzing, and visualizing data related to RNA-protein
interactions across the proteome.

> ⚠️ *Note: This README is a starting template. Please update it as your project evolves.*
>
> For inspiration on writing a comprehensive and engaging README, check out the [Awesome
> README](https://github.com/matiassingers/awesome-readme?tab=readme-ov-file) repository.

# 📚 Papers

# Reviews

-   [Sternburg et al., Global Approaches in Studying RNA-Binding Protein Interaction Networks, 2020,
    Trends in Biochemical
    Sciences.pdf](https://github.com/user-attachments/files/19981693/Sternburg.et.al.Global.Approaches.in.Studying.RNA-Binding.Protein.Interaction.Networks.2020.Trends.in.Biochemical.Sciences.pdf)
-   [Corley et al., How RNA-Binding Proteins Interact with RNA Molecules and Mechanisms, 2020,
    Molecular
    Cell.pdf](https://github.com/user-attachments/files/19981705/Corley.et.al.How.RNA-Binding.Proteins.Interact.with.RNA.Molecules.and.Mechanisms.2020.Molecular.Cell.pdf)
-   [Gebauer et al., RNA-binding proteins in human genetic disease, 2020, Nature Reviews
    Genetics.pdf](https://github.com/user-attachments/files/19981707/Gebauer.et.al.RNA-binding.proteins.in.human.genetic.disease.2020.Nature.Reviews.Genetics.pdf)

# Experimental methods

-   [Caudron-Herger et al., R-DeeP Proteome-wide and Quantitative Identification of RNA-Dependent
    Proteins by Density Gradient Ultracentrifugation, 2019, Molecular
    Cell.pdf](https://github.com/user-attachments/files/19981712/Caudron-Herger.et.al.R-DeeP.Proteome-wide.and.Quantitative.Identification.of.RNA-Dependent.Proteins.by.Density.Gradient.Ultracentrifugation.2019.Molecular.Cell.pdf)
-   [Caudron-Herger-Identification, quantification and bioinformatic analysis of RNA-dependent
    proteins by RNase treatment and density gradient ultracentrifugation using R-DeeP-2020-Nature
    Protocols_1.pdf](https://github.com/user-attachments/files/19981715/Caudron-Herger-Identification.quantification.and.bioinformatic.analysis.of.RNA-dependent.proteins.by.RNase.treatment.and.density.gradient.ultracentrifugation.using.R-DeeP-2020-Nature.Protocols_1.pdf)
-   [Rajagopal-Proteome-Wide Identification of RNA-Dependent Proteins in Lung Cancer
    Cells-2022-Cancers.pdf](https://github.com/user-attachments/files/19981723/Rajagopal-Proteome-Wide.Identification.of.RNA-Dependent.Proteins.in.Lung.Cancer.Cells-2022-Cancers.pdf)
-   [Rajagopal et al., An atlas of RNA-dependent proteins in cell division reveals the
    riboregulation of mitotic protein-protein interactions. Nat. Commun. 16, 2325
    (2025).pdf](https://github.com/user-attachments/files/19981728/Rajagopal.et.al.An.atlas.of.RNA-dependent.proteins.in.cell.division.reveals.the.riboregulation.of.mitotic.protein-protein.interactions.Nat.Commun.16.2325.2025.pdf)

# **Group 4 Data Analysis Projekt**

#Silent Signals: Using Protein Traits to Predict RNA Interactions/The Shape of Affinity: Exploring Protein Features Behind RNA Binding (Vorschläge von Jetty)

# **Projekt Überblick**

```mermaid
flowchart LR
    A[Mass Spec Daten] --> B[Data Exploration]
    B --> C[Daten Aufreinigen]
    C --> D[Daten Beschreibung]
    D --> E[Daten Normalisierung]
    E --> F[Dimensionsreduktion]
    F --> G[Lineare Regression]
```


## **Ziel:**

## **Libraries:**

## **Daten laden**

# **Daten Aufreinigung**

## **Ziel:**

# **Daten Beschreibung (Data Description)**

## **Ziel:**

# **Daten Normalisierung**

## **Ziel:**

-   Vergleichbarkeit herstellen:
    -   Methoden???
-   Verzerrungen entfernen:
    -   Batch-Effekt entfernen
-   Zahlenbereiche angleichen:
    -   Skalierung und Normierung
    -   Transformation

## **Ablauf:**

```mermaid
flowchart LR
    A[Normalisierung pro Fraktion] --> B[Mittelwertfilter]
    B --> C[Normalisierung pro Protein]
    C --> D[Mittelwertberechung]
    C --> E[Batchentfernung]
```

1.  **3 Replikate pro Fraktion werden normalisiert**:\
    Unterschiede zwischen den Replikaten werden basierend auf dem\
    ähnlichsten Replikatpaar (Normalisierungsfaktor) ausgeglichen\
    =\> Fraktions- und replikatspezifisch skaliert

    **Output**:\
    Zwei Dataframes (für Control und RNase) mit 75 Spalten (3×25)\
    mit den skalierten Intensitäten der Proteine\
    **Normed Control Dataframe (first 6 rows)**\
    ![Mein Screenshot](image/screenshot_normed_control.png)

    **Normed RNase Dataframe (first 6 rows)**\
    ![Mein Screenshot](image/screenshot_normed_rnase.png)

2.  **Glättung der Werte durch einen gleitenden Mittelwert** der gleitende Mittelwert resultiert in
    einer Glättung der mittleren\
    Fraktionen durch den Durchschnitt mit linker und rechter Nachbarfraktion =\> geglättete Werte
    pro Protein und Replikat

    **Output**:\
    Zwei Dataframes (für Control und RNase) mit 75 Spalten (3x25)\
    mit den geglätteten Werten

3.  


