# interFLOW: Maximum Flow Framework for Identifying Factors Mediating Signaling Convergence of Multiple Receptors

## Introduction

Motivation: Cell-cell crosstalk involves simultaneous interactions of multiple receptors and ligands, followed by downstream signaling cascades working through receptors converging at dominant transcription factors which then integrate and propagate multiple signals into a cellular response. Single cell RNAseq of multiple cell subsets isolated from a defined microenvironment provides us with a unique opportunity to learn about such interactions reflected in their gene expression levels.

Results: We developed the interFLOW framework to map the potential ligand-receptor interactions between different cell subsets based on a maximum flow computation in a network of protein-protein interactions (PPIs). The maximum flow approach further allows characterization of the intracellular downstream signal transduction from differentially expressed receptors towards dominant transcription factors, therefore, enabling the association between a set of receptors and their downstream activated pathways. Importantly, we were able to identify key transcription factors toward which the convergence of multiple receptor signaling occurs. These identified factors have a unique role in the integration and propagation of signaling following specific cell-cell interactions.

## Installation (Conda)

To use interFLOW with Conda, follow these steps:

1. Clone the Git repository to your local machine.
2. Create a new Conda environment using the provided `interFLOW\files\env.yaml` file or using pip requirements.txt.
3. Execute the `packages.R` file using R v4.1.1 or a later version to update all relevant packages. 

## Running the Pipeline

Before running the pipeline, ensure that your Seurat object is placed in the `InputObj` folder. Then, update the configuration CSV with the following information:

1. **Name of Seurat clusters object**
2. **Vector of signal receiving clusters (Receptors)**: Provide the names of clusters separated by spaces. The cluster names should be saved in the active ident field of the Seurat object.
3. **Vector of signal-sending clusters (Ligand)**: Provide the names of clusters separated by spaces.
4. **Project Name**

To run the pipeline, execute the `main.py` file.

 the Graphs object will be created in the `outputObj` folder, and the resulting plots will be saved in the `plotOut` folder.

## TF Analysis (Post-Basic Analysis)

The TF analysis can be performed after completing the basic analysis. It focuses on a specific interaction involving one signal-sending cluster and one signal-receiving cluster.

Follow these steps for TF analysis:

1. Run the `tf_network.py` file.
2. Enter `[i]` to generate a CSV file with the importance scores of each Transcription Factor (TF). The file will be saved in the `files` folder.
3. Run `tf_network.py` again and enter `[f]` to create the CSV edge list representing the flow to a given TF.

With the TF analysis, you can gain valuable insights into the role of individual Transcription Factors in the signaling convergence process.

