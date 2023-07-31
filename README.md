# interFLOW: maximum flow framework for the identification of factors mediating the signaling convergence of multiple receptors 

Installation:


First put your Seurat obj in the InputObj folder

Update the config csv with :

1.	name of Seurat clusters obj

3.	vector of signal receiving cluster separated by space

5.	vector of signal sending cluster separated by space

7.	Output Plot Path

9.	Project Name

the Graphs obj will be crated in the outputObj folder

Now you can run the FLOW analysis with:

python3 main.py

In order to run the TF analysis first add the cluster annotation to the Annot.csv file.

Now run -> python3 tf_network.py

first enter [i] to get a csv file (will be saved in the files folder) with the importance score of each TF 

now run tf_network.py again and enter [f] in order to crate the csv edge list of the flow to a given TF





