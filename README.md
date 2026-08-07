# driver-gene-net
## background
Cancer is the result of cells growing uncontrollably because their DNA has been damaged. When tumour DNA is sequenced and compared to pt's healthy tissue, researchers find somatic mutations, which are changes in human DNA that occur after conception and cannot be passed down to children. There are two types of somatic mutations: driver mutations and passenger mutations. 

Passenger mutations are random damage that happened to occur as a result of the rapid growth. They don't contribute to tumour growth; in fact, [a 2017 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5639691/) found that a sufficiently increased number of passenger mutations can actually slow tumour growth. 

Driver mutations cause cells to become cancerous and multiply, actively accelerating tumour growth. Identifying individual driver mutations allows doctors to precisely target and treat specific cancers, but remains challenging due to the massive diversity of the many different mutated cells that make up cancer tumours. 

Because drugs are designed to target drivers and targeting passengers does nothing, separating driver and passenger mutations remains one of the central open problems in cancer genomics today. 

## pathways 
Genes change proteins, and proteins physically interact with each other in chains and loops, like a circuit. Damaging any one gene in a circuit can produce a signal that causes the cell to divide rapidly. 

This can be modelled as a graph analysis problem, where nodes represent genes, edges represent interactions between proteins, and node signals represent how often the gene is mutated across the cohort. The task, then, is to find regions of the graph where the signal concentrates. 

## first steps 
stay tuned...
