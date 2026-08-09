# Do known cancer driver genes sit closer together in the human protein interaction network than you'd expect by chance? Does that survive controlling for degree?
## background
Cancer is the result of cells growing uncontrollably because their DNA has been damaged. When tumour DNA is sequenced and compared to pt's healthy tissue, researchers find somatic mutations, which are changes in human DNA that occur after conception and cannot be passed down to children. There are two types of somatic mutations: driver mutations and passenger mutations. 

Passenger mutations are random damage that happened to occur as a result of the rapid growth. They don't contribute to tumour growth; in fact, [a 2017 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5639691/) found that a sufficiently increased number of passenger mutations can actually slow tumour growth. 

Driver mutations cause cells to become cancerous and multiply, actively accelerating tumour growth. Identifying individual driver mutations allows doctors to precisely target and treat specific cancers, but remains challenging due to the massive diversity of the many different mutated cells that make up cancer tumours. 

Because drugs are designed to target drivers and targeting passengers does nothing, separating driver and passenger mutation data remains one of the central open problems in cancer genomics today. 

Genes code for proteins, and proteins interact with each other in chains and loops, like a circuit: A signal arrives at the cell surface, passes along through a series of proteins, and ends with the cell deciding to divide. 

<img width="1024" height="413" alt="image" src="https://github.com/user-attachments/assets/77b1d5e2-7bfe-40b7-86a4-d62e521d98cc" />

*Image via [The Baker Lab](https://www.bakerlab.org/2020/04/02/de-novo-design-protein-logic-gates/)*

Cancer doesn't need to destroy a specific gene. It only needs to break the pathway. Damaging any of several genes in one circuit produces the same resultant signal.

This can be modeled as a graph analysis problem, where nodes represent genes, edges represent interactions between proteins, and node signals represent how often the gene is mutated across the cohort. 

## methods
Every network-based method in cancer genomics relies on the assumption that driver genes sit close together in the protein interaction network, because cancer works by breaking pathways rather than individual genes. The edges in an interaction network are records of experiments we chose to run on specifically selected cancer genes, so drivers may cluster in the network partly due to bias.

Methods like [HotNet2](https://github.com/raphael-group/hotnet2) spread mutation signals across the graph because they assume this is true. 

This project asks one simple question: do the drivers for a given cancer work together, or does each one break something independently? To answer it, I used a map of which proteins physically contact each other inside human cells, marked the 49 known drivers of ovarian epithelial tumours on that map, and counted how often two drivers were directly linked. 

## results
<img width="1549" height="607" alt="image" src="https://github.com/user-attachments/assets/a897281a-7bca-484f-9d36-91b916377d06" />
*Terminal results*
<img width="2779" height="2376" alt="image" src="https://github.com/user-attachments/assets/baeaee8b-9bfb-46c0-a69f-fd2c41807ac1" />
*Each circle is a driver gene for ovarian cancer; a line means the two proteins physically interact.*
Larger circles have more interaction partners overall. Most drivers (33 of 49) form one connected group, suggesting they act on shared machinery rather than independently. The 16 below the line have no direct link to any other driver. 

The answer: driver genes cluster far more than chance allows. Most OVT drivers form a single connected group, and the effect holds up even after correcting for the fact that drivers tend to be well-connected due to the research emphasis on them. 

## bibliography 
human protein-protein interaction network data via the [Swiss Institute of Bioinformatics] (https://string-db.org/)  
ovarian cancer driver data via [IntOGen](https://www.intogen.org/search?cancer=OVT)
