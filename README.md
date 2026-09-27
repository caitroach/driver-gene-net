# are known ovarian cancer driver genes more densely interconnected than expected by chance?
(Yes! Even after controlling for network degree!!)
Among 49 ovarian epithelial tumour driver genes present in the PPI network, I observed 49 driver-driver interactions. Degree-matched random gene sets contained 12.6 ± 3.4 interactions on average across 1,000 permutations (empirical p ≤ 0.001).
<img width="2779" height="2376" alt="image" src="https://github.com/user-attachments/assets/baeaee8b-9bfb-46c0-a69f-fd2c41807ac1" />
(Each circle is a driver gene for ovarian cancer; a line means the two proteins physically interact.)
Larger circles have more interaction partners overall. Most drivers (33 of 49) form one connected group, suggesting they act on shared machinery rather than independently. The 16 below the line have no direct link to any other driver. 
## background
Cancer is the result of cells growing uncontrollably because their DNA has been damaged. When tumour DNA is sequenced and compared to pt's healthy tissue, researchers find somatic mutations, which are changes in human DNA that occur after conception and cannot be passed down to children (think of cancer caused by smoking, sunburns, etc.). There are two types of somatic mutations: driver mutations and passenger mutations. 

Passenger mutations are random damage that happened to occur as a result of the rapid growth. They don't contribute to tumour growth; in fact, [a 2017 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5639691/) found that a sufficiently increased number of passenger mutations can actually slow it. 

Driver mutations cause cells to become cancerous and multiply, actively accelerating tumour growth. Identifying individual driver mutations allows doctors to precisely target and treat specific cancers, but remains challenging due to the massive diversity of the many different mutated cells that make up cancer tumours. 

Because drugs are designed to target drivers and targeting passengers does nothing, separating driver and passenger mutation data remains one of the central open problems in cancer genomics today. 

Genes code for proteins, and proteins interact with each other in chains and loops, like a circuit: A signal arrives at the cell surface, passes along through a series of proteins, and ends with the cell deciding to divide. 

<img width="1024" height="413" alt="image" src="https://github.com/user-attachments/assets/77b1d5e2-7bfe-40b7-86a4-d62e521d98cc" />

*Image via [The Baker Lab](https://www.bakerlab.org/2020/04/02/de-novo-design-protein-logic-gates/)*

Cancer doesn't need to destroy a specific gene. It only needs to break the pathway. Damaging any of several genes in one circuit produces the same resultant signal.

This can be modeled as a graph analysis problem, where nodes represent genes, edges represent interactions between proteins, and node signals represent how often the gene is mutated across the cohort. 

## methods
Every network-based method in cancer genomics relies on the assumption that driver genes sit close together in the protein interaction network. The edges in an interaction network are records of experiments we chose to run on specifically selected cancer genes, so drivers could potentially cluster in the network partly due to bias! Methods like [HotNet2](https://github.com/raphael-group/hotnet2) spread mutation signals across the graph because they assume this is true. 

This project asks one simple question: do the drivers for a given cancer work together, or does each one break something independently? To answer it, I worked in three steps.

`build_network.py` starts by downloading human protein interaction data from STRING, keeps only the confident physical interactions, converts protein IDs to gene names, and saves the biggest connected chunk as a tsv.

`driver_clustering.py` marks the 49 drivers on the map and counts how many links there are between drivers. Highly connected proteins are more likely to have interactions with one another simply because they have many interaction partners. To distinguish genuine clustering from this effect, this script makes 1,000 fake driver sets and counts their links too. There were two types of fake sets: uniform (totally random) and degree-matched (random but just as popular as real drivers).

Finally, `draw_network.py` uses matplotlib and networkx to draw the figure.

## results
<img width="2704" height="958" alt="image" src="https://github.com/user-attachments/assets/24dd66a5-c91f-4ab1-bf2c-96ef90df3c3e" />
The observed driver set contained 49 driver-driver edges, compared with 12.6 ± 3.4 in degree-matched random sets (about 3.9x more; 1,000 permutations; p ≤ 0.001). The largest connected group of drivers was also far bigger than chance: 33 of 49, vs 5.7 ± 2.0 in degree-matched sets.

The answer: driver genes cluster far more than chance allows. Most OVT drivers form a single connected group, and the effect holds up even after correcting for the fact that drivers tend to have more interaction partners than average. 

## data 
human protein-protein interaction network data via the [Swiss Institute of Bioinformatics](https://string-db.org/)
ovarian cancer driver data via [IntOGen](https://www.intogen.org/search?cancer=OVT)
