# MSB1014-Network_Biology-Project
This repository is the final product of the project, requested by the course MSB1014 (Network Biology).

[![GitHub License](https://img.shields.io/github/license/Rrtk2/MSB1014-Network_Biology-Project)](https://github.com/Rrtk2/MSB1014-Network_Biology-Project/blob/master/LICENSE.md) ![](https://img.shields.io/badge/Status-Setting_up-orange) [![GitHub Watches](https://img.shields.io/github/watchers/Rrtk2/MSB1014-Network_Biology-Project.svg?style=social&label=Watch&maxAge=2592000)](https://github.com/Rrtk2/MSB1014-Network_Biology-Project/watchers) 


#### What is this project about
With the era of high-throughput data, the amount of gene expression, microRNA expression, and methylation data has increased tremendously while the understanding of these regulatory mechanisms is lagging behind. (Epi)genetic dysregulation can lead to altered phenotypes such as cancer, therefore, understanding the inference of these regulatory mechanisms is essential to treat or cure altered phenotypes. Little is known on how these regulating layers work together in detail. Methods to infer regulatory mechanisms from experimental and modelled datasets have been studied greatly, and many methods have been developed. However, each method has its own strengths and weaknesses. By combining these methods in an ensemble fashion, the strengths can be combined which creates a better initial prediction of the inference in the network. To further improve the prediction, statistically significant overrepresented topological patterns (network motifs) are used as an additional filter to extract biologically-derived significant patterns found in the inferred networks. 

This project aims to estimate the improvement of network inference of regulatory networks by using network motifs.
#### Usage

#### Expected output

#### Project structure
##### Where does data come from?


##### How is data shared, in what format, with what protocols?

##### Workflow
Data is obtained from a dedicated network inference challenge, the HPN-DREAM breast cancer network inference challenge (https://www.synapse.org/HPN_DREAM_Network_Challenge). This dataset was generated using in sillico simulations of a model described in Chen et al. This means that the inferred networks can always be compared to the correct, original model. The simulated model data will be analyzed in R using the following network inference techniques; Bayesian networks, correlation, regression and mutual information. Each inferred network will be analyzed on network motifs using the ‘iGraph’ package in R. The resultant motifs possibly indicate which part of the inferred network is more likely to be coherent to the true model network, this information will be used to increase network inference prediction. By combining the resulted motif-filtered inference networks using an ensemble approach, predictive power can be increased. By comparing the ensemble network with the original network, the overlap and mismatch can be determined which indicates the performance of prediction. 

#### Contact
ra.reijnders@student.maastrichtuniversity.nl


#### License and contributing guidelines
[License](/LICENSE.md) 

[Contributing guidelines](/CONTRIBUTING.md) 


#### Who is involved, and what are their roles.
RRtK2 (owner and contributor)


#### Status of project
Setting up. Large edits; adding missing files; trying concepts.


#### Copyright and authors
All code and documents in the MSB1014-Network_Biology-Project folder was created by [these author(s)](/AUTHORS.md).
