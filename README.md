Contains codes for my publications on distributed storage and private distributed computation, and related papers on distributed computation with robustness against stragglers, byzantine servers, and honest but curious servers.

# List of files in the repository:

## 1) SDMC_DFT.sage :

Implementation of the paper "Secure Distributed Matrix Computation with Discrete Fourier Transform" - Nitish Mital, Cong Ling, Deniz Gunduz. 

## Citation:
```bash
@article{DBLP:journals/corr/abs-2007-03972,
  author    = {Nitish Mital and
               Cong Ling and
               Deniz G{\"{u}}nd{\"{u}}z},
  title     = {Secure Distributed Matrix Computation with Discrete Fourier Transform},
  journal   = {CoRR},
  volume    = {abs/2007.03972},
  year      = {2020},
  url       = {https://arxiv.org/abs/2007.03972},
  eprinttype = {arXiv},
  eprint    = {2007.03972},
  timestamp = {Mon, 20 Jul 2020 14:20:39 +0200},
  biburl    = {https://dblp.org/rec/journals/corr/abs-2007-03972.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```
Open access - https://arxiv.org/abs/2007.03972

We consider the problem of secure distributed matrix computation (SDMC), where a user can query a function of data matrices generated at distributed source nodes. We assume the availability of N honest but curious computation servers, which are connected to the sources, the user, and each other through orthogonal and reliable communication links. Our goal is to minimize the amount of data that must be transmitted from the sources to the servers, called the upload cost, while guaranteeing that no T colluding servers can learn any information about the source matrices, and the user cannot learn any information beyond the computation result.

This code simulates distributed source nodes, a computing cluster with N servers, and the user wanting to compute a particular function of the input matrices generated on the source nodes using the computing cluster. 

Class Cluster emulates a cluster of N computing servers. It has cluster level functions like inter-server communication.

Class Worker emulates a computing server. It has server level functions like matrix operations, encoding and decoding.

Class Source emulates the data sources. It has source node level functions like generation of data, encoding, and uploading the data to the computing cluster.

The fixed part of the code is the classes and functions, while the main section consists of a demonstration of how to use the classes to perform desired matrix computations and operations.

This code implements matrix multiplication, addition, conversion of shares. The matrix exponentiation and matrix inversion algorithms from the paper are not included in this code. Also, unlike the paper, we employ a naive discrete Fourier transform algorithm, that is, simple evaluation of a polynomial on the roots of unity, instead of applying the fast Fourier transform algorithm.


## MatDot.sage :

Implementation of the paper "On the Optimal Recovery Threshold of Coded Matrix Multiplication" - Sanghamitra Dutta, Mohammad Fahim, Farzin Haddadpour, Haewon Jeong, Viveck Cadambe.

Open access - https://arxiv.org/abs/1801.10292

This problem considers the problem of multiplication of massive matrices, for applications such as large scale machine learning, on distributed servers, in which the variability of job completion times of different servers causes delays. Slow servers are referred to as stragglers. This approach uses coding theory to distribute parts of the computation to different servers with redundancy, so as to leverage the advantage of distributed computing as well as mitigate the delay caused by straggling servers. 

The model implemented in this code is slightly different from that in the paper, because we consider that the source matrices are generated at distributed source nodes, instead of at the user.

This code simulates distributed source nodes, a computing cluster with N computing servers, and the user wanting to compute the multiplication of the input matrices generated on the source nodes using the computing cluster. Redundancy is introduced using MatDot codes proposed in the paper, such that the result can be obtained by the user from any K out of N workers that return their computations first.

## polynomial_straggler.sage :

Implementation of the paper "Polynomial Codes: an Optimal Design for High-Dimensional Coded Matrix Multiplication" - Qian Yu, Mohammad Ali Maddah-Ali, Salman Avestimehr.

Open access - https://arxiv.org/abs/1705.10464

Considers the same problem and system model as MatDot.sage above, but employs a different coding algorithm.

This code simulates source nodes, a computing cluster with N computing servers, and the user wanting to compute the multiplication of the input matrices generated on the source nodes using the computing cluster. Redundancy is introduced using polynomial codes proposed in the paper, such that the result can be obtained by the user from any K1K2 out of N workers that return their computations first.
