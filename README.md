# Nitish
Contains my publications on coded caching, distributed storage, and recent work on distributed computation with robustness against stragglers, byzantine servers, and honest but curious servers.
I plan to put implementations of my research on this repository, and other interesting readings on deep neural networks and other topics.

# List of codes in the repository:

## 1) SDMC_DFT.sage :

Author: Nitish Mital

Implementation of the paper "Secure Distributed Matrix Computation with Discrete Fourier Transform" - Nitish Mital, Cong Ling, Deniz Gunduz. 

Open access - https://arxiv.org/abs/2007.03972

We consider the problem of secure distributed matrix computation (SDMC), where a user can query a function of data matrices generated at distributed source nodes. We assume the availability of N honest but curious computation servers, which are connected to the sources, the user, and each other through orthogonal and reliable communication links. Our goal is to minimize the amount of data that must be transmitted from the sources to the servers, called the upload cost, while guaranteeing that no T colluding servers can learn any information about the source matrices, and the user cannot learn any information beyond the computation result.

This code simulates distributed source nodes, a computing cluster with N servers, and the user wanting to compute a particular function of the input matrices generated on the source nodes using the computing cluster. 

Class Cluster emulates a cluster of N computing servers. It has cluster level functions like inter-server communication.

Class Worker emulates a computing server. It has server level functions like matrix operations, encoding and decoding.

Class Source emulates the data sources. It has source node level functions like generation of data, encoding, and uploading the data to the computing cluster.

The fixed part of the code is the classes and functions, while the main section consists of a demonstration of how to use the classes to perform desired matrix computations and operations.

This code implements matrix multiplication, addition, conversion of shares. The matrix exponentiation and matrix inversion algorithms from the paper are not included in this code. Also, unlike the paper, we employ a naive discrete Fourier transform algorithm, that is, simple evaluation of a polynomial on the roots of unity, instead of applying the fast Fourier transform algorithm.

## MBR_functional_regenerating_code.sage :

Author - Nitish Mital

Implementation of my paper - N. Mital, K. Kralevska, C. Ling, D. Gunduz, "Practical Functional Regenerating Codes for Broadcast Repair of Multiple Nodes", ISIT Paris, 2019. 
 
Open access - https://arxiv.org/abs/1904.07273 
 
The code implementation resembles centralized repair more closely because all repair packets are collected in one matrix (Yx, Yy). Theoretically, centralized repair and broadcast repair are equivalent.

The problem considers a set of n storage nodes, in which a file is stored using a distributed storage code with redundancy so that it is sufficient for any k out of n nodes to be active for a user to retrieve its desired file from them. If r nodes have become out of service, the problem considers the minimization of the amount of communication, called repair bandwidth, required to populate with content the newcomer nodes, that replace the out of service nodes, from the transmitted packets by helper nodes, a subset of the surviving nodes.
The code -
a) Generates a random message, encodes and stores it in n nodes using the code proposed in the paper.
b) Contains a function "repair()" that takes as argument a list of helper node indices and newcomer node indices, and implements the novel repair procedure proposed in the paper.

This code works only if k and d are divisible by r. An update for other values will be posted later.


## MatDot.sage :

Author - Nitish Mital

Implementation of the paper "On the Optimal Recovery Threshold of Coded Matrix Multiplication" - Sanghamitra Dutta, Mohammad Fahim, Farzin Haddadpour, Haewon Jeong, Viveck Cadambe.

Open access - https://arxiv.org/abs/1801.10292

This problem considers the problem of multiplication of massive matrices, for applications such as large scale machine learning, on distributed servers, in which the variability of job completion times of different servers causes delays. Slow servers are referred to as stragglers. This approach uses coding theory to distribute parts of the computation to different servers with redundancy, so as to leverage the advantage of distributed computing as well as mitigate the delay caused by straggling servers. 

The model implemented in this code is slightly different from that in the paper, because we consider that the source matrices are generated at distributed source nodes, instead of at the user.

This code simulates distributed source nodes, a computing cluster with N computing servers, and the user wanting to compute the multiplication of the input matrices generated on the source nodes using the computing cluster. Redundancy is introduced using MatDot codes proposed in the paper, such that the result can be obtained by the user from any K out of N workers that return their computations first.

## polynomial_straggler.sage :

Author - Nitish Mital

Implementation of the paper "Polynomial Codes: an Optimal Design for High-Dimensional Coded Matrix Multiplication" - Qian Yu, Mohammad Ali Maddah-Ali, Salman Avestimehr.

Open access - https://arxiv.org/abs/1705.10464

Considers the same problem and system model as MatDot.sage above, but employs a different coding algorithm.

This code simulates source nodes, a computing cluster with N computing servers, and the user wanting to compute the multiplication of the input matrices generated on the source nodes using the computing cluster. Redundancy is introduced using polynomial codes proposed in the paper, such that the result can be obtained by the user from any K1K2 out of N workers that return their computations first.
