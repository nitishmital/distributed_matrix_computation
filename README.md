# Nitish
Contains my publications on coded caching, distributed storage using linearized polynomials, and recent work on distributed computation with robustness against stragglers, byzantine servers, and honest but curious servers.
I plan to put implementations of my research on this repository, and other interesting readings on deep neural networks and other topics.

## 1) MatDot.sage :
Implementation of the paper "On the Optimal Recovery Threshold of Coded Matrix Multiplication" - Sanghamitra Dutta, Mohammad Fahim, Farzin Haddadpour, Haewon Jeong, Viveck Cadambe.

This problem consider the problem of multiplication of massive matrices, for applications such as large scale machine learning, on distributed servers, in which the variability of job completion times of different servers causes delays. Slow servers are referred to as stragglers. This approach uses coding theory to distribute parts of the computation to different servers with redundancy, so as to leverage the advantage of distributed computing as well as mitigate the delay caused by straggling servers. 

The model implemented in this code is slightly different from that in the paper, because we consider that the source matrices are generated at distributed source nodes, instead of at the user.

This code simulates source nodes, computing cluster nodes with N workers, and the user wanting to compute the multiplication of the input matrices generated on the source nodes using the computing cluster. Redundancy is introduced using MatDot codes proposed in the paper, such that the result can be obtained by the user from any K out of N workers that return their computations first.

## polydot.sage :

Implementation of the paper "Polynomial Codes: an Optimal Design for High-Dimensional Coded Matrix Multiplication" - Qian Yu, Mohammad Ali Maddah-Ali, Salman Avestimehr.

Considers the same problem and system model as MatDot.sage above, but employs a different coding algorithm.

This code simulates source nodes, computing cluster nodes with N workers, and the user wanting to compute the multiplication of the input matrices generated on the source nodes using the computing cluster. Redundancy is introduced using polynomial codes proposed in the paper, such that the result can be obtained by the user from any K1K2 out of N workers that return their computations first.

## MBR_functional_regenerating_code.sage :

 Implementation of my paper "Practical Functional Regenerating Codes for Broadcast Repair of Multiple Nodes", ISIT Paris, 2019. The code implementation resembles centralized repair more closely because all repair packets are collected in one matrix (Yx, Yy). Theoretically, centralized repair and broadcast repair are equivalent.

The problem considers a set of n storage nodes, in which a file is stored using a distributed storage code with redundancy so that it is sufficient for any k out of n nodes to be active for a user to retrieve its desired file from them. If r nodes have become out of service, the problem considers the minimization of the amount of communication, called repair bandwidth, required to populate with content the newcomer nodes, that replace the out of service nodes, from the transmitted packets by helper nodes, a subset of the surviving nodes.
The code -
a) Generates a random message, encodes and stores it in n nodes using the code proposed in the paper.
b) Contains a function repair() that takes as argument a list of helper node indices and newcomer node indices, and implements the novel repair procedure proposed in the paper.

This code works only if k and d are divisible by r. An update for other values will be posted later.


