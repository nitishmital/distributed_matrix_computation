## Experiments
Verify the preservation of the reconstruction property, i.e., the dimension of subspace in `k` nodes by setting:
1. ```(n,k,d,r,q) = (19,12,12,4,19)``` (q is the base field size.. must be greater than n-1).
2. ``` (n,k,d,r,q) = (19,8,12,4,19) ```
3. ``` (n,k,d,r,q) = (24,8,16,8,29) ```
4. ``` (n,k,d,r,q) = (13,8,8,4,13) ```
5. ``` (n,k,d,r,q) = (14,9,9,3,17) ```
6. ``` (n,k,d,r,q) = (16,9,12,3,17) ```

## Installation of SageMath:
1. ``` pip install sagemath ```
2. Run sage in terminal with the command - ```> sage ```
3. ``` sage: load("~/Distributed-matrix-computations/MBR.sage") ```. Put appropriate path to the file instead of this path.  

## Sample run:
```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 16
recovery threshold: 9
number of helpers: 12
number of nodes repaired: 3
size of base field: 17
('initial evaluation points are linearly independent - ', True)
('The chosen packets are linearly independent: ', True)
('dimension of subspace: ', 86)
('Number of subpackets: ', 81)
 
