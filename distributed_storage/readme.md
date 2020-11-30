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

## Sample runs:

```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 20
recovery threshold: 7
number of helpers: 13
number of nodes repaired: 7
size of base field: 23
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 91.0)
('minimum dimension of k-subspace: ', 91.0)
('Optimal desired dimension: ', 91)
```

```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 20
recovery threshold: 7
number of helpers: 12
number of nodes repaired: 7
size of base field: 23
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 84.0)
('minimum dimension of k-subspace: ', 84.0)
('Optimal desired dimension: ', 84)
```

```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 20
recovery threshold: 7
number of helpers: 11
number of nodes repaired: 7
size of base field: 23
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 77.0)
('minimum dimension of k-subspace: ', 77.0)
('Optimal desired dimension: ', 77)

```

```sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage"
....: )
number of nodes: 20
recovery threshold: 7
number of helpers: 10
number of nodes repaired: 7
size of base field: 23
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 70.0)
('minimum dimension of k-subspace: ', 70.0)
('Optimal desired dimension: ', 70)

```

```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage"
....: )
number of nodes: 20
recovery threshold: 7
number of helpers: 9
number of nodes repaired: 7
size of base field: 23
('initial evaluation points are linearly independent - ', False)
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 63.0)
('minimum dimension of k-subspace: ', 63.0)
('Optimal desired dimension: ', 63)

```

```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage"
....: )
number of nodes: 20
recovery threshold: 7
number of helpers: 8
number of nodes repaired: 7
size of base field: 23
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 56.0)
('minimum dimension of k-subspace: ', 56.0)
('Optimal desired dimension: ', 56)

```
```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage"
....: )
number of nodes: 20
recovery threshold: 7
number of helpers: 7
number of nodes repaired: 7
size of base field: 23
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 49.0)
('minimum dimension of k-subspace: ', 49.0)
('Optimal desired dimension: ', 49)

```

 ```
 sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 20
recovery threshold: 10
number of helpers: 12
number of nodes repaired in a single round: 3
size of base field: 23
('initial evaluation points are linearly independent - ', True)
('The chosen packets after the first repair round are linearly independent: ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 105.22)
('minimum dimension of k-subspace: ', 85.0)
('Optimal desired dimension: ', 85)
```
```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 19
recovery threshold: 13
number of helpers: 13
number of nodes repaired: 4
size of base field: 19
('initial evaluation points are linearly independent - ', False)
('initial evaluation points are linearly independent - ', True)
('The chosen packets after the first repair round are linearly independent: ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 110.0)
('minimum dimension of k-subspace: ', 110.0)
('Optimal desired dimension: ', 110)


```

```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 16
recovery threshold: 9
number of helpers: 12
number of nodes repaired in a single round: 3
size of base field: 17
('initial evaluation points are linearly independent - ', True)
('The chosen packets after the first repair round are linearly independent: ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 93.02)
('minimum dimension of k-subspace: ', 85.0)
('Optimal desired dimension: ', 81)
```
```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 20
recovery threshold: 10
number of helpers: 12
number of nodes repaired in a single round: 5
size of base field: 23
('initial evaluation points are linearly independent - ', True)
('The chosen packets after the first repair round are linearly independent: ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 118.68)
('minimum dimension of k-subspace: ', 113.0)
('Optimal desired dimension: ', 95)
```
```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 19
recovery threshold: 9
number of helpers: 11
number of nodes repaired in a single round: 4
size of base field: 19
('initial evaluation points are linearly independent - ', True)
('The chosen packets after the first repair round are linearly independent: ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 87.78)
('minimum dimension of k-subspace: ', 76.0)
('Optimal desired dimension: ', 76)
```

```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 10
recovery threshold: 4
number of helpers: 6
number of nodes repaired: 4
size of base field: 11
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 24.0)
('minimum dimension of k-subspace: ', 24.0)
('Optimal desired dimension: ', 24)

```
```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 12
recovery threshold: 8
number of helpers: 8
number of nodes repaired: 4
size of base field: 13
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 49.0)
('minimum dimension of k-subspace: ', 49.0)
('Optimal desired dimension: ', 48)

```
```
sage: load("~/Distributed-matrix-computations/MBR_functional_regenerating_code.sage")
number of nodes: 16
recovery threshold: 7
number of helpers: 9
number of nodes repaired: 7
size of base field: 17
('initial evaluation points are linearly independent - ', True)
Performing 30 repair rounds of random failures and randomly chosen helper nodes..
Randomly checking the dimension of various k-subspaces..
('average dimension of 50 k-subspaces: ', 63.0)
('minimum dimension of k-subspace: ', 63.0)
('Optimal desired dimension: ', 63)

```
