# Shape Metrics
Point Cloud Utils has functions to compute a number of commonly used metrics between point clouds.


## Chamfer Distance
The Chamfer distance between two point clouds $P_1 = \{x_i \in \mathbb{R}^3\}_{i=1}^n$ and $P_2 = \{x_j \in \mathbb{R}^3\}_{j=1}^m$ is defined as the average distance between pairs of nearest neighbors between $P_1$ and $P_2$ *i.e.*
$$
\text{chamfer}(P_1, P_2) = \frac{1}{2n} \sum_{i=1}^n \|x_i - \text{NN}(x_i, P_2)\| + \frac{1}{2m} \sum_{j=1}^n \|x_j - \text{NN}(x_j, P_1)\|
$$
and $\text{NN}(x, P) = \text{argmin}_{x' \in P} \|x - x'\|$ is the nearest neighbor function.

The following code computes the Chamfer distance between two point clouds:
```python
import point_cloud_utils as pcu

# p1 is an (n, 3)-shaped numpy array containing one point per row
p1 = pcu.load_mesh_v("point_cloud_1.ply")

# p2 is an (m, 3)-shaped numpy array containing one point per row
p2 = pcu.load_mesh_v("point_cloud_2.ply")

# Compute the chamfer distance between p1 and p2
cd = pcu.chamfer_distance(p1, p2)
```

## Hausdorff distance
The Hausdorff distance between two point clouds $P_1 = \{x_i \in \mathbb{R}^3\}_{i=1}^n$ and $P_2 = \{x_j \in \mathbb{R}^3\}_{j=1}^m$ is defined as the maxmimum distance between any pair of nearest neighbors between $P_1$ and $P_2$ *i.e.*
$$
\text{hausdorff}(P_1, P_2) = \frac{1}{2} \max_{x \in P_1} \|x - \text{NN}(x, P_2)\| + \frac{1}{2} \max_{x' \in P_2}  \|x' - \text{NN}(x', P_1)\|
$$
and $\text{NN}(x, P) = \text{argmin}_{x' \in P} \|x - x'\|$ is the nearest neighbor function.

The following code computes the Hausdorff distance between two point clouds:
```python
import point_cloud_utils as pcu

# p1 is an (n, 3)-shaped numpy array containing one point per row
p1 = pcu.load_mesh_v("point_cloud_1.ply")

# p2 is an (m, 3)-shaped numpy array containing one point per row
p2 = pcu.load_mesh_v("point_cloud_2.ply")

# Compute the chamfer distance between p1 and p2
hd = pcu.hausdorff_distance(p1, p2)
```

### One sided Hausdorff distance
In some applications, one only needs the *one-sided Hausdorff distance* between $P_1$ and $P_2$, *i.e.*

$$
\text{hausdorff}_{P_1 \rightarrow P_2}(P_1, P_2) = \max_{x \in P_1} \|x - \text{NN}(x, P_2)\|
$$

The following code computes the one-sided Hausdorff distance between two point clouds:
```python
import point_cloud_utils as pcu

# p1 is an (n, 3)-shaped numpy array containing one point per row
p1 = pcu.load_mesh_v("point_cloud_1.ply")

# p2 is an (m, 3)-shaped numpy array containing one point per row
p2 = pcu.load_mesh_v("point_cloud_2.ply")

# Compute the chamfer distance between p1 and p2
hd_p1_to_p2 = pcu.one_sided_hausdorff_distance(p1, p2)
```

!!! note
    To get the $P_2 \rightarrow P_1$ Hausdorff distance, just swap the arguments to `pcu.one_sided_hausdorff_distance`


## Earth-Mover's (Sinkhorn) distance
The [Earth Mover's distance](https://en.wikipedia.org/wiki/Earth_mover%27s_distance) between two point clouds $P = \{p_i \in \mathbb{R}^3\}_{i=1}^n$ and $Q = \{q_j \in \mathbb{R}^3\}_{j=1}^m$ is computed as the average distance between pairs of points according to an optimal correspondence $\pi \in \Pi(P, Q)$, where $\Pi(P, Q)$ is the set of $n \times m$ matrices where the rows and columns sum to one. The assignment $\pi$ is thus a matrix where $\Pi_{i,j}$ is a number between $0$ and $1$ denoting how much point $p_i$ and $q_j$ correspond. We can write the EMD formally as:
$$
\text{EMD}(P, Q) = \min_{\pi \in \Pi(P, Q)} \sum_{i = 1}^n \sum_{j = 1}^m \pi_{i,j} \|p_i - q_j\| % \langle \pi, D \rangle \qquad D_{ij} = \|p_i - q_j\|
$$

Point Cloud Utils implements the sinkhorn algorithm for computing the (approximate) Earth Mover's Distance. To compute the EMD, run:
```python
import point_cloud_utils as pcu

# p1 is an (n, 3)-shaped numpy array containing one point per row
p1 = pcu.load_mesh_v("point_cloud_1.ply")

# p2 is an (m, 3)-shaped numpy array containing one point per row
p2 = pcu.load_mesh_v("point_cloud_2.ply")

# Compute the chamfer distance between p1 and p2
emd, pi = pcu.earth_movers_distance(p1, p2)
```
