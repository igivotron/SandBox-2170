


# POISSON RECONSTRUCTION

**Authors:** Lionel Peduzzi, Myrine Msallem, Igor Grégoire


## Overview

This project's goal is the reconstruct a continuous, watertight surface from a discrete set of points. 
The Poisson Surface Reconstruction method provides
a global, smooth and robust solution based on solving a Poisson equation that integrates the normal field
of an oriented point cloud.

## Project Structure

```
.
├── data/          
│   ├── bunny.ply          # Sample point dataset
│   └── PLYgen.py          # Generates simple PLY file      
├── marchingCubes/         
│   ├── MC.py              # [NOT USED] Python implementation of marching cubes algorithm
│   └── triangleTable.py   # [NOT USED] Triangles table for the python implementation
│   ├── Mcc.c              # C implementation of marching cubes algorithm
│   └── Mcc.h              # Headers for the C implementation
├── output/                # Output folder for the surface
├── shared_lib/            # Python-C binding files
├── main.py                # [MAIN] File to execute
├── poisson.py             # Python file containing the vector field construction and the resolution of Poisson
└── vector_field.c         # C file containing the vector field construction
```

## Usage
### ⚠️ Important Note for Cross-Platform Usage

**<span style="color: red;">If you use different Operating Systems, you need to change line 13 in `poisson.py` and line 21 in `main.py`</span>**

### Binding the C-files
Bind the C-files with python to execute the code.


- For MacOS:
```gcc
gcc -shared -o shared_lib/vector_field.so -fPIC -O3 vector_field.c
gcc -shared -o shared_lib/Mcc.so -fPIC ./marchingCubes/Mcc.c -O3  
```
- For Linux:
```gcc
gcc -shared -o shared_lib/vector_field.so -fPIC -O3 vector_field.c
gcc -shared -o shared_lib/Mcc.so -fPIC ./marchingCubes/Mcc.c -O3  
```

- For Windows:
```gcc
gcc -shared -o shared_lib/vector_field.dll -fPIC -O3 vector_field.c
gcc -shared -o shared_lib/Mcc.so -fPIC ./marchingCubes/Mcc.c -O3  
```

### Running Surface reconstruction

```python
python main.py -i <input_file> -o <output_file> -N <int> -k <int> -skfmm <0 or 1>
```
- `-i`: Input file path
- `-o`: Output file path
- `-N`: Grid size (NxNxN)
- `-k`: Number of closest neighbours to take account for the spatial research knn-tree
- `-skfmm`: Boolean, utilisation of skfmm for the construction of the iso-surface. Might have negative impact on the execution time.

## Algorithm

### Computing the normals
1. Build a knn-graph with the points
2. Compute the PCA for each points using the k closest neighbours
3. Compute the normal by taking the eigenvector associated with the smallest eigenvalue
4. Weight the graph edges:  $W_{i,j} = n_{i} \cdot n_{j}$
5. Compute the minimum spanning tree
6. Orient the normals by propaging a normal chosen arbitrarily to his neighbours

### Solve the Poisson equation
![Go look wikipedia](explications/part1.png)
![Go look wikipedia](explications/part2.png)

### Marching Cubes




## Complexity

### Execution time repartition
![Computation times](plots/timing_plot.svg)

### Total execution time complexity
![Computation times](plots/complex.svg)

### Marching cubes complexity
![Computation times](plots/marching_cubes_complexity.svg)


## Task's done by members
- **Igor** :
    - Normals

- **Myrine** : 
    - Marching cubes

- **Lionel** : 
    - Vector field
