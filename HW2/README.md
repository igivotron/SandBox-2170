


# Spatial Data Structure - BVH

**Authors:** Lionel Peduzzi, Myrine Msallem, Igor Grégoire


## Overview

This project implements bounding volume hierarchies constructed top-down to a software simulating spheres collisions inside a box. The aim of the implementation is to improve the performance of geometric operation and GPU based ray-tracing.

## Project Structure

```
.
├── shaders/               # Shaders files (GLSL)
├── logs/                  # FPS logfiles
├── bvh_C.dll              # Python-C binding file (Windows)
├── bvh_C.so               # Python-C binding file (Linux/MacOS)
├── bvh.c                  # BVH - implementation
├── bvh.h                  # BVH - header
├── bvh.py                 # BVH - python implementation
├── homework.py            # Main Python interface
└── simulator.py           # Physic simulation
```

## Usage
### ⚠️ Important Note for Cross-Platform Usage

**<span style="color: red;">If you use different Operating Systems, you need to change line 15 in `homework.py`</span>**

### Binding the C-files
Bind the C-files with python to execute the code.


- For MacOS:
```gcc
gcc -shared -o bvh_C.so -fPIC bvh.c -lm -O3
```
- For Linux:
```gcc
gcc -shared -o bvh_C.so -fPIC bvh.c -lm -O3  
```

- For Windows:
```gcc
gcc -shared -o bvh_C.dll -fPIC bvh.c -lm -O3 
```

### Running Simulation

```python
python homework.py
```
Run the simulation. Use the interface to add/remove balls from the simulation and activate ray-tracing.


## Input/Output

- **Input:** N/A
- **Output:** FPS for each step: `logs/trianglefps_log.txt`

## Dependencies

See  `requirement.txt`
- python 3.x
- glfw 2.10.0
- imgui-bundle 1.92.4
- numpy 2.3.4
- rendercanvas 2.3.0
- wgpu 0.27.0

## Algorithm
We've implemented the algorithm throught a tree.
1. Tree construction:
    - The root node contains all the index.
    - We look for the best cut on each axis minimizing the heuristic.
    - We separate the points into two subnodes with respect to the best cut.
    - We continue iteratively until reaching all the leaves.
2. Tree update:
    - Every iterations, we update the hit boxes with wrt the new points positions iteratively from leaves to root.
    - Every 6 iterations, we do a rotation on the nodes to minimize the heurisitc
3. Itersection finding:
    - To find intersection between two boxes, we compare their extremities.
    - We iteratively go from parent to children following the intersections between boxes until reaching leaves.
4. RayTracing:
    - We flatten the BVH into a 1D vector containing for each node:

    ```[left, right, item, padding, bbox_min(x,y,z), padding, bbox_max(x,y,z), n_items]```

    The vector contains padding to facilitate its utilization by the GPU
    - We send the vector to the GPU
    - In the GPU, we compute the intersection by iteratively go from parent to children following the ray-box intersection until reaching leaves.

## Performance
All the performance have been retrieved on 

    - MacOS Sequoi
    - RAM: 8go
    - CPU: M1

- NO-RTX
    - 20 balls: 200fps
    - 100 balls: 120fps
    - 250 balls: 25fps

- RTX
    - 1 balls: 1fps

    This results is not representative, the computer not having a GPU.


## Task's done by members
- **Igor** :
    - Static bvh implementation on python
    - Translation of the code in C
    - Flatten bvh

- **Myrine** : 
    - Node rotation in python
    - Finding intersection
    - LBVH : Morton code in python and C


- **Lionel** :
    - Raytracing
    - Intersection
    - C implementation



