


# Delaunay Triangulation

**Authors:** Lionel Peduzzi, Myrine Msallem, Igor Grégoire


## Overview

This project implements 2D Delaunay triangulation using a robust Bowyer Watson algorithm with optimized point ordering via Hilbert Curve combining C implementations with Python bindings for performance and ease to use. Using Delaunay triangulation, we are also capable of implementing its dual graph, the Voronoi Cells Algorithm.

## Project Structure

```
.
├── header/                 # C header files
├── input_file/            # Sample points cloud files for triangulation
├── output_file/           # Generated triangulation results
│   ├── triangles.dat      # Output triangulation data
│   └── [plots]            # Visualizations from plot.py
├── shared_lib/            # Python-C binding files
├── BowyerWatson.c         # Bowyer Watson triangulation algorithm
├── HalfEdge.c             # Half-edge data structure routines
├── Hilbert.c              # Hilbert curve point sorting
├── predicates.c           # Routines for Fast Robust Geometric Predicates
├── del2d.py               # Main Python interface for triangulation
├── plot.py                # Visualization and comparison with scipy
├── gents.py               # Random point cloud generator
└── Interactive.py         # Interactive Delaunay Triangulation 
```

## Usage
### ⚠️ Important Note for Cross-Platform Usage

**<span style="color: red;">If you use different Operating Systems, you need to change line 23 in `del2d.py` and line 154 in `voronoi.py`</span>**

### Binding the C-files
Bind the C-files with python to execute the code.


- For MacOS:
```gcc
gcc -shared -o shared_lib/BowyerWatsonMACOS.so -fPIC BowyerWatson.c Hilbert.c HalfEdge.c predicates.c   
```
- For Linux:
```gcc
gcc -shared -o shared_lib/BowyerWatson.dll -fPIC BowyerWatson.c Hilbert.c HalfEdge.c predicates.c   
```

- For Windows:
```gcc
gcc -shared -o shared_lib/BowyerWatson.so -fPIC BowyerWatson.c Hilbert.c HalfEdge.c predicates.c   
```

### Running Triangulation

```python
python del2d.py -i <input_file> -o <output_file> -p <0 or 1>
```
Run the triangulation from the points in the input file and save the triangles inside the output file. If -p 1 then plot the result and save the figure inside output_file/

### Visualizing Results

```python
python plot.py
```

The visualization will display your triangulation and compare it with scipy's Delaunay implementation.

### Interactive Triangulation

```python
python voronoi.py
```
Interactive visualization of the Voronoi Algorithm. 

`r` : Switch between Voronoï and Delaunay visualisation

`e` : Relax the graph using Lloyd

- <span style="color: green;"> Green lines</span> are finite Voronoï boundaries
- <span style="color: red;"> Red doted lines</span> are infinite Voronoï boundaries
- <span style="color: lightblue;"> blue lines</span> are the triangles edges
- <span style="color: white;"> Black dots</span> are vertices
- <span style="color: red;"> Red dots</span> are the circumcenter of the triangles




## Input/Output

- **Input:** Point cloud files in `input_file/`
- **Output:** Triangulation data written to `output_file/triangles.dat`

## Dependencies

- Python 3.x
- NumPy
- SciPy (for comparison)
- Matplotlib (for plotting)
- C compiler (GCC or compatible)

## Algorithm

The implementation uses the Bowyer-Watson incremental algorithm:

1. Points are sorted along a Hilbert curve for spatial locality
2. Super-triangles containing all points are created
3. Points are inserted incrementally, maintaining the Delaunay property and using robust routines
4. The super-triangles are removed

## To infinite and beyond
We've implemented the following improvements:

1. Removing infinite points
2. Robustness
3. An implementation in O(n log n)
4. Go Interactive

## Task's done by members
- **Igor** :
    - Implementation of the Hilbert curve
    - Robustness
    - Bowyer Watson algorithm
    - Removing infinite points
    - Lloyd relaxation 

- **Myrine** : 
    - Removing infinite points
    - Voronoï algorithm
    - Interactiveness of the plot

- **Lionel** : 
    - Bowyer Watson algorithm in $\mathcal{O}(n\log(n))$
    - Binding C-files with python files

