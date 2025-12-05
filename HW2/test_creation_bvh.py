import bvh
import numpy as np
import time

n=100
positions=np.random.rand(n,3)
radii=np.random.rand(n)*0.1+0.01

bvh1 = bvh.BVH(positions, radii)
bvh2 = bvh.BVH(positions, radii)

start=time.time()
bvh1.build()
stop=time.time()
print("Build time (SAH):",stop-start)
start=time.time()
bvh2.build_quick(splits=4)
stop=time.time()
print("Build time (quick):",stop-start)

print(bvh1.root.update_cost())
print(bvh2.root.update_cost())

positions+=np.random.rand(n,3)*0.08

bvh1.positions = positions
bvh2.positions = positions
bvh1.update()
bvh2.update()

print("After update:")
print(bvh1.root.update_cost())
print(bvh2.root.update_cost())

bvh1.optimize_w_rotation(bvh1.root, n)
bvh2.optimize_w_rotation(bvh2.root, n)

print("After optimization:")
print(bvh1.root.update_cost())
print(bvh2.root.update_cost())