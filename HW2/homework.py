"""LMECA2170 — Homework 2 — Spatial Datastructures."""

from simulator import Simulator, Contact, copy_buffer
import wgpu
import numpy as np
import bvh
import time
import ctypes
import os

i=0

# Machin qui fait le lien entre python et C
# lib = ctypes.CDLL(os.path.abspath("bvh_C.so"))

lib = ctypes.CDLL(os.path.abspath("bvh_C.dll"))

class BVHNode(ctypes.Structure):
    _fields_ = [
        ("left", ctypes.c_int),                         # int
        ("right", ctypes.c_int),                        # int
        ("parent", ctypes.c_int),                       # int
        ("bbox", ctypes.POINTER(ctypes.c_double)),      # double*
        ("items", ctypes.POINTER(ctypes.c_int)),        # int*
        ("n_items", ctypes.c_int),                      # int
        ("index", ctypes.c_int),                        # int
        ("state", ctypes.c_int),                        # int
    ]

class bvh_C(ctypes.Structure):
    _fields_ = [
        ("nodes", ctypes.POINTER(BVHNode)),     # BVHNode*
        ("n_nodes", ctypes.c_int),              # int
        ("positions", ctypes.POINTER(ctypes.c_double)),  # double*
        ("radii", ctypes.POINTER(ctypes.c_double)),      # double*
        ("NperLeaf", ctypes.c_int),             # int
        ("root", ctypes.c_int),                 # int
    ]


lib.create_bvh.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int]
lib.create_bvh.restype = ctypes.POINTER(bvh_C)
lib.update.argtypes = [ctypes.POINTER(bvh_C), ctypes.POINTER(BVHNode)]
lib.update.restype = None
lib.free_bvh.argtypes = [ctypes.POINTER(bvh_C)]
lib.free_bvh.restype = None
lib.update_positions.argtypes = [ctypes.POINTER(bvh_C), ctypes.POINTER(ctypes.c_double)]
lib.update_positions.restype = None
lib.build_bvh.argtypes = [ctypes.POINTER(bvh_C)]
lib.build_bvh.restype = None
lib.is_leaf.argtypes = [ctypes.POINTER(BVHNode), ctypes.c_int]
lib.is_leaf.restype = ctypes.c_bool
lib.find_pot_inter.argtypes = [ctypes.POINTER(bvh_C), ctypes.POINTER(ctypes.c_int)]
lib.find_pot_inter.restype = ctypes.c_int



class Homework:
    def __init__(self):
        """Caled at the start of the simulation.

        You can use this to set the initial state of the datastructures that you will need.

        The simulator object contains the following attributes:
            - radii: array of shape (N,)
            - positions: array of shape (N, 3)
            - rotations: array of shape (N, 4), storing rotations as quaternions
              (each row contains the scalar first, followed by the
               x, y, and z components).
            - angular_velocities: array of shape (N, 3)
            - velocities: array of shape (N, 3)

        where N is the number of objects being simulated.

        You may use all of these attributes throughout the homework.
        """
        self.bvh_tree = None
        self.data_structure_staging = None
        self.data_structure_gpu = None

    def on_size_changed(self, simulator):
        """Called every time objects are added to or removed from the simulation.

        This can be used to update the data structures you use, if they change
        significantly when new objects are added to the simulation.
        """
        if len(simulator.positions)==0:
            self.bvh_tree = None
            return
        
        if self.bvh_tree is not None:
            lib.free_bvh(self.bvh_tree)
            self.bvh_tree = None
        
        # flatten positions and radii for C compatibility
        pos = simulator.positions.astype(np.float64).flatten()
        radii = simulator.radii.astype(np.float64).flatten()
        pos_c = pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        radii_c = radii.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        N = len(simulator.positions)
        bvh_tree_c = lib.create_bvh(pos_c, radii_c, N, 1)
        lib.build_bvh(bvh_tree_c)
        self.bvh_tree = bvh_tree_c

        # bvh_tree = bvh.BVH(simulator.positions, simulator.radii)
        # bvh_tree.build()
        # # bvh.print_bvh(bvh_tree.root)
        # self.data_structure_staging = bvh_tree
        global i
        i=0
        return

    def find_intersections(self, simulator) -> list[Contact]:
        """Return all contacts using the BVH."""
        global i
        i+=1
        bvh_tree = self.bvh_tree
        if bvh_tree is None:
            return []
        # bvh_tree.positions = simulator.positions
        pos = simulator.positions.astype(np.float64).flatten()
        pos_c = pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lib.update_positions(bvh_tree, pos_c)

        if i==100:
            #every 100 frames, rebuild the BVH
            lib.free_bvh(bvh_tree)
            radii = simulator.radii.astype(np.float64).flatten()
            radii_c = radii.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            bvh_tree = lib.create_bvh(pos_c, radii_c, len(simulator.positions), 1)
            self.bvh_tree = bvh_tree
            lib.build_bvh(bvh_tree)
            # bvh_tree.build()
            i=0
        else:
            #otherwise, just update the positions
            # bvh_tree.update()
            current = ctypes.pointer(bvh_tree.contents.nodes[bvh_tree.contents.root])
            lib.update(bvh_tree, current)
        
        potential_intersections = np.zeros((len(simulator.positions)**2), dtype=np.int32)
        potential_intersections_c = potential_intersections.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        n=lib.find_pot_inter(bvh_tree, potential_intersections_c)
        contacts = []
        for i in range(n//2):
            a = potential_intersections[2*i]
            b = potential_intersections[2*i+1]
            contact = simulator.intersect(a, b)
            if contact:
                contacts.append(contact)
        return contacts

    def gpu_bind_group_layouts(self, device: wgpu.GPUDevice) -> list[wgpu.GPUBindGroupLayout]:
        """For the GPU part of the project. This allows you to describe the memory resources you'll access from your GPU implementation.

        The first binding group returned will have correspond to descriptor set
        2 in your GPU code, the second to descriptor set 3, and so on.
        """
        return []

        # This example code lets you access one buffer (at descriptor set 2, binding 0) in your GPU code
        self.example_layout = device.create_bind_group_layout(
            label="Example Datastructure",
            entries=[wgpu.BindGroupLayoutEntry(binding=0, visibility=wgpu.ShaderStage.COMPUTE,
                                               buffer=wgpu.BufferBindingLayout(
                                                   type=wgpu.BufferBindingType.storage))]
        )

        return [self.example_layout]

    def update_gpu_data(self, simulator: Simulator, device: wgpu.GPUDevice, encoder: wgpu.GPUCommandEncoder):
        """For the GPU part of the project. This allows you to send your spatial data structure to the GPU.

        You may either copy the data structure that you implemented in part 1 as
        is (see the example code below), or you may create use your own GPU code
        to build the data structure directly on the GPU.

        If the number of accumulation frames is greater than 1, this is only called
        when the simulation has actually being updated.
        """
        return

        # 0. Make sure to do nothing if the simulation is empty, otherwise you
        #    will get errors due to the empty buffers.
        if len(simulator.radii) == 0:
            return

        # 1. For the purposes of the example, an array of random numbers (This
        # could be a numpy array containing the data structure you used for the
        # first part).
        data_structure_cpu = np.random.rand(len(simulator.radii), 2).astype(np.float32)

        # 2. Make sure we have a GPUBuffer to send the data from the CPU to the
        #    GPU We use a buffer that's double the size that we need, so that
        #    the GPU can use one half of it while we write to the other half.
        #
        #    Due to the different types of memory available to the GPU, the
        #    WebGPU specification does not allow us to use this buffer directly
        #    in our GPU code.
        if not self.data_structure_staging or self.data_structure_staging.size < 2 * data_structure_cpu.nbytes:
            self.data_structure_staging = device.create_buffer(
                label="My Data Structure (Staging Buffer)",
                size=2 * data_structure_cpu.nbytes,
                usage=wgpu.BufferUsage.COPY_SRC | wgpu.BufferUsage.MAP_WRITE)

        # 3. Make sure we have a GPUBuffer large enough to use the data
        #    structure on the GPU. This is the one that we'll use from our GPU
        #    code.
        if not self.data_structure_gpu or self.data_structure_gpu.size < data_structure_cpu.nbytes:
            self.data_structure_gpu = device.create_buffer(
                label="My Data Structure",
                size=data_structure_cpu.nbytes,
                usage=wgpu.BufferUsage.COPY_DST | wgpu.BufferUsage.STORAGE)

            # 4. Create a bind group so that we can use this buffer in our GPU code
            self.data_structure_binding = device.create_bind_group(
                layout=self.example_layout,
                entries=[
                    wgpu.BindGroupEntry(
                        binding=0,
                        resource=wgpu.BufferBinding(buffer=self.data_structure_gpu))
                ])

        # 5. Schedule the copy from the CPU to the GPU
        copy_buffer(encoder, data_structure_cpu, self.data_structure_staging, self.data_structure_gpu,
                    # On even frames, we use the first half of the staging
                    # buffer on the CPU, while the GPU might still be using the
                    # other half of the buffer, since it might still be
                    # executing the previous frame.
                    slot_id=simulator.frame_id % 2)

        # If you wish to build the data structure directly on the GPU, you'll
        # need to create your own compute pipelines and dispatch compute tasks.
        # You can see how this was done in the implementation of the Simulator
        # class.

    def gpu_bind_groups(self, simulation, device: wgpu.GPUDevice) -> list[wgpu.GPUBindGroup]:
        """Return the list of data that you want to access during the Ray Tracing computation.

        The array returned should match the one returned by gpu_bind_group_layouts.
        """
        return []

        # Return the bind group created when we created the GPU buffer
        return [self.data_structure_binding]

if __name__ == "__main__":
    from rendercanvas.auto import RenderCanvas, loop
    canvas = RenderCanvas(
        title="LMECA2170 — Homework 2",
        size=(1920, 1080),
        update_mode="continuous",
        max_fps=240,
        vsync=False)
    simulator = Simulator(canvas, loop, Homework())
    loop.run()
