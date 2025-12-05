"""LMECA2170 — Homework 2 — Spatial Datastructures."""

from simulator import Simulator, Contact, copy_buffer
import wgpu
import numpy as np
import bvh
import time

i=0

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
        bvh_tree = bvh.BVH(simulator.positions, simulator.radii)
        bvh_tree.build_quick(splits=4)
        self.bvh_tree = bvh_tree
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
        bvh_tree.positions = simulator.positions

        # if i==100:
        #     #every 10000 frames, rebuild the BVH
        #     bvh_tree.build_quick(splits=4)
        #     # print("cost after rebuild:",bvh_tree.root.update_cost())
        #     i=0

        # update the positions and optimize with rotations
        bvh_tree.update()
        N = len(bvh_tree.positions)
        # print("cost before:",bvh_tree.root.update_cost())
        bvh_tree.optimize_w_rotation(bvh_tree.root, N)
        # print("cost after:",bvh_tree.root.cost,"\n")
        
        contacts = []

        # Stack of node pairs to test
        stack = [(bvh_tree.root, bvh_tree.root)]

        def bbox_intersect(bb1, bb2):
            """AABB vs AABB intersection test."""
            return not (
                bb1[1][0] < bb2[0][0] or bb1[0][0] > bb2[1][0] or
                bb1[1][1] < bb2[0][1] or bb1[0][1] > bb2[1][1] or
                bb1[1][2] < bb2[0][2] or bb1[0][2] > bb2[1][2]
            )

        while stack:
            A, B = stack.pop()
            if A is B:
                if A.is_leaf():
                    continue
                else:
                    stack.append((A.left, A.right))
                    stack.append((A.left, A.left))
                    stack.append((A.right, A.right))
                    continue
            #Then A != B
            if A.is_leaf() and B.is_leaf():
                contact = simulator.intersect(A.item[0], B.item[0])
                if contact:
                    contacts.append(contact)
                continue
            
            if bbox_intersect(A.bbox, B.bbox):
                # Expand children
                if A.is_leaf():
                    # A leaf vs B internal
                    stack.append((A, B.left))
                    stack.append((A, B.right))
                elif B.is_leaf():
                    # B leaf vs A internal
                    stack.append((A.left, B))
                    stack.append((A.right, B))
                else:
                    # Both internal → 4 combinations
                    stack.append((A.left,  B.left))
                    stack.append((A.left,  B.right))
                    stack.append((A.right, B.left))
                    stack.append((A.right, B.right))
        return contacts

    def gpu_bind_group_layouts(self, device: wgpu.GPUDevice) -> list[wgpu.GPUBindGroupLayout]:
        """For the GPU part of the project. This allows you to describe the memory resources you'll access from your GPU implementation.

        The first binding group returned will have correspond to descriptor set
        2 in your GPU code, the second to descriptor set 3, and so on.
        """

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

        # 0. Make sure to do nothing if the simulation is empty, otherwise you
        #    will get errors due to the empty buffers.
        if len(simulator.radii) == 0:
            return

        # 1. For the purposes of the example, an array of random numbers (This
        # could be a numpy array containing the data structure you used for the
        # first part).

            # give the BVH flat representation
        bvh_tree = self.bvh_tree
        flat_bvh = self.bvh_tree.flat_bvh().flatten().astype(np.float32)
        data_structure_cpu = np.concatenate(([bvh_tree.nNode], flat_bvh)).astype(np.float32)

            # exemple given in the template
        # data_structure_cpu = np.random.rand(len(simulator.radii),2).astype(np.float32)

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
        # Return the bind group created when we created the GPU buffer
        return [self.data_structure_binding]

if __name__ == "__main__":
    from rendercanvas.auto import RenderCanvas, loop
    canvas = RenderCanvas(
        title="LMECA2170 — Homework 2",
        size=(1920, 1080),
        update_mode="continuous",
        max_fps=240,
        vsync=True)
    simulator = Simulator(canvas, loop, Homework())
    loop.run()
