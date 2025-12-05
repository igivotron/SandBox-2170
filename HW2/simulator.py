"""Simulation of a spheres bouncing around in a box."""

from rendercanvas import BaseRenderCanvas
import wgpu
from imgui_bundle import imgui
from wgpu.utils.imgui import ImguiWgpuBackend
import time
import numpy as np
from numpy.typing import ArrayLike
from pathlib import Path
import math
from typing import Optional
from dataclasses import dataclass

TIMESTEP = 1.0 / 240.0
PERFORMANCE_NUM_STEPS = 128
MAX_GPU_PERFORMANCE_QUERIES = 32

frame_count = 0
fps_start_time = time.time()

class PerformanceCounter:
    """Convenient interface to measure the time it takes to run a piece of code on the CPU."""

    def __init__(self, simulation, name: str):
        """Create a performance counter without starting the clock yet."""
        self.simulation = simulation
        self.name = name
        if self.name in simulation.measured_counts:
            self.slot_id = simulation.measured_counts[name] % PERFORMANCE_NUM_STEPS
        else:
            simulation.measured_times[self.name] = np.zeros(PERFORMANCE_NUM_STEPS, dtype=np.float32)
            simulation.measured_counts[self.name] = 0
            self.slot_id = 0

    def __enter__(self):
        """Start measuring the wall-clock execution time of a piece of code."""
        self.start = time.monotonic()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Store the result of the measurement into the Simulation."""
        if not exc_value:
            end = time.monotonic()
            self.simulation.measured_times[self.name][self.slot_id] = end - self.start
            self.simulation.measured_counts[self.name] += 1
        return False

class GPUPerformanceQueryHandler:
    """Utilities to compute the execution time of compute and render passes."""

    def __init__(self, device: wgpu.GPUDevice):
        self.query_set = device.create_query_set(
            label="GPU Performance Queries",
            type=wgpu.QueryType.timestamp,
            count=2 * MAX_GPU_PERFORMANCE_QUERIES)

        self.query_buffer = device.create_buffer(
            label="Query Result Buffer",
            size=2 * 8 * MAX_GPU_PERFORMANCE_QUERIES,
            usage=wgpu.BufferUsage.COPY_SRC | wgpu.BufferUsage.QUERY_RESOLVE)

        self.query_readback_buffer = device.create_buffer(
            label="Query Readback BUffer",
            size=4 * 8 * MAX_GPU_PERFORMANCE_QUERIES,
            usage=wgpu.BufferUsage.MAP_READ | wgpu.BufferUsage.COPY_DST)

        self.query_results = np.zeros(2 * MAX_GPU_PERFORMANCE_QUERIES, dtype=np.uint64)

        self.measured_keys = []

        # Assume nanoseconds because we don't have access to queue.get_timestamp_period
        # https://github.com/pygfx/wgpu-py/issues/656
        #
        # It seems to be that timestamp_period = 1 on NVIDIA GPUs, 10 on
        # AMD GPUs, and crazy numbers on Intel GPUs (52.08?).
        self.timestamp_period = 1.0 # Replace with the timestamp period for your GPU according to vulkaninfo
        if "intel" in device.adapter.info.device.lower():
            self.timestamp_perido = 52.0833
        if "AMD" in device.adapter.info.device:
            self.timestamp_period = 10.0

    def measure_compute(self, key: str) -> Optional[wgpu.ComputePassTimestampWrites]:
        """Return parameter to pass to begin_compute_pass to measure its execution time."""
        if len(self.measured_keys) == MAX_GPU_PERFORMANCE_QUERIES:
            return None

        start_id = 2 * len(self.measured_keys)
        self.measured_keys.append(key)

        return wgpu.ComputePassTimestampWrites(
            query_set=self.query_set,
            beginning_of_pass_write_index=start_id,
            end_of_pass_write_index=start_id + 1)

    def measure_graphics(self, key: str) -> Optional[wgpu.RenderPassTimestampWrites]:
        """Return a parameter to pass to begin_render_pass to measure its execution time."""
        if len(self.measured_keys) == MAX_GPU_PERFORMANCE_QUERIES:
            return None

        start_id = 2 * len(self.measured_keys)
        self.measured_keys.append(key)

        return wgpu.RenderPassTimestampWrites(
            query_set=self.query_set,
            beginning_of_pass_write_index=start_id,
            end_of_pass_write_index=start_id + 1)

    def submit(self, encoder: wgpu.GPUCommandEncoder, simulation):
        """Submit commands to resolve queries and read them into a CPU-accesible buffer."""
        if self.query_readback_buffer.map_state != "unmapped":
            return

        slot_id = simulation.frame_id % 2

        encoder.resolve_query_set(
            query_set=self.query_set,
            first_query=0, query_count=2 * len(self.measured_keys),
            destination=self.query_buffer,
            destination_offset=0)
        encoder.copy_buffer_to_buffer(
            source=self.query_buffer, source_offset=0,
            destination=self.query_readback_buffer,
            destination_offset=slot_id * self.query_buffer.size,
            size=self.query_buffer.size)

    def schedule_readback(self, simulation):
        """Schedule to asynchronously read back the time measurements submited in this frame."""
        old_keys = self.measured_keys
        self.measured_keys = []

        queue: wgpu.GPUQueue = simulation.device.queue
        slot_id = simulation.frame_id % 2

        promise = queue.on_submitted_work_done_async()

        @promise.then
        def on_readback(_none):
            self.query_readback_buffer.map_sync(
                wgpu.MapMode.READ,
                offset=slot_id * self.query_buffer.size,
                size=self.query_buffer.size)

            view = self.query_readback_buffer.read_mapped(
                slot_id * self.query_buffer.size,
                self.query_buffer.size,
                copy=False).cast("Q")

            for i, key in enumerate(old_keys):
                duration = (view[2 * i + 1] - view[2 * i]) * (1e-9 * self.timestamp_period)

                if key not in simulation.measured_counts:
                    simulation.measured_counts[key] = 0
                    simulation.measured_times[key] = \
                        np.zeros(PERFORMANCE_NUM_STEPS, dtype=np.float32)

                index = simulation.measured_counts[key] % PERFORMANCE_NUM_STEPS
                simulation.measured_times[key][index] = duration
                simulation.measured_counts[key] += 1

            self.query_readback_buffer.unmap()

@dataclass
class Contact:
    """Contact between two spheres, or a sphere and a wall."""

    obj_a: int
    obj_b: Optional[int]

    force: np.ndarray
    torque: np.ndarray

@dataclass
class Ball:
    """Information about a sphere."""

    pos: np.ndarray
    velocity: np.ndarray
    radius: float
    rotation: np.ndarray
    angular_velocity: np.ndarray
    index: Optional[int]

    def intersect_wall(self, pos: np.ndarray):
        """
        Check intersection with a wall.

        pos is the projection of the ball onto the plane containing the wall.
        """
        other = Ball(pos, np.zeros(3, dtype='float32'), 1.0,
                     np.array([0.0, 0.0, 0.0, 1.0], dtype='float32'),
                     np.zeros(3, dtype='float32'), None)
        return self.intersect(other)

    def intersect(self, other) -> Optional[Contact]:
        """Check if two balls intersect."""
        combined_radii = self.radius + other.radius
        delta_pos = self.pos - other.pos
        squared_dist = np.dot(delta_pos, delta_pos)

        if squared_dist >= combined_radii * combined_radii:
            return None

        dist = math.sqrt(squared_dist)
        fc = 1e4 * (combined_radii / dist - 1.0)

        delta_vel = self.velocity - other.velocity
        h = np.dot(delta_pos, delta_vel) / squared_dist

        fc = max(fc - 5.0 * h, 0.0)
        p = (self.radius * self.angular_velocity + \
             other.radius * other.angular_velocity) / \
             combined_radii

        delta_vel -= h * delta_pos + np.cross(p, delta_pos)
        ft = min(0.5, 0.05 * abs(fc) * dist / max(1e-3, np.linalg.norm(delta_vel)))

        torque = (ft / dist) * np.cross(delta_pos, delta_vel)

        return Contact(self.index, other.index,
                       fc * delta_pos - ft * delta_vel,
                       torque)

CAMERA_TYPE = np.dtype([
    ('view', 'f4', (4, 4)),
    ('projection', 'f4', (4, 4))
])

MATERIAL_TYPE = np.dtype([
    ('albedo', 'f4', (3,)),
    ('roughness', 'f4', (1,)),
    ('emittance', 'f4', (3,)),
    ('metalness', 'f4', (1,))
])

VERTEX_LAYOUT = wgpu.VertexBufferLayout(
    array_stride=8 * 4,
    attributes=[
        # location
        wgpu.VertexAttribute(
            format=wgpu.VertexFormat.float32x3,
            offset=0,
            shader_location=0),
        # normal
        wgpu.VertexAttribute(
            format=wgpu.VertexFormat.float32x3,
            offset=4*4,
            shader_location=1),
    ])

KEY_MAP = {
    "ArrowDown": imgui.Key.down_arrow,
    "ArrowUp": imgui.Key.up_arrow,
    "ArrowLeft": imgui.Key.left_arrow,
    "ArrowRight": imgui.Key.right_arrow,
    "Backspace": imgui.Key.backspace,
    "CapsLock": imgui.Key.caps_lock,
    "Delete": imgui.Key.delete,
    "End": imgui.Key.end,
    "Enter": imgui.Key.enter,
    "Escape": imgui.Key.escape,
    "F1": imgui.Key.f1,
    "F2": imgui.Key.f2,
    "F3": imgui.Key.f3,
    "F4": imgui.Key.f4,
    "F5": imgui.Key.f5,
    "F6": imgui.Key.f6,
    "F7": imgui.Key.f7,
    "F8": imgui.Key.f8,
    "F9": imgui.Key.f9,
    "F10": imgui.Key.f10,
    "F11": imgui.Key.f11,
    "F12": imgui.Key.f12,
    "Home": imgui.Key.home,
    "Insert": imgui.Key.insert,
    "Alt": imgui.Key.left_alt,
    "Control": imgui.Key.left_ctrl,
    "Shift": imgui.Key.left_shift,
    "Meta": imgui.Key.left_super,
    "NumLock": imgui.Key.num_lock,
    "PageDown": imgui.Key.page_down,
    "PageUp": imgui.Key.page_up,
    "Pause": imgui.Key.pause,
    "PrintScreen": imgui.Key.print_screen,
    "ScrollLock": imgui.Key.scroll_lock,
    "Tab": imgui.Key.tab,
}

KEY_MAP_MOD = {
    "Shift": imgui.Key.mod_shift,
    "Control": imgui.Key.mod_ctrl,
    "Alt": imgui.Key.mod_alt,
    "Meta": imgui.Key.mod_super,
}

def camera_data(view, projection):
    """Return the data for a camera, with the right memory layout to send it to the GPU."""
    return np.array((np.transpose(view), np.transpose(projection)), dtype=CAMERA_TYPE)

def pad_size(n):
    """Pad a number to the next multiple of 8. Used for memory alignment purposes."""
    return (n + 7) & ~7

def create_sphere(device: wgpu.GPUDevice, n_azimuth=40, n_elevation=40):
    """Build a mesh of sphere, and upload it to the GPU."""
    vertices = []
    indices  = []

    for i in range(n_elevation):
        elevation = (math.pi / 2.0) + (math.pi * i) / (n_elevation - 1.0)
        for j in range(n_azimuth):
            azimuth = 2 * math.pi * j / (n_azimuth - 1.0)
            x = math.cos(azimuth) * math.cos(elevation)
            y = math.sin(elevation)
            z = math.sin(azimuth) * math.cos(elevation)

            vertices += [x, y, z, 0.0, x, y, z, 0.0]

    for i in range(n_elevation - 1):
        for j in range(n_azimuth):
            next_j = (j + 1) % n_azimuth
            next_i = n_elevation - 2 if i == n_elevation else i + 1

            indices += [
                j + (i * n_azimuth),
                next_j + (i * n_azimuth),
                j + (next_i * n_azimuth),

                j + (next_i * n_azimuth),
                next_j + (i * n_azimuth),
                next_j + (next_i * n_azimuth),
            ]

    vbo = device.create_buffer_with_data(
        label="Sphere VBO",
        data=np.array(vertices, dtype=np.float32),
        usage=wgpu.BufferUsage.VERTEX,
    )

    ibo = device.create_buffer_with_data(
        label="Sphere IBO",
        data=np.array(indices, dtype=np.uint16),
        usage=wgpu.BufferUsage.INDEX,
    )

    return vbo, ibo, len(indices)

def create_box(device: wgpu.GPUDevice, size=(1.0, 1.0, 1.0)) -> (wgpu.GPUBuffer, wgpu.GPUBuffer):
    """Create a mesh of sphere and upload it to the GPU."""
    sx, sy, sz = size

    vertices = [
        # -X
        -sx, -sy, -sz, 0.0, +1.0, 0.0, 0.0, 0.0,
        -sx, +sy, -sz, 0.0, +1.0, 0.0, 0.0, 0.0,
        -sx, +sy, +sz, 0.0, +1.0, 0.0, 0.0, 0.0,
        -sx, -sy, +sz, 0.0, +1.0, 0.0, 0.0, 0.0,

        # +X
        +sx, -sy, -sz, 0.0, -1.0, 0.0, 0.0, 0.0,
        +sx, +sy, -sz, 0.0, -1.0, 0.0, 0.0, 0.0,
        +sx, +sy, +sz, 0.0, -1.0, 0.0, 0.0, 0.0,
        +sx, -sy, +sz, 0.0, -1.0, 0.0, 0.0, 0.0,

        # -Y
        -sx, -sy, -sz, 0.0, 0.0, +1.0, 0.0, 0.0,
        +sx, -sy, -sz, 0.0, 0.0, +1.0, 0.0, 0.0,
        +sx, -sy, +sz, 0.0, 0.0, +1.0, 0.0, 0.0,
        -sx, -sy, +sz, 0.0, 0.0, +1.0, 0.0, 0.0,

        # +Y
        -sx, +sy, -sz, 0.0, 0.0, -1.0, 0.0, 0.0,
        +sx, +sy, -sz, 0.0, 0.0, -1.0, 0.0, 0.0,
        +sx, +sy, +sz, 0.0, 0.0, -1.0, 0.0, 0.0,
        -sx, +sy, +sz, 0.0, 0.0, -1.0, 0.0, 0.0,

        # -Z
        -sx, -sy, -sz, 0.0, 0.0, 0.0, +1.0, 0.0,
        +sx, -sy, -sz, 0.0, 0.0, 0.0, +1.0, 0.0,
        +sx, +sy, -sz, 0.0, 0.0, 0.0, +1.0, 0.0,
        -sx, +sy, -sz, 0.0, 0.0, 0.0, +1.0, 0.0,

        # +Z
        -sx, -sy, +sz, 0.0, 0.0, 0.0, -1.0, 0.0,
        +sx, -sy, +sz, 0.0, 0.0, 0.0, -1.0, 0.0,
        +sx, +sy, +sz, 0.0, 0.0, 0.0, -1.0, 0.0,
        -sx, +sy, +sz, 0.0, 0.0, 0.0, -1.0, 0.0,
    ]

    indices = [
        # -X
        0, 1, 2, 0, 2, 3,

        # +X
        4, 6, 5, 4, 7, 6,

        # -Y
        9, 8, 10, 8, 11, 10,

        # +Y
        12, 13, 14, 12, 14, 15,

        # -Z
        16, 17, 18, 16, 18, 19,

        # +Z
        20, 22, 21, 20, 23, 22
    ]

    vertices = np.array(vertices, dtype=np.float32)
    indices = np.array(indices, dtype=np.uint16)

    vbo = device.create_buffer_with_data(label="Box VBO", data=vertices,
                                         usage=wgpu.BufferUsage.VERTEX)

    ibo = device.create_buffer_with_data(label="Box IBO", data=indices,
                                         usage=wgpu.BufferUsage.INDEX)

    return vbo, ibo

def perspective_matrix(fovy: float, aspect_ratio: float, near_z: float) -> np.ndarray:
    """Build a 4×4 projection matrix.

    Uses the recommendations from 'Tightening the Precision of Perspective
    Rendering', Upchurch and Desbrun, 2012, to set the far plane infinitely far
    away to improve precision.
    """
    tan_up = math.tan(fovy / 2)
    tan_down = -tan_up

    tan_right = aspect_ratio * tan_up
    tan_left = -tan_right

    tan_width = tan_right - tan_left
    tan_height = tan_up - tan_down

    return np.array([
        [2.0 / tan_width, 0.0, 0.0, 0.0],
        [0.0, 2.0 / tan_height, 0.0, 0.0],
        [0.0, 0.0, -1.0, -near_z],
        [0.0, 0.0, -1.0, 0.0],
    ])

def random_materials(n):
    """Return an array of n random materials."""
    out = np.zeros(n, dtype=MATERIAL_TYPE)

    out['albedo'] = np.random.rand(n, 3).astype(np.float32)
    out['roughness'] = np.random.rand(n, 1).astype(np.float32)
    out['emittance'] = np.zeros((n, 3), dtype=np.float32)
    out['metalness'] = np.random.choice([0.0, 1.0], size=(n, 1)).astype(np.float32)

    return out

def copy_buffer(encoder: wgpu.GPUCommandEncoder,
                array: np.ndarray, staging: wgpu.GPUBuffer, gpu_buffer: wgpu.GPUBuffer,
                slot_id: int = 0):
    """Copy data from a numpy array to the GPU.

    Data is directly copied into the staging buffer, and a command is written
    for the GPU to copy from the staging buffer to the destination buffer. The
    copy will be completed when the command buffer is executed.
    """
    if array.nbytes == 0:
        # Need a special case because wgpu doesn't allow mapping an empty range
        return

    staging.map_sync(
        mode=wgpu.MapMode.WRITE,
        offset=slot_id * pad_size(array.nbytes),
        size=array.nbytes)
    staging.write_mapped(
        data=array,
        buffer_offset=slot_id * pad_size(array.nbytes))
    staging.unmap()

    encoder.copy_buffer_to_buffer(
        source=staging,
        source_offset=slot_id * pad_size(array.nbytes),
        destination=gpu_buffer,
        destination_offset=0,
        size=array.nbytes)

class Simulator:
    """Simulation of spheres bouncing in a box."""

    def __init__(self, canvas: BaseRenderCanvas, loop, user_data):
        """Initialize the GPU resources for an empty simulation."""
        self.user_data = user_data

        self.adapter: wgpu.GPUAdapter = wgpu.gpu.request_adapter_async(
            power_preference="high-performance",
            loop=loop
        ).sync_wait()

        self.supports_barriers = "subgroup-barrier" in self.adapter.features

        timestamp_features = [
            'timestamp-query',
            'timestamp-query-inside-encoders',
            'timestamp-query-inside-passes'
        ]

        self.supports_gpu_time_queries = all(feature in self.adapter.features for feature in timestamp_features)

        features = timestamp_features if self.supports_gpu_time_queries else []
        if self.supports_barriers:
            features.append("subgroup")
            features.append("subgroup-barrier")

        self.device: wgpu.GPUDevice = self.adapter.request_device_sync(
            required_features=features)

        self.gpu_timer = GPUPerformanceQueryHandler(self.device) if self.supports_gpu_time_queries else None

        self.canvas: BaseRenderCanvas = canvas

        self.measured_times = {}
        self.measured_counts = {}

        self.present_context = self.canvas.get_context("wgpu")
        self.color_format: wgpu.TextureFormat = self.present_context.get_preferred_format(self.adapter)
        self.present_context.configure(device=self.device, format=self.color_format)

        imgui.create_context()
        self.imgui_backend = ImguiWgpuBackend(self.device, self.color_format)

        self.uniform_layout = self.device.create_bind_group_layout(
            entries=[
                wgpu.BindGroupLayoutEntry(
                    binding=0,
                    visibility=wgpu.ShaderStage.COMPUTE | wgpu.ShaderStage.VERTEX | wgpu.ShaderStage.COMPUTE,
                    buffer=wgpu.BufferBindingLayout(type=wgpu.BufferBindingType.uniform)
                )
            ]
        )

        self.camera_staging = self.device.create_buffer(
            label="Camera Staging Buffer",
            size=2 * CAMERA_TYPE.itemsize,
            usage=wgpu.BufferUsage.MAP_WRITE |  wgpu.BufferUsage.COPY_SRC,
        )

        self.camera_ubo = self.device.create_buffer(
            label="Camera UBO",
            size=CAMERA_TYPE.itemsize,
            usage=wgpu.BufferUsage.UNIFORM |  wgpu.BufferUsage.COPY_DST,
        )

        self.rt_camera_staging = self.device.create_buffer(
            label="Camera Staging Buffer (RT)",
            size=2 * CAMERA_TYPE.itemsize,
            usage=wgpu.BufferUsage.MAP_WRITE |  wgpu.BufferUsage.COPY_SRC,
        )

        self.rt_camera_ubo = self.device.create_buffer(
            label="Camera UBO (RT)",
            size=CAMERA_TYPE.itemsize,
            usage=wgpu.BufferUsage.UNIFORM |  wgpu.BufferUsage.COPY_DST,
        )

        self.rt_reset_staging = self.device.create_buffer(
            label="Ray Tracing — Reset (Staging Buffer) ",
            size=4,
            usage=wgpu.BufferUsage.MAP_WRITE |  wgpu.BufferUsage.COPY_SRC,
        )

        self.rt_reset_ubo = self.device.create_buffer(
            label="Ray Tracing — Reset (UBO)",
            size=4,
            usage=wgpu.BufferUsage.UNIFORM |  wgpu.BufferUsage.COPY_DST,
        )

        self.uniform_binding = self.device.create_bind_group(
            label="Camera Binding",
            layout=self.uniform_layout,
            entries=[
                wgpu.BindGroupEntry(binding=0, resource=wgpu.BufferBinding(
                    buffer=self.camera_ubo,
                    size=CAMERA_TYPE.itemsize
                ))
            ]
        )

        self.ray_tracing_camera_binding = self.device.create_bind_group(
            label="Camera Binding (RT)",
            layout=self.uniform_layout,
            entries=[
                wgpu.BindGroupEntry(binding=0, resource=wgpu.BufferBinding(
                    buffer=self.rt_camera_ubo,
                    size=CAMERA_TYPE.itemsize
                ))
            ]
        )

        self.raster_pipeline_layout = self.device.create_pipeline_layout(bind_group_layouts=[
            self.uniform_layout
        ])

        self._create_box_pipeline()
        self._create_ball_pipeline()
        self._create_path_tracing_pipeline()

        self.depth_texture: Optional[wgpu.GPUTexture] = None
        self.depth_view: Optional[wgpu.GPUTextureView] = None

        self.rt_texture: Optional[wgpu.GPUTexture] = None
        self.rt_view: Optional[wgpu.GPUTextureView] = None

        self.frame_id = 0
        self.last_frame_time: Optional[float] = None

        self.box_vbo, self.box_ibo = create_box(self.device, size=(10.0, 10.0, 10.0))
        self.sphere_vbo, self.sphere_ibo, self.sphere_num_indices = create_sphere(self.device)

        self._load_sobol_directions()

        self.clear()

        self.elapsed_time = 0.0

        self.input_min_radius = 0.5
        self.input_max_radius = 1.5
        self.input_num_balls = 20
        self.ray_tracing = False
        self.rt_accumulation_frames = 1
        self.rt_step_id = 0
        self.force_rt_update = False

        @self.canvas.add_event_handler("pointer_move")
        def on_mouse_move(ev):
            self.imgui_backend.io.add_mouse_pos_event(ev["x"], ev["y"])

        @self.canvas.add_event_handler("pointer_up", "pointer_down")
        def on_mouse(ev):
            self.imgui_backend.io.add_mouse_button_event(
                ev["button"] - 1,
                ev["event_type"] == "pointer_down"
            )

        @self.canvas.add_event_handler("resize")
        def on_resize(event):
            self.imgui_backend.io.display_size = (event["width"], event["height"])

        @self.canvas.add_event_handler("key_up", "key_down")
        def on_key(event):
            # Copied from wgpu.utils.imgui
            down = event["event_type"] == "key_down"

            key_name = event["key"]
            if key_name in KEY_MAP:
                key = KEY_MAP[key_name]
            else:
                key = ord(key_name.lower())
                if key >= 48 and key <= 57:  # numbers 0-9
                    key = imgui.Key(imgui.Key._0.value + (key - 48))
                elif key >= 97 and key <= 122:  # letters a-z
                    key = imgui.Key(imgui.Key.a.value + (key - 97))
                else:
                    return  # Unknown key: {key_name}

            self.imgui_backend.io.add_key_event(key, down)

            if key_name in KEY_MAP_MOD:
                key = KEY_MAP_MOD[key_name]
                self.imgui_backend.io.add_key_event(key, down)

            if self.imgui_backend.io.want_capture_keyboard:
                event["stop_propagation"] = True

        @self.canvas.add_event_handler("char")
        def on_text(event):
            self.imgui_backend.io.add_input_characters_utf8(event["char_str"])

        @self.canvas.request_draw
        def draw():
            self.render()

    def clear(self):
        """Reset the simulation to contain zero spheres."""
        self.positions: np.ndarray = np.zeros((0, 3), order='C', dtype=np.float32)
        self.rotations: np.ndarray = np.zeros((0, 4), order='C', dtype=np.float32)
        self.velocities: np.ndarray =  np.zeros((0, 3), order='C', dtype=np.float32)
        self.angular_velocities: np.ndarray = np.zeros((0, 3), order='C', dtype=np.float32)
        self.radii: np.ndarray = np.array([], dtype=np.float32)
        self.materials: np.ndarray = np.array([], dtype=MATERIAL_TYPE)
        self._create_physics_buffers()
        self.user_data.on_size_changed(self)

    def spawn_balls(self, positions, rotations, radii):
        """Add objects at the given locations to the simulation."""
        pos_array = np.array(positions, order='C', dtype='float32')
        rot_array = np.array(rotations, order='C', dtype='float32')
        velocities = np.zeros((len(radii), 3), dtype='float32', order='C')
        angular_velocities = np.zeros((len(radii), 3), dtype='float32', order='C')
        radii_array = np.array(radii, order='C', dtype='float32')
        materials = random_materials(len(radii))

        if len(self.radii) == 0:
            self.positions = pos_array
            self.rotations = rot_array
            self.velocities = velocities
            self.angular_velocities = angular_velocities
            self.radii = radii_array
            self.materials = materials
        else:
            self.positions = np.concat((self.positions, pos_array))
            self.rotations = np.concat((self.rotations, rot_array))
            self.velocities = np.concat((self.velocities, velocities))
            self.angular_velocities = np.concat((self.angular_velocities, angular_velocities))
            self.radii = np.concat((self.radii, radii_array))
            self.materials = np.concat((self.materials, materials))

        self._create_physics_buffers()
        self.user_data.on_size_changed(self)

    def spawn_random_balls(self, n, min_radius=0.1, max_radius=0.5):
        """Add random objects to the simulation."""
        radii = min_radius + (max_radius - min_radius) * np.random.rand(n).astype(np.float32)
        positions = np.random.random((n, 3)) * 20 - 10.0
        rotations = np.repeat(np.array([[0, 0, 0, 1]], dtype=np.float32), n, axis=0)
        self.spawn_balls(positions, rotations, radii)

    def measure_cpu(self, name: str) -> PerformanceCounter:
        """Return an object to track the execution time of a piece of code on the CPU."""
        return PerformanceCounter(self, name)

    def measure_gpu_compute(self, name: str) -> Optional[wgpu.RenderPassTimestampWrites]:
        """Return a parameter to pass to begin_render_pass to measure its execution time."""
        if self.gpu_timer:
            return self.gpu_timer.measure_compute(name)
        else:
            return None

    def measure_gpu_graphics(self, name: str) -> Optional[wgpu.RenderPassTimestampWrites]:
        """Return a parameter to pass to begin_compute_pass to measure its execution time."""
        if self.gpu_timer:
            return self.gpu_timer.measure_graphics(name)
        else:
            return None

    def measurements_for(self, key) -> np.ndarray:
        """Return an ordered numpy array of the measured execution times with a given key."""
        times = self.measured_times[key]
        cur_id = self.measured_counts[key] % PERFORMANCE_NUM_STEPS
        return np.concat((times[cur_id:], times[:cur_id]))

    def average_time(self, key) -> float:
        """Return the average execution time for the code blocks associated with a key."""
        return np.sum(self.measured_times[key]) / min(PERFORMANCE_NUM_STEPS, self.measured_counts[key])

    def average_time_string(self, key) -> str:
        """Return a string describing average execution time along with units."""
        time = self.average_time(key)
        if time > 1.0:
            return f"{time:.2f} s"
        elif time > 1e-3:
            return f"{(time * 1e3):.2f} ms"
        elif time > 1e-6:
            return f"{(time * 1e6):.2f} us"
        else:
            return f"{(time * 1e9):.2f} ns"

    def draw_gui(self):
        """Draw the GUI to a buffer, and update the settings according to user input."""
        self.imgui_backend.io.display_size = self.canvas.get_logical_size()
        self.imgui_backend.io.display_framebuffer_scale = (
            self.canvas.get_pixel_ratio(),
            self.canvas.get_pixel_ratio(),
        )

        imgui.new_frame()

        imgui.begin("Debugging Information")
        imgui.text(f"# Balls: {len(self.radii)}")
        new, value = imgui.slider_float("Min Radius", self.input_min_radius, 0.05, 3.0)
        if new:
            self.input_min_radius = min(self.input_max_radius, value)

        new, value = imgui.slider_float("Max radius", self.input_max_radius, 0.05, 3.0)
        if new:
            self.input_max_radius = max(self.input_min_radius, value)

        new, value = imgui.input_int("Num balls", self.input_num_balls, step=1, step_fast=50)
        if new:
            self.input_num_balls = max(value, 1)

        if imgui.button("Spawn Balls"):
            self.spawn_random_balls(
                self.input_num_balls,
                min_radius=self.input_min_radius,
                max_radius=self.input_max_radius)

        if imgui.button("Clear"):
            self.clear()

        new, value = imgui.checkbox("Ray Tracing", self.ray_tracing)
        if new and self.ray_tracing != value:
            self.force_rt_update
            self.ray_tracing = value

        imgui.text("# Accumulation Frames")
        new, value = imgui.input_int(" ", self.rt_accumulation_frames, 1, 1)
        if new:
            self.rt_accumulation_frames = max(value, 1)

        for key in self.measured_times:
            imgui.plot_lines(key, self.measurements_for(key), graph_size=(200, 50), scale_min=0.0)
            imgui.text(f"Avg. {key}: {self.average_time_string(key)}")

        imgui.end()

        imgui.end_frame()
        imgui.render()
        return imgui.get_draw_data()

    def find_intersections(self):
        """Compute a list of all pairs of balls that intersect."""
        with self.measure_cpu("intersection"):
            return self.user_data.find_intersections(self)

    def ball(self, i: int) -> Ball:
        """Return all data relative to ball at index i."""
        return Ball(self.positions[i], self.velocities[i],
                    self.radii[i], self.rotations[i], self.angular_velocities[i],
                    i)

    def intersect(self, i: int, j: int):
        """Return a Contact object if two spheres intersect. None otherwise."""
        return self.ball(i).intersect(self.ball(j))

    def apply_contact(self, contact: Contact):
        """Apply contact forces to resolve a contact."""
        i, j = contact.obj_a, contact.obj_b
        m_a = self.radii[i] * self.radii[i]
        self.velocities[i] += TIMESTEP * contact.force / m_a
        self.angular_velocities[i] += TIMESTEP * contact.torque / (0.4 * m_a * self.radii[i])

        if j:
            m_b = self.radii[j] * self.radii[j]
            self.velocities[j] -= TIMESTEP * contact.force / m_b
            self.angular_velocities[j] -= TIMESTEP * contact.torque / (0.4 * m_b * self.radii[j])

    def update(self, dt):
        """Simulate rigidbody objects for dt seconds. Actual updates only occur at fixed interval."""


        global frame_count, fps_start_time

        # === Mesure des FPS ===
        frame_count += 1
        now = time.time()
        elapsed = now - fps_start_time
        if elapsed >= 1.0:  # chaque seconde
            fps = frame_count / elapsed
            print(f"FPS: {fps:.2f}")
            # Optionnel : enregistrer dans un fichier CSV
            with open("./logs/fps_log.txt", "a") as f:
                f.write(f"{now},{fps:.2f}\n")
            frame_count = 0
            fps_start_time = now



        self.elapsed_time += dt
        i = 0

        MAX_STEPS = 4
        while self.elapsed_time >= TIMESTEP:
            self.step()
            self.elapsed_time -= TIMESTEP
            i += 1
            if i >= MAX_STEPS:
                self.elapsed_time = self.elapsed_time % TIMESTEP
                break

    def step(self):
        """Simulate objects for TIMESTEP seconds."""
        if len(self.radii) == 0:
            return

        for i in range(len(self.radii)):
            ball = self.ball(i)
            domain_size = 10.0
            walls = np.array([
                [-domain_size - 1.0, ball.pos[1], ball.pos[2]],
                [+domain_size + 1.0, ball.pos[1], ball.pos[2]],
                [ball.pos[0], -domain_size - 1.0, ball.pos[2]],
                [ball.pos[0], +domain_size + 1.0, ball.pos[2]],
                [ball.pos[0], ball.pos[1], -domain_size - 1.0],
                [ball.pos[0], ball.pos[1], +domain_size + 1.0],
            ], dtype='float32')

            for wall in walls:
                contact = ball.intersect_wall(wall)
                if contact:
                    self.apply_contact(contact)

        for contact in self.find_intersections():
            self.apply_contact(contact)

        gravity = np.array([0, -9.81, 0])

        self.positions += TIMESTEP * self.velocities
        self.velocities += TIMESTEP * gravity

        drotation = 0.5 * np.einsum('ij,oj->oi',
                                    np.array([[0.0, 0.0, 0.0],
                                              [1.0, 0.0, 0.0],
                                              [0.0, 1.0, 0.0],
                                              [0.0, 0.0, 1.0]]),
                                    self.angular_velocities)

        weights_a = np.array([[[-1.0, 0.0, 0.0, 0.0], [0.0, -1.0, 0.0, 0.0], [0.0, 0.0, -1.0, 0.0], [0.0, 0.0, 0.0, -1.0]],
                              [[+1.0, 0.0, 0.0, 0.0], [0.0, +1.0, 0.0, 0.0], [0.0, 0.0, +1.0, 0.0], [0.0, 0.0, 0.0, -1.0]],
                              [[+1.0, 0.0, 0.0, 0.0], [0.0, -1.0, 0.0, 0.0], [0.0, 0.0, +1.0, 0.0], [0.0, 0.0, 0.0, +1.0]],
                              [[+1.0, 0.0, 0.0, 0.0], [0.0, +1.0, 0.0, 0.0], [0.0, 0.0, -1.0, 0.0], [0.0, 0.0, 0.0, -1.0]]])

        weights_b = np.array([[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
                              [[0.0, 1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0]],
                              [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]],
                              [[0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]]])

        self.rotations += TIMESTEP * \
            np.einsum('qij,qij->qi',
                      np.einsum('ijk,qk->qij', weights_a, drotation),
                      np.einsum('ijk,qk->qij', weights_b, self.rotations))

        inv_norms = 1.0 / np.linalg.norm(self.rotations, axis=1)
        self.rotations = np.einsum('i,ij->ij', inv_norms, self.rotations)

    def render(self):
        """Update the simulation and draw the next frame."""
        self.frame_id += 1

        run_update = not self.ray_tracing or (self.rt_step_id + 1 >= self.rt_accumulation_frames)
        if run_update:
            current_time = time.monotonic()
            if self.last_frame_time:
                self.update(current_time - self.last_frame_time)
                self.last_frame_time = current_time
            else:
                self.last_frame_time = current_time

            self.rt_step_id = 0
        else:
            self.rt_step_id += 1

        self._update_depth_texture()
        self._update_rt_texture()

        w, h = self.canvas.get_physical_size()
        aspect_ratio = w / h
        view = np.array([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, -20.0],
            [0.0, 0.0, 0.0, 1.0]
        ], dtype=np.float32)

        projection = perspective_matrix(math.pi / 2, aspect_ratio, 1e-1)

        camera = camera_data(view, projection)

        enc = self.device.create_command_encoder()

        copy_buffer(enc, camera, self.camera_staging, self.camera_ubo, slot_id=self.frame_id % 2)

        inv_camera = camera_data(np.linalg.inv(view), np.linalg.inv(projection))
        copy_buffer(enc, inv_camera, self.rt_camera_staging, self.rt_camera_ubo, slot_id=self.frame_id % 2)

        if run_update or self.force_rt_update:
            copy_buffer(
                enc,
                self.positions, self.positions_staging, self.positions_attr,
                slot_id=self.frame_id % 2)
            copy_buffer(
                enc,
                self.rotations, self.rotations_staging, self.rotations_attr,
                slot_id=self.frame_id % 2)
            copy_buffer(
                enc,
                self.velocities, self.velocities_staging, self.velocities_attr,
                slot_id=self.frame_id % 2)
            copy_buffer(
                enc,
                self.angular_velocities, self.angular_velocities_staging, self.angular_velocities_attr,
                slot_id=self.frame_id % 2)

        if self.ray_tracing:
            copy_buffer(enc, np.array([0 if self.rt_step_id == 0 else 1], dtype='i4'),
                        self.rt_reset_staging, self.rt_reset_ubo)
            if run_update or self.force_rt_update:
                self.user_data.update_gpu_data(self, self.device, enc)

        self.force_rt_update = False

        current_texture: wgpu.GPUTexture = self.present_context.get_current_texture()
        current_view = current_texture.create_view()

        if not self.ray_tracing or len(self.radii) == 0:
            # WGPU won't let you attach empty buffers, so this path is used
            # until a sphere has been spawned
            render_pass = enc.begin_render_pass(
                color_attachments=[
                    wgpu.RenderPassColorAttachment(
                        view=current_view,
                        resolve_target=None,
                        clear_value=(0.0, 0.0, 0.0, 1.0),
                        load_op=wgpu.LoadOp.clear,
                        store_op=wgpu.StoreOp.store,
                    ),
                ],
                depth_stencil_attachment=wgpu.RenderPassDepthStencilAttachment(
                    view=self.depth_view,
                    depth_clear_value=1.0,
                    depth_load_op=wgpu.LoadOp.clear,
                    depth_store_op=wgpu.StoreOp.store),
                timestamp_writes=self.measure_gpu_graphics("rasterizer"))

            render_pass.set_pipeline(self.box_pipeline)
            render_pass.set_bind_group(0, self.uniform_binding)
            render_pass.set_index_buffer(self.box_ibo, wgpu.IndexFormat.uint16)
            render_pass.set_vertex_buffer(0, self.box_vbo)
            render_pass.draw_indexed(6*6)

            render_pass.set_pipeline(self.ball_pipeline)
            render_pass.set_index_buffer(self.sphere_ibo, wgpu.IndexFormat.uint16)
            render_pass.set_vertex_buffer(0, self.sphere_vbo)
            render_pass.set_vertex_buffer(1, self.positions_attr)
            render_pass.set_vertex_buffer(2, self.rotations_attr)
            render_pass.set_vertex_buffer(3, self.radii_buffer)
            render_pass.set_vertex_buffer(4, self.materials_buffer)
            render_pass.draw_indexed(
                index_count=self.sphere_num_indices,
                instance_count=len(self.radii),
            )

            render_pass.end()
        else:
            compute_pass = enc.begin_compute_pass(
                label="Ray Tracing",
                timestamp_writes=self.measure_gpu_compute("ray tracing"))
            compute_pass.set_pipeline(self.ray_tracing_pipeline)
            compute_pass.set_bind_group(0, self.ray_tracing_camera_binding)
            compute_pass.set_bind_group(1, self.ray_tracing_bindings)
            for i, group in enumerate(self.user_data.gpu_bind_groups(self, self.device)):
                compute_pass.set_bind_group(2 + i, group)

            w, h, _ = self.rt_texture.size
            if self.supports_barriers:
                compute_pass.dispatch_workgroups(w // 4, h, 1)
            else:
                compute_pass.dispatch_workgroups(((w // 4) + 31) // 32, h, 1)

            compute_pass.end()

            blit_pass = enc.begin_render_pass(
                color_attachments=[
                    wgpu.RenderPassColorAttachment(
                        view=current_view,
                        resolve_target=None,
                        clear_value=(0.0, 0.0, 0.0, 1.0),
                        load_op=wgpu.LoadOp.clear,
                        store_op=wgpu.StoreOp.store,
                    ),
                ]
            )

            blit_pass.set_pipeline(self.blit_pipeline)
            blit_pass.set_bind_group(0, self.blit_bindings)
            blit_pass.draw(3)

            blit_pass.end()

        gui_pass = enc.begin_render_pass(
            color_attachments=[
                wgpu.RenderPassColorAttachment(
                    view=current_view,
                    resolve_target=None,
                    clear_value=(0.0, 0.0, 0.0, 1.0),
                    load_op=wgpu.LoadOp.load,
                    store_op=wgpu.StoreOp.store,
                ),
            ]
        )

        self.imgui_backend.render(self.draw_gui(), gui_pass)
        gui_pass.end()

        if self.gpu_timer:
            self.gpu_timer.submit(enc, self)

        self.device.queue.submit([enc.finish()])

        if self.gpu_timer:
            self.gpu_timer.schedule_readback(self)

    def _create_box_pipeline(self):
        vs = self.device.create_shader_module(label="Box_vert", code=Path("shaders/box.vert").read_text())
        fs = self.device.create_shader_module(label="Box_frag", code=Path("shaders/box.frag").read_text())

        pipeline_layout = self.device.create_pipeline_layout(bind_group_layouts=[
            self.uniform_layout
        ])

        self.box_pipeline = self.device.create_render_pipeline(
            layout=self.raster_pipeline_layout,
            vertex=wgpu.VertexState(module=vs, entry_point="main", buffers=[
                VERTEX_LAYOUT
            ]),
            primitive=wgpu.PrimitiveState(
                cull_mode=wgpu.CullMode.back
            ),
            depth_stencil=wgpu.DepthStencilState(
                depth_write_enabled=True,
                format=wgpu.TextureFormat.depth16unorm,
                depth_compare=wgpu.CompareFunction.less,
            ),
            multisample=wgpu.MultisampleState(count=1),
            fragment=wgpu.FragmentState(module=fs, entry_point="main", targets=[
                wgpu.ColorTargetState(format=self.color_format)
            ]))

    def _create_ball_pipeline(self):
        vs = self.device.create_shader_module(label="Ball_vert", code=Path("shaders/sphere.vert").read_text())
        fs = self.device.create_shader_module(label="Ball_frag", code=Path("shaders/box.frag").read_text())

        self.ball_pipeline = self.device.create_render_pipeline(
            layout=self.raster_pipeline_layout,
            vertex=wgpu.VertexState(module=vs, entry_point="main", buffers=[
                VERTEX_LAYOUT,

                # sphere center
                wgpu.VertexBufferLayout(
                    array_stride=3 * 4,
                    step_mode=wgpu.VertexStepMode.instance,
                    attributes=[
                        wgpu.VertexAttribute(format=wgpu.VertexFormat.float32x3,
                                             offset=0,
                                             shader_location=2)
                    ]
                ),

                # rotation
                wgpu.VertexBufferLayout(
                    array_stride=4 * 4,
                    step_mode=wgpu.VertexStepMode.instance,
                    attributes=[
                        wgpu.VertexAttribute(format=wgpu.VertexFormat.float32x4,
                                             offset=0,
                                             shader_location=3)
                    ]
                ),

                # radius
                wgpu.VertexBufferLayout(
                    array_stride=4,
                    step_mode=wgpu.VertexStepMode.instance,
                    attributes=[
                        wgpu.VertexAttribute(format=wgpu.VertexFormat.float32,
                                             offset=0,
                                             shader_location=4)
                    ]
                ),

                # Material
                wgpu.VertexBufferLayout(
                    array_stride=MATERIAL_TYPE.itemsize,
                    step_mode=wgpu.VertexStepMode.instance,
                    attributes=[
                        # albedo
                        wgpu.VertexAttribute(
                            format=wgpu.VertexFormat.float32x3,
                            offset=0,
                            shader_location=5),

                        # roughness
                        wgpu.VertexAttribute(
                            format=wgpu.VertexFormat.float32,
                            offset=3 * 4,
                            shader_location=6),

                        # emittance
                        wgpu.VertexAttribute(
                            format=wgpu.VertexFormat.float32x3,
                            offset=4 * 4,
                            shader_location=7),

                        # metalness
                        wgpu.VertexAttribute(
                            format=wgpu.VertexFormat.float32,
                            offset=7 * 4,
                            shader_location=8),
                    ]
                )
            ]),
            primitive=wgpu.PrimitiveState(
                cull_mode=wgpu.CullMode.back
            ),
            depth_stencil=wgpu.DepthStencilState(
                depth_write_enabled=True,
                format=wgpu.TextureFormat.depth16unorm,
                depth_compare=wgpu.CompareFunction.less,
            ),
            multisample=wgpu.MultisampleState(count=1),
            fragment=wgpu.FragmentState(module=fs, entry_point="main", targets=[
                wgpu.ColorTargetState(format=self.color_format)
            ]))

    def _update_depth_texture(self):
        w, h = self.canvas.get_physical_size()
        if not self.depth_texture or self.depth_texture.size != (w, h, 1):
            self.depth_texture = self.device.create_texture(
                label="depth buffer",
                size=(w, h, 1),
                sample_count=1,
                format=wgpu.TextureFormat.depth16unorm,
                usage=wgpu.TextureUsage.RENDER_ATTACHMENT,
            )

            self.depth_view = self.depth_texture.create_view()

    def _update_rt_texture(self):
        w, h = self.canvas.get_physical_size()
        if not self.rt_texture or self.rt_texture.size != (4 * w, h, 1):
            self.rt_texture = self.device.create_texture(
                label="Ray Tracing Texture",
                size=(4 * w, h, 1),
                sample_count=1,
                format=wgpu.TextureFormat.r32float,
                usage=wgpu.TextureUsage.STORAGE_BINDING,
            )

            self.rt_view = self.rt_texture.create_view()

            self._write_ray_tracing_bindings()

    def _create_physics_buffers(self):
        self.positions_staging = self.device.create_buffer(
            label="Positions (Staging Buffer)",
            size=2 * pad_size(self.positions.nbytes),
            usage=wgpu.BufferUsage.COPY_SRC | wgpu.BufferUsage.MAP_WRITE
        )

        self.positions_attr = self.device.create_buffer(
            label="Positions (Attribute)",
            size=self.positions.nbytes,
            usage=wgpu.BufferUsage.COPY_DST | wgpu.BufferUsage.STORAGE | wgpu.BufferUsage.VERTEX
        )

        self.rotations_staging = self.device.create_buffer(
            label="Rotations (Staging Buffer)",
            size=2 * pad_size(self.rotations.nbytes),
            usage=wgpu.BufferUsage.COPY_SRC | wgpu.BufferUsage.MAP_WRITE
        )

        self.rotations_attr = self.device.create_buffer(
            label="Rotations (Attribute)",
            size=self.rotations.nbytes,
            usage=wgpu.BufferUsage.COPY_DST | wgpu.BufferUsage.STORAGE | wgpu.BufferUsage.VERTEX
        )

        self.velocities_staging = self.device.create_buffer(
            label="Velocities (Staging Buffer)",
            size=2 * pad_size(self.velocities.nbytes),
            usage=wgpu.BufferUsage.COPY_SRC | wgpu.BufferUsage.MAP_WRITE
        )

        self.velocities_attr = self.device.create_buffer(
            label="Velocities (Attribute)",
            size=self.velocities.nbytes,
            usage=wgpu.BufferUsage.COPY_DST | wgpu.BufferUsage.STORAGE | wgpu.BufferUsage.VERTEX
        )

        self.angular_velocities_staging = self.device.create_buffer(
            label="Angular Velocities (Staging Buffer)",
            size=2 * pad_size(self.angular_velocities.nbytes),
            usage=wgpu.BufferUsage.COPY_SRC | wgpu.BufferUsage.MAP_WRITE
        )

        self.angular_velocities_attr = self.device.create_buffer(
            label="Angular Velocities (Attribute)",
            size=self.angular_velocities.nbytes,
            usage=wgpu.BufferUsage.COPY_DST | wgpu.BufferUsage.STORAGE | wgpu.BufferUsage.VERTEX
        )

        if len(self.radii) != 0:
            self.radii_buffer = self.device.create_buffer_with_data(
                label="Radii (Attribute)",
                data=self.radii,
                usage=wgpu.BufferUsage.STORAGE | wgpu.BufferUsage.VERTEX
            )
        else:
            self.radii_buffer = self.device.create_buffer(
                label="Radii (Attribute)",
                size=self.radii.nbytes,
                usage=wgpu.BufferUsage.STORAGE | wgpu.BufferUsage.VERTEX
            )

        if len(self.materials) != 0:
            self.materials_buffer = self.device.create_buffer_with_data(
                label="Materials",
                data=self.materials,
                usage=wgpu.BufferUsage.VERTEX | wgpu.BufferUsage.STORAGE
            )
        else:
            self.materials_buffer = self.device.create_buffer(
                label="Materials",
                size=self.materials.nbytes,
                usage=wgpu.BufferUsage.VERTEX | wgpu.BufferUsage.STORAGE
            )

        self._write_ray_tracing_bindings()

    def _load_sobol_directions(self):
        with open("sobol_directions.data", "rb") as sobol_file:
            self.sobol_directions = self.device.create_buffer_with_data(
                label="Sobol Directions",
                data=np.fromfile(sobol_file, '<u4'),
                usage=wgpu.BufferUsage.STORAGE)

    def _create_path_tracing_pipeline(self):
        header = "#version 450\n"
        if self.supports_barriers:
            header += "#define HAS_BARRIER 1\n"
            header += "#define NUM_SAMPLES 1\n"
        else:
            header += "#define NUM_SAMPLES 8\n"

        shader = self.device.create_shader_module(
            label="path_trace_comp",
            code=header + "\n" + \
                 Path("shaders/path_trace_header.glsl").read_text() + "\n" + \
                 Path("shaders/homework.glsl").read_text() + "\n" + \
                 Path("shaders/path_trace.comp").read_text())

        self.ray_tracing_set_layout = self.device.create_bind_group_layout(
            label="Ray Tracing Bind Group",
            entries=[
                wgpu.BindGroupLayoutEntry(binding=0, visibility=wgpu.ShaderStage.COMPUTE,
                                          buffer=wgpu.BufferBindingLayout(
                                              type=wgpu.BufferBindingType.storage)),
                wgpu.BindGroupLayoutEntry(binding=1, visibility=wgpu.ShaderStage.COMPUTE,
                                          buffer=wgpu.BufferBindingLayout(
                                              type=wgpu.BufferBindingType.storage)),
                wgpu.BindGroupLayoutEntry(binding=2, visibility=wgpu.ShaderStage.COMPUTE,
                                          buffer=wgpu.BufferBindingLayout(
                                              type=wgpu.BufferBindingType.storage)),
                wgpu.BindGroupLayoutEntry(binding=3, visibility=wgpu.ShaderStage.COMPUTE,
                                          buffer=wgpu.BufferBindingLayout(
                                              type=wgpu.BufferBindingType.storage)),
                wgpu.BindGroupLayoutEntry(binding=4, visibility=wgpu.ShaderStage.COMPUTE,
                                          buffer=wgpu.BufferBindingLayout(
                                              type=wgpu.BufferBindingType.storage)),
                wgpu.BindGroupLayoutEntry(binding=5, visibility=wgpu.ShaderStage.COMPUTE,
                                          storage_texture=wgpu.StorageTextureBindingLayout(
                                              access=wgpu.StorageTextureAccess.read_write,
                                              format=wgpu.TextureFormat.r32float)),
                wgpu.BindGroupLayoutEntry(binding=6, visibility=wgpu.ShaderStage.COMPUTE,
                                          buffer=wgpu.BufferBindingLayout(
                                              type=wgpu.BufferBindingType.uniform))
            ]
        )

        self.blit_image_layout = self.device.create_bind_group_layout(
            label="Ray Tracing Blit Group",
            entries=[
                wgpu.BindGroupLayoutEntry(binding=0, visibility=wgpu.ShaderStage.FRAGMENT,
                                          storage_texture=wgpu.StorageTextureBindingLayout(
                                              access=wgpu.StorageTextureAccess.read_only,
                                              format=wgpu.TextureFormat.r32float)),
            ]
        )

        self.ray_tracing_layout = self.device.create_pipeline_layout(
            label="Ray Tracing Layout",
            bind_group_layouts=[
                self.uniform_layout,
                self.ray_tracing_set_layout,
                *self.user_data.gpu_bind_group_layouts(self.device)
            ]
        )

        self.blit_layout = self.device.create_pipeline_layout(
            label="Ray Tracing Layout",
            bind_group_layouts=[
                self.blit_image_layout,
            ]
        )

        self.ray_tracing_pipeline = self.device.create_compute_pipeline(
            label="Ray Tracing Pipeline",
            layout=self.ray_tracing_layout,
            compute=wgpu.ProgrammableStage(module=shader, entry_point="main"))

        vs = self.device.create_shader_module(label="Blit_vert", code=Path("shaders/rt_blit.vert").read_text())
        fs = self.device.create_shader_module(label="Blit_frag", code=Path("shaders/rt_blit.frag").read_text())

        self.blit_pipeline = self.device.create_render_pipeline(
            layout=self.blit_layout,
            vertex=wgpu.VertexState(module=vs, entry_point="main", buffers=[]),
            primitive=wgpu.PrimitiveState(),
            multisample=wgpu.MultisampleState(count=1),
            fragment=wgpu.FragmentState(module=fs, entry_point="main", targets=[
                wgpu.ColorTargetState(format=self.color_format)
            ]))

    def _write_ray_tracing_bindings(self):
        if not self.rt_view or self.radii_buffer.size == 0:
            return

        self.ray_tracing_bindings = self.device.create_bind_group(
            label="Ray Tracing Bindings",
            layout=self.ray_tracing_set_layout,
            entries=[
                wgpu.BindGroupEntry(binding=0, resource=wgpu.BufferBinding(buffer=self.sobol_directions)),
                wgpu.BindGroupEntry(binding=1, resource=wgpu.BufferBinding(buffer=self.positions_attr)),
                wgpu.BindGroupEntry(binding=2, resource=wgpu.BufferBinding(buffer=self.rotations_attr)),
                wgpu.BindGroupEntry(binding=3, resource=wgpu.BufferBinding(buffer=self.radii_buffer)),
                wgpu.BindGroupEntry(binding=4, resource=wgpu.BufferBinding(buffer=self.materials_buffer)),
                wgpu.BindGroupEntry(binding=5, resource=self.rt_view),
                wgpu.BindGroupEntry(binding=6, resource=wgpu.BufferBinding(buffer=self.rt_reset_ubo)),
            ])

        self.blit_bindings = self.device.create_bind_group(
            label="Ray Tracing Bindings (Blit)",
            layout=self.blit_image_layout,
            entries=[
                wgpu.BindGroupEntry(binding=0, resource=self.rt_view)
            ])
