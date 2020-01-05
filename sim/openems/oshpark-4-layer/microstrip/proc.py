#!/usr/bin/env python

from CSXCAD import CSXCAD as csxcad
import openEMS as openems
from automesh import AutoMesh
from params import *


# ============================= structure ============================
csx = csxcad.ContinuousStructure()
# simulation box [x_min, x_max, y_min, y_max, z_min, z_max]
simbox = [
    [-2 * substrate_len / 2, 2 * substrate_len / 2],
    [-2.5 * substrate_width / 2, 2.5 * substrate_width / 2],
    [-50 * substrate_thick, 100 * substrate_thick],
]
# microstrip
microstrip = csx.AddMetal("PEC")
# microstrip = csx.AddMaterial("Cu", kappa=56e6)
start = [-microstrip_len / 2, -microstrip_width / 2, 0]
stop = [microstrip_len / 2, microstrip_width / 2, 0 + outer_cu_thick]
microstrip.AddBox(priority=10, start=start, stop=stop)
# substrate
substrate = csx.AddMaterial(
    "FR408HR", epsilon=substrate_epsr, kappa=substrate_kappa
)
start = [-substrate_len / 2, -substrate_width / 2, -substrate_thick]
stop = [substrate_len / 2, substrate_width / 2, 0]
substrate.AddBox(priority=0, start=start, stop=stop)
# ground
gnd = csx.AddMetal("PEC")
start = [
    -substrate_len / 2,
    -substrate_width / 2,
    -substrate_thick - inner_cu_thick,
]
stop = [substrate_len / 2, substrate_width / 2, -substrate_thick]
gnd.AddBox(priority=10, start=start, stop=stop)

# ============================ simulation ============================
fdtd = openems.openEMS()
fdtd.SetCSX(csx)
fdtd.SetGaussExcite(f0, fc)
fdtd.SetBoundaryCond(["PML_8", "PML_8", "PML_8", "PML_8", "PML_8", "PML_8"])
# add ports and excitation
ports = [None, None]
start = [
    # -microstrip_len / 2,
    -microstrip_len / 2 + (microstrip_len / 200),
    0,
    # -microstrip_width / 2,
    -substrate_thick,
]
stop = [
    # -microstrip_len / 2,
    -microstrip_len / 2 + (microstrip_len / 200),
    0,
    # microstrip_width / 2,
    0,
]
ports[0] = fdtd.AddLumpedPort(
    port_nr=0, R=50, start=start, stop=stop, p_dir="z", excite=1, priority=999,
)

start = [
    microstrip_len / 2 - (microstrip_len / 200),
    # microstrip_len / 2,
    0,
    # -microstrip_width / 2,
    -substrate_thick,
]
stop = [
    microstrip_len / 2 - (microstrip_len / 200),
    # microstrip_len / 2,
    0,
    # microstrip_width / 2,
    0,
]
ports[1] = fdtd.AddLumpedPort(
    port_nr=1, R=50, start=start, stop=stop, p_dir="z", excite=0, priority=999,
)

# start = [0, 0, -substrate_thick]
# stop = [0, 0, 0]
# vprobe = csx.AddProbe("vprobe", 0)
# vprobe.AddBox(start=start, stop=stop)
# start = [0, -microstrip_width / 2, 0]
# stop = [0, microstrip_width / 2, 0 + outer_cu_thick]
# iprobe = csx.AddProbe("iprobe", 1, norm_dir=0)
# iprobe.AddBox(start=start, stop=stop)

vprobe = [None, None, None]
iprobe = [None, None, None]

for i, probe in enumerate(vprobe):
    start = [-microstrip_len / 2 + i * microstrip_len / 4, 0, -substrate_thick]
    stop = [-microstrip_len / 2 + i * microstrip_len / 4, 0, 0]
    vprobe[i] = csx.AddProbe("vprobe_" + str(i), 0)
    vprobe[i].AddBox(start=start, stop=stop)

for i, probe in enumerate(iprobe):
    start = [
        -microstrip_len / 2 + ((i + 1) * microstrip_len / 4),
        -microstrip_width / 2,
        0,
    ]
    stop = [
        -microstrip_len / 2 + ((i + 1) * microstrip_len / 4),
        microstrip_width / 2,
        0 + outer_cu_thick,
    ]
    iprobe[i] = csx.AddProbe("iprobe_" + str(i), 1, norm_dir=0)
    iprobe[i].AddBox(start=start, stop=stop)

# =============================== mesh ===============================
auto_mesh = AutoMesh(
    csx,
    lmin,
    simbox,
    mres=1 / 20,
    sres=1 / 10,
    smooth=1.4,
    unit=1e-3,
    min_lines=5,
)
auto_mesh.AutoGenMesh()

# E-field recording
Et = csx.AddDump(os.path.join(sim_path, "Et_"), file_type=0)
start = [-substrate_len / 2, -substrate_width / 2, -substrate_thick / 2]
stop = [substrate_len / 2, substrate_width / 2, -substrate_thick / 2]
Et.AddBox(start, stop)

csx.Write2XML(csx_file)

fdtd.Run(sim_path, verbose=3, cleanup=True)
