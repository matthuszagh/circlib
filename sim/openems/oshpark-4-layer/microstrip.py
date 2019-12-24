#!/usr/bin/env python

import os, tempfile
from CSXCAD import CSXCAD as csxcad
import openEMS as openems
from openEMS.physical_constants import EPS0, C0
import numpy as np
import matplotlib.pyplot as plt
from automesh import AutoMesh


# ============================ parameters ============================
# sim_path = os.path.join(tempfile.gettempdir(), "sim_microstrip")
# sim_path = os.getcwd()
sim_path = os.path.abspath("sim_microstrip")

# dimensions
unit = 1e-3
microstrip_width = 0.38
microstrip_len = 50
outer_cu_thick = 0.03556
inner_cu_thick = 0.01778
substrate_len = 1.2 * microstrip_len
substrate_width = substrate_len / 2
substrate_thick = 1.6

# frequency
f0 = 5.6e9
fc = 4e8  # 20dB corner frequency
lmin = C0 / ((f0 + fc) * unit)

# substrate parameters
substrate_epsr = 3.69
substrate_kappa = unit * 2 * np.pi * 2.45e9 * EPS0 * substrate_epsr

# # feed
# feed_pos_x = -microstrip_len / 2
# feed_res = 50

# ============================ primitives ============================
csx = csxcad.ContinuousStructure()
# simulation box [x_min, x_max, y_min, y_max, z_min, z_max]
simbox = [
    [-1.2 * substrate_len / 2, 1.2 * substrate_len / 2],
    [-1.2 * substrate_width, 1.2 * substrate_width],
    [-5 * substrate_thick, 10 * substrate_thick],
]
# microstrip
microstrip = csx.AddMetal("PEC")
# microstrip = csx.AddMaterial("Cu", kappa=56e6)
start = [-microstrip_len / 2, -microstrip_width / 2, 0]
stop = [microstrip_len / 2, microstrip_width / 2, 0]
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
start = [-substrate_len / 2, -substrate_width / 2, -substrate_thick]
stop = [substrate_len / 2, substrate_width / 2, -substrate_thick]
# gnd.AddBox(priority=10, start=start, stop=stop)

# =============================== mesh ===============================
# mesh = csx.GetGrid()
# mesh.SetDeltaUnit(unit)
auto_mesh = AutoMesh(
    csx,
    lmin,
    simbox,
    mres=1 / 20,
    sres=1 / 10,
    smooth=1.5,
    unit=1e-3,
    min_lines=10,
)
auto_mesh.AutoGenMesh()
# cres = lmin / 20  # conductor mesh resolution
# sres = lmin / 10  # substrate mesh resolution
# ares = lmin  # air mesh resolution
# # initialize mesh with outer boundaries of simulation box dimensions
# mesh.AddLine("x", [simbox[0], simbox[1]])
# mesh.AddLine("y", [simbox[2], simbox[3]])
# mesh.AddLine("z", [simbox[4], simbox[5]])
# # microstrip mesh
# mesh.AddLine(
#     "x",
#     np.linspace(
#         -microstrip_len / 2,
#         microstrip_len / 2,
#         int(np.ceil(microstrip_len / cres)),
#     ),
# )
# micro_sep = microstrip_width / (9 + (2 / 3))
# # mesh.AddLine(
# #     "y",
# #     [
# #         microstrip_width / 2 - micro_sep / 3,
# #         -microstrip_width / 2 + micro_sep / 3,
# #     ],
# # )
# # mesh.AddLine(
# #     "y",
# #     [
# #         [
# #             microstrip_width / 2 - micro_sep / 3 - (i * micro_sep)
# #             for i in range(1, 4)
# #         ],
# #         [
# #             -microstrip_width / 2 + micro_sep / 3 + (i * micro_sep)
# #             for i in range(1, 4)
# #         ],
# #     ],
# # )
# # mesh.AddLine(
# #     "z", np.linspace(0, outer_cu_thick, int(np.ceil(outer_cu_thick / cres))),
# # )
# mesh.AddLine("z", [0])

# # substrate mesh
# mesh.AddLine(
#     "y",
#     np.linspace(
#         -substrate_width / 2,
#         -microstrip_width / 2,
#         int(np.ceil((substrate_width / 2 - (microstrip_width / 2)) / sres)),
#     ),
# )
# mesh.AddLine(
#     "y",
#     np.linspace(
#         microstrip_width / 2,
#         substrate_width / 2,
#         int(np.ceil((substrate_width / 2 - (microstrip_width / 2)) / sres)),
#     ),
# )
# mesh.AddLine(
#     "z",
#     np.linspace(
#         -substrate_thick, 0, int(10 * np.ceil(substrate_thick / sres)),
#     ),
# )

# # air mesh
# # mesh.AddLine(
# #     "x",
# #     np.linspace(
# #         simbox[0],
# #         -substrate_len / 2,
# #         int(np.ceil((-substrate_len / 2 - simbox[0]) / ares)),
# #     ),
# # )
# mesh.AddLine(
#     "y",
#     np.linspace(
#         simbox[2],
#         -substrate_width / 2,
#         int(np.ceil((-substrate_width / 2 - simbox[2]) / ares)),
#     ),
# )
# mesh.AddLine(
#     "y",
#     np.linspace(
#         substrate_width / 2,
#         simbox[3],
#         int(np.ceil((simbox[3] - (substrate_width / 2)) / ares)),
#     ),
# )
# mesh.AddLine(
#     "z",
#     np.linspace(
#         outer_cu_thick,
#         simbox[5],
#         int(np.ceil((simbox[5] - outer_cu_thick) / ares)),
#     ),
# )

# mesh.SmoothMeshLines("all", max_res=cres, ratio=1.5)

# ============================ simulation ============================
fdtd = openems.openEMS()
fdtd.SetCSX(csx)
fdtd.SetGaussExcite(f0, fc)
fdtd.SetBoundaryCond(["PML_8", "PML_8", "MUR", "MUR", "PEC", "MUR"])
# # add ports
# ports = [None, None]
# start = [0, microstrip_width / 2, 0]
# stop = [-microstrip_len / 2, -microstrip_width / 2, -substrate_thick]
# ports[0] = fdtd.AddMSLPort(
#     port_nr=1,
#     metal_prop=microstrip,
#     start=start,
#     stop=stop,
#     prop_dir="x",
#     exc_dir="z",
#     excite=-1,
#     FeedShift=10 * cres,  # TODO what is this
#     MeasPlaneShift=microstrip_len / 3,  # TODO what is this
#     priority=10,
# )
# start = [microstrip_len / 2, microstrip_width / 2, 0]
# stop = [0, -microstrip_width / 2, -substrate_thick]
# ports[1] = fdtd.AddMSLPort(
#     port_nr=2,
#     metal_prop=microstrip,
#     start=start,
#     stop=stop,
#     prop_dir="x",
#     exc_dir="z",
#     excite=0,
#     MeasPlaneShift=microstrip_len / 3,  # TODO what is this
#     priority=10,
# )
# # add excitation
# start = [-microstrip_len / 2, -microstrip_width / 2, -substrate_thick]
# stop = [-microstrip_len / 2, microstrip_width / 2, 0]
# port = fdtd.AddLumpedPort(
#     1, feed_res, start, stop, "z", 1, priority=5, edges2grid="xy"
# )

# E-field recording
Et = csx.AddDump(os.path.join(sim_path, "Et_"), file_type=0)
start = [-substrate_len / 2, -substrate_width / 2, -substrate_thick / 2]
stop = [substrate_len / 2, substrate_width / 2, -substrate_thick / 2]
Et.AddBox(start, stop)

csx_file = os.path.join(sim_path, "csx.xml")
csx.Write2XML(csx_file)
# os.system(r"AppCSXCAD {}".format(csx_file))

# fdtd.Run(sim_path, verbose=3, cleanup=True)

# # ========================== post-processing =========================
# freq = np.linspace(f0 - fc, f0 + fc, 401)
# # excitation port
# for p in ports:
#     p.CalcPort(sim_path, freq, ref_impedance=50)

# # S11/S21
# s11 = ports[0].uf_ref / ports[0].uf_inc
# s11_db = 20 * np.log10(np.abs(s11))
# s21 = ports[1].uf_ref / ports[0].uf_inc
# s21_db = 20 * np.log10(np.abs(s21))
# plt.figure()
# plt.plot(freq / 1e9, s11_db, "k-", linewidth=2, label="$S_{11}$")
# plt.plot(freq / 1e9, s21_db, "r--", linewidth=2, label="$S_{21}$")
# plt.grid()
# plt.legend()
# plt.ylabel("S-Parameter (dB)")
# plt.xlabel("Frequency (GHz)")
# # # Impedance
# # zin = port.uf_tot / port.if_tot
# # plt.figure()
# # plt.plot(freq / 1e9, np.real(zin), "k-", linewidth=2, label="$\Re{Z_{in}}$")
# # plt.plot(freq / 1e9, np.imag(zin), "r--", linewidth=2, label="$\Im{Z_{in}}$")
# # plt.grid()
# # plt.title("Excitation port")
# # plt.legend()
# # plt.ylabel("Zin (Ohm)")
# # plt.xlabel("Frequency (GHz)")

# plt.show()
