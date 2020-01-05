#!/usr/bin/env python

from params import sim_path, f0, fc
import numpy as np
import matplotlib.pyplot as plt
from openEMS.ports import Port

freq = np.linspace(f0 - fc, f0 + fc, 2000)

# ===================== characteristic impedance =====================
plt.figure()
z0_ports = [None, None, None]
for n, port in enumerate(z0_ports):
    port = Port(
        None,
        None,
        None,
        None,
        None,
        U_filenames=["vprobe_" + str(n)],
        I_filenames=["iprobe_" + str(n)],
    )
    port.CalcPort(sim_path, freq, ref_impedance=50)
    plt.plot(
        freq / 1e9,
        np.absolute(port.uf_tot) / np.absolute(port.if_tot),
        linewidth=2,
        label="Port " + str(n),
    )

plt.legend()
plt.grid()
plt.ylabel("Characteristic Impedance ($Z_0$)")
plt.xlabel("Frequency (GHz)")

# =========================== S-parameters ===========================
plt.figure()
s_ports = [None, None]
for n, _ in enumerate(s_ports):
    s_ports[n] = Port(
        None,
        None,
        None,
        None,
        None,
        U_filenames=["port_ut_" + str(n)],
        I_filenames=["port_it_" + str(n)],
    )
    s_ports[n].CalcPort(sim_path, freq, ref_impedance=50)

s11 = s_ports[0].uf_ref / s_ports[0].uf_inc
# s11 = s_ports[1].uf_inc / s_ports[1].uf_ref
s11_db = 20 * np.log10(np.abs(s11))
# s21 = s_ports[1].uf_inc / s_ports[0].uf_inc
s21 = s_ports[1].uf_ref / s_ports[0].uf_inc
s21_db = 20 * np.log10(np.abs(s21))
plt.plot(freq / 1e9, s11_db, linewidth=2, label="$S_{11}$")
plt.plot(freq / 1e9, s21_db, linewidth=2, label="$S_{21}$")

plt.legend()
plt.grid()
plt.ylabel("S (dB)")
plt.xlabel("Frequency (GHz)")

# =============================== power ==============================
plt.figure()
plt.plot(freq / 1e9, s_ports[1].P_acc, linewidth=2, label="P_acc")
plt.legend()
plt.grid()
plt.ylabel("Power Accepted (dB)")
plt.xlabel("Frequency (GHz)")

plt.show()
