import os
from openEMS.physical_constants import C0


sim_path = os.path.abspath("sim")

# dimensions
unit = 1e-3
microstrip_width = 0.4
# microstrip_width = 0.3429
microstrip_len = 50
# outer_cu_thick = 0.0356
# inner_cu_thick = 0.0178
outer_cu_thick = 0
inner_cu_thick = 0
substrate_len = 1.2 * microstrip_len
substrate_width = substrate_len / 2
substrate_thick = 0.1702

# frequency
f0 = 5.6e9
fc = f0
# fc = 4e8  # 20dB corner frequency
lmin = C0 / ((f0 + fc) * unit)

# substrate parameters
substrate_epsr = 3.64
substrate_resistivity = 4.4e14  # in mm
# substrate_kappa = unit * 2 * np.pi * 2.45e9 * EPS0 * substrate_epsr
substrate_kappa = 1 / substrate_resistivity

# csx_file = os.path.join(sim_path, "csx.xml")
csx_file = "csx.xml"
