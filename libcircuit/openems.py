from bisect import bisect_left
from automesh import Mesh
from pathos.multiprocessing import ProcessingPool as Pool
from functools import partial
from subprocess import run
import tempfile
import os
from shutil import rmtree
from typing import List, Tuple
from CSXCAD import CSXCAD as csxcad
import openEMS as openems
from openEMS.physical_constants import C0
from openEMS.ports import UI_data
import numpy as np


class PCB:
    """
    A PCB structure and material properties for use in an openems simulation.
    """

    def __init__(
        self,
        layers: int,
        sub_epsr: List[Tuple[float, float]],
        sub_rho: float,
        layer_sep: List[float],
        layer_thickness: List[float],
    ):
        """
        :param layers: number of conductive layers
        :param sub_epsr: substrate dielectric constant.  dictionary of
            frequency (Hz) and associated dielectric.
        :param sub_rho: volume resistivity (ohm*mm)
        :param layer_sep: separations (in mm) between adjacent copper
            layers.  A list where the first value is the separation
            between the top layer and second layer, etc.  This is
            equivalently the substrate thickness.
        :param layer_thickness: thickness of each conductive layer (in
            mm).  Again proceeds from top to bottom layer.
        """
        self.layers = layers
        self.sub_epsr = sorted(sub_epsr, key=lambda epsr: epsr[0])
        self.sub_rho = sub_rho
        self.sub_kappa = 1 / sub_rho
        self.layer_sep = layer_sep
        self.layer_thickness = layer_thickness

    def epsr_at_freq(self, freq: float):
        """
        Approximate the dielectric at a given frequency given the
        provided epsr values.

        :param freq: frequency of interest (Hz)

        :param returns: dielectric constant
        """
        if freq <= self.sub_epsr[0][0]:
            return self.sub_epsr[0][1]
        elif freq >= self.sub_epsr[-1][0]:
            return self.sub_epsr[-1][1]

        # perform linear interpolation
        tup_low = bisect_left([x[0] for x in self.sub_epsr], freq)
        if self.sub_epsr[tup_low][0] == freq:
            return self.sub_epsr[tup_low][1]

        tup_high = tup_low
        tup_low -= 1
        xlow = self.sub_epsr[tup_low][0]
        xhigh = self.sub_epsr[tup_high][0]
        ylow = self.sub_epsr[tup_low][1]
        yhigh = self.sub_epsr[tup_high][1]
        slope = (yhigh - ylow) / (xhigh - xlow)
        return ylow + (slope * (freq - xlow))


common_pcbs = {
    "oshpark4": PCB(
        layers=4,
        sub_epsr=[
            (100e6, 3.72),
            (1e9, 3.69),
            (2e9, 3.68),
            (5e9, 3.64),
            (10e9, 3.65),
        ],
        sub_rho=4.4e14,
        layer_sep=[0.1702, 1.1938, 0.1702],
        layer_thickness=[0.0356, 0.0178, 0.0178, 0.0356],
    )
}


def wheeler_z0(w: float, t: float, er: float, h: float) -> float:
    """
    Calculate the microstrip characteristic impedance for a given
    width using Wheeler's equation.  Wheeler's equation can be found
    at:

    https://en.wikipedia.org/wiki/Microstrip#Characteristic_impedance

    :param w: microstrip trace width (mm)
    :param t: trace thickness (mm)
    :param er: substrate relative permittivity
    :param h: substrate height (thickness) (mm)

    :returns: characteristic impedance
    """
    z0 = 376.730313668
    w *= 1e-3
    t *= 1e-3
    h *= 1e-3
    weff = w + (
        (t * ((1 + (1 / er)) / (2 * np.pi)))
        * np.log(
            (4 * np.e)
            / (
                np.sqrt(
                    ((t / h) ** 2)
                    + (((1 / np.pi) * ((1 / ((w / t) + (11 / 10))))) ** 2)
                )
            )
        )
    )
    tmp1 = 4 * h / weff
    tmp2 = (14 + (8 / er)) / 11
    zm = (z0 / (2 * np.pi * np.sqrt(2 * (1 + er)))) * np.log(
        1
        + (
            tmp1
            * (
                (tmp2 * tmp1)
                + (
                    np.sqrt(
                        (
                            (tmp2 * tmp1) ** 2
                            + ((np.pi ** 2) * ((1 + (1 / er)) / 2))
                        )
                    )
                )
            )
        )
    )
    return zm


def wheeler_z0_width(
    z0: float,
    t: float,
    er: float,
    h: float,
    tol: float = 0.01,
    guess: float = 0.3,
) -> float:
    """
    Calculate the microstrip width for a given characteristic
    impedance using Wheeler's formula.

    :param z0: characteristic impedance (ohm)
    :param t: trace thickness (mm)
    :param er: substrate relative permittivity
    :param h: substrate height (thickness) (mm)
    :param tol: acceptable impedance tolerance (ohm)
    :param guess: an initial guess for the width (mm).  This can
        improve convergence time when the approximate width is known.

    :returns: trace width (mm)
    """
    width = guess
    zm = wheeler_z0(w=width, t=t, er=er, h=h)
    wlow = width / 10
    zlow = wheeler_z0(w=wlow, t=t, er=er, h=h)
    # inverse relation between width and z0
    while zlow < z0:
        wlow /= 10
        zlow = wheeler_z0(w=wlow, t=t, er=er, h=h)

    whigh = width * 10
    zhigh = wheeler_z0(w=whigh, t=t, er=er, h=h)
    while zhigh > z0:
        whigh *= 10
        zhigh = wheeler_z0(w=whigh, t=t, er=er, h=h)

    while np.absolute(zm - z0) > tol:
        if zm > z0:
            m = (zhigh - zm) / (whigh - width)
            wlow = width
            zlow = zm
        else:
            m = (zm - zlow) / (width - wlow)
            whigh = width
            zhigh = zm

        # use linear interpolation to update guess
        width = width + ((z0 - zm) / m)
        zm = wheeler_z0(w=width, t=t, er=er, h=h)

    return width


class Probe:
    """
    Wrapper class for openems probes.  The main additional feature is
    this class holds onto data.
    """

    def __init__(
        self, name: str, csx: csxcad.ContinuousStructure(), box, p_type=0,
    ):
        """
        :param name: probe name, which becomes filename where data is
            stored.
        :param csx: csx structure to add probe to.
        :param box: 2D list, outer dim=2, inner dim=3.  Outer
            dimensions is start and stop of box, and inner dimension
            are the x, y, and z coordinates, respectively.
        :param p_type: 0 if voltage probe, 1 if current probe.
        """
        self.data = None
        self.f_data = None
        self.box = box
        self.csx = csx
        self.name = name
        self.p_type = p_type
        self.csxProbe = self.csx.AddProbe(self.name, self.p_type)
        self.csxProbe.AddBox(start=self.box[0], stop=self.box[1])

    def calc_frequency_data(self, sim_path, freq, signal_type="pulse"):
        self.data = UI_data([self.name], sim_path, freq, signal_type)
        self.f_data = self.data.ui_f_val[0]


class Microstrip:
    """
    Microstrip transmission line.
    """

    def __init__(
        self,
        pcb: PCB,
        f0: float,
        fc: float,
        z0_ref: float = 50,
        microstrip_width: float = None,
        microstrip_len: float = 100,
        substrate_width: float = 20,
        fcsx: str = None,
        fvtr_dir: str = None,
        efield: bool = False,
    ):
        """
        :param pcb: PCB object
        :param f0: center frequency (Hz)
        :param fc: corner frequency, given as difference from f0 (Hz)
        :param z0_ref: desired impedance (ohm)
        :param microstrip_width: microstrip trace width (mm).  If
            ommitted, an analytical best guess for z0_ref and f0 will
            be used (see `wheeler_z0_width`).
        :param microstrip_len: microstrip trace length (mm).  This
            value isn't critical but does need to be large enough for
            the openems simulation to work.  If unsure, stick with the
            default.
        :param substrate_width: substrate width (mm)
        :param substrate_len: substrate length (mm).  Must be at least
            as large as the microstrip length.  If unsure, use the
            default.
        :param fcsx: CSX file name to write CSXCAD structure to.  If
            ommitted, do not save the file.  Omitting this does not
            prevent you from viewing the CSX structure, just from
            saving it for later.  Relative to current directory.
        :param fvtr_dir: directory for VTR E-field dump files.  Files
            will be prefixed with 'Et_'.  If ommitted, E-field can
            still be viewed unless `efield` is set to false in which
            case this parameter has no effect.
        :param efield: dump E-field time values for viewing.
        """
        self.pcb = pcb
        self.f0 = f0
        self.fc = fc
        self.z0_ref = z0_ref
        if microstrip_width is None:
            self.microstrip_width = wheeler_z0_width(
                z0=z0_ref,
                t=self.pcb.layer_thickness[0],
                er=self.pcb.epsr_at_freq(f0),
                h=self.pcb.layer_sep[0],
            )
        else:
            self.microstrip_width = microstrip_width
        self.microstrip_len = microstrip_len
        self.substrate_width = substrate_width
        if fcsx is None:
            self.fcsx = tempfile.mkstemp()[1]
        else:
            self.fcsx = os.path.abspath(fcsx)

        self.efield = efield
        if fvtr_dir is None:
            if efield:
                self.fvtr_dir = tempfile.mkdtemp()
        else:
            self.fvtr_dir = os.path.abspath(fvtr_dir)
            if os.path.exists(self.fvtr_dir):
                rmtree(self.fvtr_dir)
            os.mkdir(self.fvtr_dir)

        self.mesh = None
        self.csx = None
        self.vprobes = [None]
        self.iprobes = [None]
        # track simulation steps so that we don't reperform steps
        self.sim_done = False
        self.csx_done = False
        self.freq_bins = None
        self.fdtd = None

    def sim(
        self, num_freq_bins: int = 501, zero_trace_height: bool = False,
    ):
        """
        Run an openems simulation for the current microstrip structure
        and write the results to a file.

        :param num_freq_bins: number of frequency bins.  More
            frequency bins means longer simulation time but possibly
            greater accuracy.
        :param zero_trace_height: approximate trace height as 0.  This
            can greatly reduce simulation time but may reduce
            accuracy.
        """
        self.gen_csx(zero_trace_height=zero_trace_height)
        tmpdir = tempfile.mkdtemp()
        self.fdtd.Run(tmpdir, cleanup=True)

        self.freq_bins = np.linspace(
            self.f0 - self.fc, self.f0 + self.fc, num_freq_bins
        )

        for i, _ in enumerate(self.vprobes):
            self.vprobes[i].calc_frequency_data(
                sim_path=tmpdir, freq=self.freq_bins
            )

        for i, _ in enumerate(self.iprobes):
            self.iprobes[i].calc_frequency_data(
                sim_path=tmpdir, freq=self.freq_bins
            )

        self.calc_params()
        self.sim_done = True
        return np.absolute(self.Z_ref)

    def gen_csx(self, zero_trace_height: bool = False) -> None:
        """
        Generate CSX structure.

        :param zero_trace_height: approximate trace height as 0.  This
            can greatly reduce simulation time but may reduce
            accuracy.  Args: num_probes: zero_trace_height:
        """
        # dimensions and parameters
        unit = 1e-3
        lmin = C0 / ((self.f0 + self.fc) * unit)
        if zero_trace_height:
            trace_height = 0
        else:
            trace_height = self.pcb.layer_thickness[0]

        # structures
        self.csx = csxcad.ContinuousStructure()
        microstrip = self.csx.AddMetal("PEC")
        microstrip.AddBox(
            priority=10,
            start=[-self.microstrip_len / 2, -self.microstrip_width / 2, 0],
            stop=[
                self.microstrip_len / 2,
                self.microstrip_width / 2,
                trace_height,
            ],
        )
        substrate = self.csx.AddMaterial(
            "substrate",
            epsilon=self.pcb.epsr_at_freq(self.f0),
            kappa=self.pcb.sub_kappa,
        )
        substrate.AddBox(
            priority=0,
            start=[
                -self.microstrip_len / 2,
                -self.substrate_width / 2,
                -self.pcb.layer_sep[0],
            ],
            stop=[self.microstrip_len / 2, self.substrate_width / 2, 0],
        )

        # simulation
        self.fdtd = openems.openEMS(EndCriteria=1e-5)
        self.fdtd.SetCSX(self.csx)
        self.fdtd.SetGaussExcite(self.f0, self.fc)
        self.fdtd.SetBoundaryCond(
            ["PML_8", "PML_8", "MUR", "MUR", "PEC", "MUR"]
        )

        self.mesh = Mesh(
            self.csx,
            lmin,
            mres=1 / 20,
            sres=1 / 10,
            smooth=1.4,
            unit=unit,
            min_lines=9,
            expand_bounds=[0, 0, 10, 10, 0, 10],
        )
        self.mesh.generate_mesh()

        _, port_xpos = self.mesh.nearest_mesh_line(0, -self.microstrip_len / 4)
        self.fdtd.AddLumpedPort(
            port_nr=0,
            R=self.z0_ref,
            start=[
                port_xpos,
                -self.microstrip_width / 2,
                -self.pcb.layer_sep[0],
            ],
            stop=[port_xpos, self.microstrip_width / 2, 0],
            p_dir="z",
            excite=1,
            priority=999,
        )

        self.gen_probes(trace_height=trace_height)

        # E-field recording
        if self.efield:
            Et = self.csx.AddDump(
                os.path.join(self.fvtr_dir, "Et_"), file_type=0
            )
            start = [
                -self.microstrip_len / 2,
                -self.substrate_width / 2,
                -self.pcb.layer_sep[0] / 2,
            ]
            stop = [
                self.microstrip_len / 2,
                self.substrate_width / 2,
                -self.pcb.layer_sep[0] / 2,
            ]
            Et.AddBox(start, stop)

        # write CSX file
        self.csx.Write2XML(self.fcsx)
        self.csx_done = True

    def gen_probes(self, trace_height):
        mid_idx, mid_xpos = self.mesh.nearest_mesh_line(0, 0)
        vprobe_x_pos = [
            self.mesh.mesh_lines[0][mid_idx - 1],
            mid_xpos,
            self.mesh.mesh_lines[0][mid_idx + 1],
        ]
        iprobe_x_pos = [
            (vprobe_x_pos[0] + vprobe_x_pos[1]) / 2,
            (vprobe_x_pos[1] + vprobe_x_pos[2]) / 2,
        ]

        self.vprobes = [None for i in range(len(vprobe_x_pos))]
        self.iprobes = [None for i in range(len(iprobe_x_pos))]

        for i, _ in enumerate(self.vprobes):
            self.vprobes[i] = Probe(
                name="ut_" + str(i),
                csx=self.csx,
                box=[
                    [vprobe_x_pos[i], 0, -self.pcb.layer_sep[0]],
                    [vprobe_x_pos[i], 0, 0],
                ],
                p_type=0,
            )

        for i, _ in enumerate(self.iprobes):
            self.iprobes[i] = Probe(
                name="it_" + str(i),
                csx=self.csx,
                box=[
                    [iprobe_x_pos[i], -self.microstrip_width / 2, 0],
                    [
                        iprobe_x_pos[i],
                        self.microstrip_width / 2,
                        trace_height,
                    ],
                ],
                p_type=1,
            )

    def calc_params(self):
        Et = self.vprobes[1].f_data
        Ht = 0.5 * (self.iprobes[0].f_data + self.iprobes[1].f_data)
        unit = 1e-3
        dEt = (self.vprobes[2].f_data - self.vprobes[0].f_data) / (
            unit * (self.vprobes[2].box[0][0] - self.vprobes[0].box[0][0])
        )
        dHt = (self.iprobes[1].f_data - self.iprobes[0].f_data) / (
            unit * (self.iprobes[1].box[0][0] - self.iprobes[0].box[0][0])
        )
        beta = np.sqrt(-dEt * dHt / (Ht * Et))
        beta[
            np.real(beta) < 0
        ] *= -1  # determine correct sign (unlike the paper)
        self.beta = beta

        # determine ZL
        self.Z_ref = np.sqrt(Et * dEt / (Ht * dHt))

    def view_csx(self):
        """
        View CSX structure. This blocks the calling thread and AppCSXCAD
        must be closed before proceeding.
        """
        if not self.csx_done:
            raise RuntimeWarning(
                "Must generate CSX before viewing it. Doing nothing."
            )
        else:
            run(["AppCSXCAD", self.fcsx])

    def view_efield(self):
        """
        View E-field time dump. This blocks the calling thread and
        Paraview must be closed before proceeding.
        """
        if not self.sim_done:
            raise RuntimeWarning(
                "Must run simulation before viewing E-field dump. "
                "Doing nothing."
            )
        else:
            run(
                [
                    "paraview",
                    "--data={}".format(
                        os.path.join(self.fvtr_dir, "Et__..vtr")
                    ),
                ]
            )


def microstrip_sweep_width(
    pcb: PCB,
    f0: float,
    fc: float,
    z0_ref: float,
    width: float = None,
    width_dev_factor: float = 0.1,
    num_points: int = 11,
) -> List[Tuple[float, float]]:
    """
    Calculate microstrip characteristic impedance as a function of
    trace width.

    :param pcb: PCB object
    :param f0: center frequency (Hz)
    :param fc: corner frequency, given as difference from f0 (Hz)
    :param z0_ref: desired impedance (ohm)
    :param width: center width (mm).  If ommitted, find the best guess
        width for f0 analytically.
    :param width_dev_factor: determines width sweep bounds
        [width*(1-wdf),width*(1+wdf)]
    :param num_points: number of width values to calculate.  This
        value has a significant effect on simulation time since each
        point is calculated in its own thread.  Therefore, to minimize
        computation time, its recommended to choose some multiple of
        the number of cores available on the simulation machine.  The
        total simulation time will be roughly equal to the time it
        takes to compute 1 point times the ratio of the number of
        points to number of machine cores.
    :param fout: filename where results should be written.  relative
        to current dir.

    :returns: A list of tuples, where the first tuple element is a
        width and the second is the corresponding impedance value.
    """
    if width is None:
        width = wheeler_z0_width(
            z0=z0_ref,
            t=pcb.layer_thickness[0],
            er=pcb.epsr_at_freq(f0),
            h=pcb.layer_sep[0],
        )

    widths = np.linspace(
        width * (1 - width_dev_factor),
        width * (1 + width_dev_factor),
        num_points,
    )
    microstrips = [None for i in range(num_points)]
    for i, width in enumerate(widths):
        microstrips[i] = Microstrip(
            pcb=pcb, f0=f0, fc=fc, z0_ref=z0_ref, microstrip_width=width
        )

    pool = Pool(nodes=11)
    freq_bins = 501
    func = partial(
        Microstrip.sim, num_freq_bins=freq_bins, zero_trace_height=False,
    )
    z0s = [None for i in range(num_points)]
    z0s = list(pool.map(func, microstrips))

    ret_vals = [[None, None] for i in range(num_points)]
    for i, _ in enumerate(z0s):
        f0_bin_idx = int(freq_bins / 2)
        z0 = z0s[i]
        ret_vals[i][0] = widths[i]
        ret_vals[i][1] = z0[f0_bin_idx]

    return ret_vals
