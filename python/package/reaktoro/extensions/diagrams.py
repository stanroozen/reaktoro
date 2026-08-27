"""
reaktoro.extensions.diagrams
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CHNOSZ-inspired thermodynamic diagram utilities for Reaktoro.

Classes
-------
PredominancePlot
    Eh–pH (or any 2-D grid) predominance / stability-field diagram.

ActivityDiagram
    Stability-field diagram on log(aX) vs log(aY) axes.

SolubilityPlot
    Filled-contour diagram of total dissolved element concentration.

MosaicPlot
    Layered predominance diagram for multiple species groups (e.g.,
    aqueous + minerals).

SpeciationPlot
    1-D speciation plot of log activity or fractional abundance vs a sweep
    variable.

Helper functions
----------------
format_species_label(formula)
    Convert a Reaktoro species formula string to a matplotlib mathtext label.

axis_label(variable)
    Return a formatted axis-label string for a recognised variable name.

water_lines(ax, T_K=298.15, P_Pa=1e5, pH_range=None, **line_kw)
    Draw the upper O₂ / lower H₂ water-stability boundaries on *ax*.

saturation_curve(ax, xvalues, yvalues, si_grid, label=None, **line_kw)
    Draw the SI = 0 iso-saturation line of a mineral on *ax*.
"""

from __future__ import annotations

import re
import warnings
from typing import List, Optional, Sequence, Union

import numpy as np

# ---------------------------------------------------------------------------
# Optional matplotlib import – kept lazy so the module can be imported in
# environments that lack matplotlib (e.g. headless CI) without crashing.
# ---------------------------------------------------------------------------
try:
    import matplotlib as _mpl
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.patches as mpatches
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    _HAS_MPL = True
except ImportError:  # pragma: no cover
    _HAS_MPL = False


# ===========================================================================
# Item 8 – format_species_label
# ===========================================================================


def format_species_label(formula: str) -> str:
    """Convert a Reaktoro species formula to a matplotlib mathtext string.

    Examples
    --------
    >>> format_species_label("Fe+2")
    'Fe$^{2+}$'
    >>> format_species_label("HCO3-")
    'HCO$_3^-$'
    >>> format_species_label("Fe")
    'Fe'
    >>> format_species_label("CO2")
    'CO$_2$'
    >>> format_species_label("H2O")
    'H$_2$O'
    >>> format_species_label("SO4-2")
    'SO$_4^{2-}$'
    """
    # -----------------------------------------------------------------------
    # Step 1: separate charge suffix from the rest.
    #   Reaktoro convention: "+" / "-" alone, or sign-then-digits "+N"/"-N".
    #   We do NOT treat digit-then-sign ("3-") as charge because in standard
    #   geochemical notation the digit is usually a stoichiometric subscript
    #   (e.g. "HCO3-" → body "HCO3", charge "-").
    # -----------------------------------------------------------------------
    charge_suffix = ""
    body = formula

    # Match: optional digit(s) at end acting as charge magnitude when preceded
    # by sign, e.g. "+2", "-2", "+" or "-" alone.
    m = re.search(r"([+-]\d+|[+-])$", formula)
    if m:
        raw_charge = m.group(0)
        body = formula[: m.start()]
        if re.match(r"^[+-]$", raw_charge):
            sign = raw_charge
            n = ""
        else:  # sign-first, e.g. "+2" or "-2"
            sign = raw_charge[0]
            n = raw_charge[1:]
        charge_suffix = f"^{{{n}{sign}}}" if n else f"^{sign}"

    # -----------------------------------------------------------------------
    # Step 2: subscript digits that follow letters inside the body.
    #   "HCO3" → "HCO$_3$"  /  "H2O" → "H$_2$O"
    # -----------------------------------------------------------------------
    # Replace runs of digits preceded by a letter with LaTeX subscripts.
    # No braces for single digits to keep the output compact.
    body_tex = re.sub(r"(?<=[A-Za-z])(\d+)", r"$_\1$", body)

    # -----------------------------------------------------------------------
    # Step 3: Combine body and charge.
    # -----------------------------------------------------------------------
    if not charge_suffix:
        return body_tex

    # If body_tex ends with "$" (already in math mode) we can merge into
    # the last math segment; otherwise wrap charge in its own math group.
    if body_tex.endswith("$"):
        # E.g. "HCO$_3$" + "^{2-}" → "HCO$_3^{2-}$"
        return body_tex[:-1] + charge_suffix + "$"
    else:
        # E.g. "Fe" + "^{2+}" → "Fe$^{2+}$"
        return body_tex + "$" + charge_suffix + "$"


# ===========================================================================
# Item 9 – axis_label
# ===========================================================================

_AXIS_LABELS = {
    "pH": "pH",
    "Eh": r"$E_\mathrm{h}$ (V)",
    "pe": r"$p\varepsilon$",
    "T": r"$T$ (°C)",
    "TK": r"$T$ (K)",
    "P": r"$P$ (bar)",
    "PPa": r"$P$ (Pa)",
    "logfO2": r"$\log f_{\mathrm{O_2}}$",
    "logfH2": r"$\log f_{\mathrm{H_2}}$",
    "logfCO2": r"$\log f_{\mathrm{CO_2}}$",
    "logaH2O": r"$\log a_{\mathrm{H_2O}}$",
    "IS": r"Ionic strength (mol/kg)",
}


def axis_label(variable: str) -> str:
    """Return a formatted axis-label string for *variable*.

    Parameters
    ----------
    variable:
        One of ``"pH"``, ``"Eh"``, ``"pe"``, ``"T"``, ``"TK"``, ``"P"``,
        ``"PPa"``, ``"logfO2"``, ``"logfH2"``, ``"logfCO2"``, ``"logaH2O"``,
        ``"IS"``.  Falls back to *variable* unchanged if not recognised.
    """
    return _AXIS_LABELS.get(variable, variable)


def _anchored_label_position(px, py, x_min, x_max, y_min, y_max):
    """Clamp text position inside axes and select edge-aware text alignment."""
    xr = x_max - x_min
    yr = y_max - y_min
    mx = 0.03 * xr if xr > 0 else 0.0
    my = 0.03 * yr if yr > 0 else 0.0

    px = min(max(px, x_min + mx), x_max - mx)
    py = min(max(py, y_min + my), y_max - my)

    if px <= x_min + 1.5 * mx:
        ha = "left"
    elif px >= x_max - 1.5 * mx:
        ha = "right"
    else:
        ha = "center"

    if py <= y_min + 1.5 * my:
        va = "bottom"
    elif py >= y_max - 1.5 * my:
        va = "top"
    else:
        va = "center"

    return px, py, ha, va


# ===========================================================================
# Item 6 – water_lines (analytic Eh boundaries)
# ===========================================================================

# Faraday constant, universal gas constant
_F = 96485.0  # C/mol
_R = 8.31446  # J/(mol·K)


def _nernst_slope(T_K: float) -> float:
    """Return (RT/F) * ln(10) = 2.303 RT/F in volts."""
    return 2.302585 * _R * T_K / _F


def water_lines(
    ax: "Axes",
    T_K: float = 298.15,
    P_Pa: float = 1e5,
    pH_range: Optional[Sequence[float]] = None,
    **line_kw,
) -> None:
    """Draw the O₂ and H₂ water-stability lines on *ax*.

    The standard potentials include a pressure correction for activities of
    dissolved O₂ and H₂ at the given pressure *P_Pa*.

    Parameters
    ----------
    ax:
        Matplotlib axes object.
    T_K:
        Temperature in kelvin (default 298.15 K = 25 °C).
    P_Pa:
        Pressure in Pascal (default 1e5 Pa = 1 bar).
    pH_range:
        ``[pH_min, pH_max]``; if *None* uses the current x-axis limits.
    **line_kw:
        Extra keyword arguments forwarded to ``ax.plot()``.
    """
    if not _HAS_MPL:
        raise ImportError("matplotlib is required for water_lines()")

    slope = _nernst_slope(T_K)  # 0.05916 V at 25 °C

    # E°(O₂/H₂O) = 1.2291 V at 25 °C (Pourbaix convention).
    # Pressure correction: ΔE = (RT/4F) ln(P/P°).
    E0_O2 = 1.2291 + (_R * T_K / (4 * _F)) * np.log(P_Pa / 1e5)
    # Upper line:  O₂ + 4H⁺ + 4e⁻ → 2H₂O   Eh = E°_O2 - slope·pH
    # Lower line:  2H⁺ + 2e⁻ → H₂             Eh = -slope·pH  (E°=0 by definition)
    E0_H2 = 0.0 + (_R * T_K / (2 * _F)) * np.log(P_Pa / 1e5)

    if pH_range is None:
        pH_range = ax.get_xlim()

    pH = np.array([pH_range[0], pH_range[1]], dtype=float)

    kw = dict(color="blue", linestyle="--", linewidth=1.0, alpha=0.6)
    kw.update(line_kw)

    ax.plot(pH, E0_O2 - slope * pH, **kw, label=r"O$_2$/H$_2$O")
    ax.plot(pH, E0_H2 - slope * pH, **kw, label=r"H$_2$O/H$_2$")


# ===========================================================================
# Internal colour-palette helper
# ===========================================================================


def _make_colormap(n: int, palette: str = "tab20") -> list:
    """Return *n* RGBA colours from a matplotlib colormap."""
    cmap = _mpl.colormaps.get_cmap(palette)
    return [cmap(i / max(n - 1, 1)) for i in range(n)]


# ===========================================================================
# Item 5 / 6 / 7 – PredominancePlot
# ===========================================================================


class PredominancePlot:
    """CHNOSZ-style Eh–pH (or any 2-D grid) predominance / stability diagram.

    Typical workflow::

        grid_result = solver.sweepPHEhGrid(state, pH_vals, Eh_vals, "V")
        species = ["Fe+2", "Fe+3", "FeO", "Fe2O3", "Fe3O4"]
        predominance = grid_result.predominantSpeciesGrid(species)

        pp = PredominancePlot(grid_result.xvalues, grid_result.yvalues,
                              predominance, species)
        fig, ax = pp.plot()
        pp.add_water_lines(ax, T_K=298.15, P_Pa=1e5)
        pp.add_contours(ax, grid_result.saturationIndexGrid("Hematite"),
                        levels=[0.0], colors="black")
        plt.show()

    Parameters
    ----------
    xvalues:
        1-D array of x-axis values (pH).
    yvalues:
        1-D array of y-axis values (Eh, V).
    predominance:
        2-D integer array of shape ``(nx, ny)`` – the index into *species*
        for the predominant species at each grid point.  ``-1`` or ``nan``
        indicates a failed solve.
    species:
        Sequence of species-name strings whose indices correspond to the
        integer values in *predominance*.
    xlabel, ylabel:
        Axis labels; defaults to ``axis_label("pH")`` and
        ``axis_label("Eh")``.
    palette:
        Name of a matplotlib colormap to sample species colours from.
    """

    def __init__(
        self,
        xvalues,
        yvalues,
        predominance,
        species: Sequence[str],
        xlabel: str = "pH",
        ylabel: str = "Eh",
        palette: str = "tab20",
    ):
        if not _HAS_MPL:
            raise ImportError("matplotlib is required for PredominancePlot")
        self.xvalues = np.asarray(xvalues, dtype=float)
        self.yvalues = np.asarray(yvalues, dtype=float)
        self.predominance = np.asarray(predominance, dtype=float)
        self.species = list(species)
        self.xlabel = axis_label(xlabel)
        self.ylabel = axis_label(ylabel)
        self.palette = palette

    # ------------------------------------------------------------------
    # Item 5 – plot()
    # ------------------------------------------------------------------

    def plot(
        self,
        ax: Optional["Axes"] = None,
        figsize: tuple = (7, 5),
        label_min_fraction: float = 0.01,
        label_fontsize: int = 10,
        boundary_color: str = "black",
        boundary_linewidth: float = 0.8,
    ) -> tuple:
        """Render the predominance diagram.

        Parameters
        ----------
        ax:
            Existing axes to draw on; if *None* a new figure is created.
        figsize:
            Figure size passed to ``plt.subplots()`` when *ax* is *None*.
        label_min_fraction:
            Minimum fractional area (0–1) a stability field must occupy
            before its label is drawn.  Prevents tiny-field clutter.
        label_fontsize:
            Font size for species labels.
        boundary_color:
            Line colour for stability-field boundaries.
        boundary_linewidth:
            Line width for stability-field boundaries.

        Returns
        -------
        fig, ax:
            The :class:`~matplotlib.figure.Figure` and
            :class:`~matplotlib.axes.Axes` objects.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        n_species = len(self.species)
        colours = _make_colormap(max(n_species, 2), self.palette)

        nx = len(self.xvalues)
        ny = len(self.yvalues)

        # predominance has shape (nx, ny); imshow expects (rows=y, cols=x).
        # Transpose so rows = Eh, cols = pH; origin='lower' puts low-Eh at bottom.
        Z = self.predominance.T  # shape (ny, nx), float (may contain nan/-1)

        # Build a masked integer array – mask failed points.
        Zint = np.ma.masked_invalid(Z)
        Zint = np.ma.masked_less(Zint, 0)

        # Build discrete colour map of exactly n_species colours.
        cmap = mcolors.ListedColormap(colours[:n_species])
        norm = mcolors.BoundaryNorm(
            boundaries=np.arange(-0.5, n_species + 0.5, 1),
            ncolors=n_species,
        )

        x_min, x_max = self.xvalues[0], self.xvalues[-1]
        y_min, y_max = self.yvalues[0], self.yvalues[-1]

        im = ax.imshow(
            Zint,
            origin="lower",
            extent=[x_min, x_max, y_min, y_max],
            aspect="auto",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
        )

        # ---- Boundary lines: exactly one line per species transition ------
        # Draw a separate binary contour for each pair of neighbouring species
        # so there is no risk of spurious iso-crossings from blurring a
        # non-monotone integer field.
        X, Y = np.meshgrid(self.xvalues, self.yvalues)  # (ny, nx)
        Zfloat = np.ma.filled(Zint.astype(float), np.nan)
        if not np.all(np.isnan(Zfloat)):
            try:
                from scipy.ndimage import gaussian_filter as _gf

                _has_scipy = True
            except ImportError:
                _has_scipy = False
            seen_pairs = set()
            # Collect all horizontally or vertically adjacent differing pairs.
            for drow, dcol in [(0, 1), (1, 0)]:
                a = Zfloat[: -drow or None, : -dcol or None]
                b = Zfloat[drow:, dcol:]
                mask = np.isfinite(a) & np.isfinite(b) & (a != b)
                pairs = set(zip(a[mask].astype(int), b[mask].astype(int)))
                seen_pairs |= {tuple(sorted(p)) for p in pairs}
            for lo, hi in seen_pairs:
                binary = np.where(
                    Zfloat == hi, 1.0, np.where(Zfloat == lo, 0.0, np.nan)
                )
                if np.all(np.isnan(binary)):
                    continue
                if _has_scipy:
                    binary = _gf(
                        np.nan_to_num(binary, nan=0.5), sigma=0.6, mode="nearest"
                    )
                try:
                    ax.contour(
                        X,
                        Y,
                        binary,
                        levels=[0.5],
                        colors=boundary_color,
                        linewidths=boundary_linewidth,
                    )
                except Exception:
                    pass

        # ---- Species labels at centroid of each field ---------------------
        total_valid = np.sum(~Zint.mask) if np.ma.is_masked(Zint) else Zint.size
        for idx, name in enumerate(self.species):
            mask = Zint == idx
            count = np.sum(mask)
            if total_valid == 0 or count == 0:
                continue
            fraction = count / total_valid
            if fraction < label_min_fraction:
                continue
            rows, cols = np.where(mask)
            cy = np.mean(rows)  # row → y axis
            cx = np.mean(cols)  # col → x axis
            # Map pixel indices back to data coordinates.
            px = x_min + cx / (nx - 1) * (x_max - x_min) if nx > 1 else x_min
            py = y_min + cy / (ny - 1) * (y_max - y_min) if ny > 1 else y_min
            px, py, ha, va = _anchored_label_position(
                px, py, x_min, x_max, y_min, y_max
            )
            ax.text(
                px,
                py,
                format_species_label(name),
                ha=ha,
                va=va,
                fontsize=label_fontsize,
                fontweight="bold",
                color="black",
                clip_on=False,
            )

        # ---- Legend patches -----------------------------------------------
        patches = []
        for idx, name in enumerate(self.species):
            mask = Zint == idx
            if np.sum(mask) == 0:
                continue
            patches.append(
                mpatches.Patch(
                    facecolor=colours[idx],
                    label=format_species_label(name),
                )
            )
        if patches:
            ax.legend(
                handles=patches,
                loc="upper left",
                fontsize=8,
                framealpha=0.7,
            )

        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.tick_params(which="both", direction="in", top=True, right=True)
        try:
            ax.minorticks_on()
        except Exception:
            pass

        return fig, ax

    # ------------------------------------------------------------------
    # Item 6 – add_water_lines()
    # ------------------------------------------------------------------

    def add_water_lines(
        self,
        ax: "Axes",
        T_K: float = 298.15,
        P_Pa: float = 1e5,
        **line_kw,
    ) -> None:
        """Overlay O₂ / H₂ water-stability lines on *ax*.

        Delegates to :func:`water_lines`; pH range taken from current axes
        x-limits.
        """
        water_lines(ax, T_K=T_K, P_Pa=P_Pa, **line_kw)

    # ------------------------------------------------------------------
    # Item 7 – add_contours()
    # ------------------------------------------------------------------

    def add_contours(
        self,
        ax: "Axes",
        grid: "np.ndarray",
        levels: Union[int, Sequence[float]] = 5,
        fmt: Optional[str] = None,
        label_contours: bool = True,
        **contour_kw,
    ) -> None:
        """Overlay contour lines of a scalar grid on *ax*.

        Parameters
        ----------
        ax:
            Axes to draw on (should be the same as returned by :meth:`plot`).
        grid:
            2-D array of shape ``(nx, ny)`` – e.g. from
            ``grid_result.saturationIndexGrid(mineral)`` or
            ``grid_result.logActivityGrid(species)``.
        levels:
            Contour levels; int or explicit list.
        fmt:
            Format string for contour labels (e.g. ``"%.1f"``).
        label_contours:
            Whether to draw inline labels on the contour lines.
        **contour_kw:
            Extra keyword arguments forwarded to ``ax.contour()``.
        """
        g = np.asarray(grid, dtype=float)
        # grid has shape (nx, ny); meshgrid expects (ny, nx) for contour.
        X, Y = np.meshgrid(self.xvalues, self.yvalues)
        cs = ax.contour(X, Y, g.T, levels=levels, **contour_kw)
        if label_contours:
            ax.clabel(cs, inline=True, fmt=fmt or "%.2g", fontsize=8)


# ===========================================================================
# Item 10 – SpeciationPlot
# ===========================================================================


class SpeciationPlot:
    """1-D speciation plot of log activity or fractional abundance (α).

    Parameters
    ----------
    xvalues:
        1-D array of x-axis sweep values.
    xlabel:
        Variable name (used by :func:`axis_label`).
    palette:
        Matplotlib colormap name for line colours.
    """

    def __init__(
        self,
        xvalues,
        xlabel: str = "pH",
        palette: str = "tab10",
    ):
        if not _HAS_MPL:
            raise ImportError("matplotlib is required for SpeciationPlot")
        self.xvalues = np.asarray(xvalues, dtype=float)
        self.xlabel = axis_label(xlabel)
        self.palette = palette
        self._series: list = []  # list of dicts

    def add_log_activity(self, species: str, log_activity: "np.ndarray") -> None:
        """Add a log₁₀ activity series for *species*.

        Parameters
        ----------
        species:
            Species name (used for legend label).
        log_activity:
            1-D array with the same length as *xvalues*.
        """
        self._series.append(
            {
                "name": species,
                "y": np.asarray(log_activity, dtype=float),
                "mode": "log_activity",
            }
        )

    def add_fraction(self, species: str, fraction: "np.ndarray") -> None:
        """Add a fractional-abundance (α) series for *species*.

        Parameters
        ----------
        species:
            Species name.
        fraction:
            1-D array of fractional abundance values (0–1), same length as
            *xvalues*.
        """
        self._series.append(
            {
                "name": species,
                "y": np.asarray(fraction, dtype=float),
                "mode": "fraction",
            }
        )

    def add_element_molality(
        self,
        element: str,
        molalities: "np.ndarray",
        log10: bool = True,
    ) -> None:
        """Add a total dissolved element-molality series.

        Parameters
        ----------
        element:
            Element symbol used for the legend label.
        molalities:
            1-D array of molality values (mol/kg water), same length as
            *xvalues*.  The values may already be log10-transformed if
            *log10* is *False*.
        log10:
            When *True* (default) the values are log₁₀-transformed before
            storage so the plot y-axis is consistent with log-activity mode.
        """
        y = np.asarray(molalities, dtype=float)
        if log10:
            with np.errstate(divide="ignore"):
                y = np.log10(y)
        self._series.append(
            {
                "name": element,
                "y": y,
                "mode": "log_activity",  # same axes as log activity
            }
        )

    def plot(
        self,
        ax: Optional["Axes"] = None,
        figsize: tuple = (8, 5),
        mode: str = "log_activity",
        ylim: Optional[tuple] = None,
        min_log_activity: float = -10.0,
        linewidth: float = 1.5,
    ) -> tuple:
        """Render the speciation diagram.

        Parameters
        ----------
        ax:
            Existing axes; a new figure is created if *None*.
        figsize:
            Figure size when creating a new figure.
        mode:
            ``"log_activity"`` or ``"fraction"``.  Determines which series
            are drawn and the y-axis label.
        ylim:
            Explicit y-axis limits; auto-scaled when *None*.
        min_log_activity:
            Series with maximum log activity below this value are hidden
            (only applies in ``"log_activity"`` mode).
        linewidth:
            Line width for all series.

        Returns
        -------
        fig, ax:
            Figure and axes objects.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        series = [s for s in self._series if s["mode"] == mode]
        if not series:
            warnings.warn(
                f"SpeciationPlot: no series with mode='{mode}' to plot.",
                stacklevel=2,
            )
            return fig, ax

        colours = _make_colormap(max(len(series), 2), self.palette)

        for i, s in enumerate(series):
            y = s["y"]
            name = s["name"]
            if mode == "log_activity" and np.nanmax(y) < min_log_activity:
                continue
            ax.plot(
                self.xvalues,
                y,
                color=colours[i],
                linewidth=linewidth,
                label=format_species_label(name),
            )

        ax.set_xlabel(self.xlabel)

        if mode == "log_activity":
            ax.set_ylabel(r"$\log_{10} a$")
        else:
            ax.set_ylabel(r"Mole fraction $\alpha$")

        if ylim is not None:
            ax.set_ylim(*ylim)

        ax.tick_params(which="both", direction="in", top=True, right=True)
        try:
            ax.minorticks_on()
        except Exception:
            pass

        ax.legend(fontsize=8, framealpha=0.7)
        return fig, ax

    # ------------------------------------------------------------------
    # Convenience: build directly from an EquilibriumSweepResult
    # ------------------------------------------------------------------

    @classmethod
    def from_sweep_result(
        cls,
        sweep_result,
        species: Sequence[str],
        xlabel: str = "pH",
        xvalues=None,
        palette: str = "tab10",
        elements: Sequence[str] = (),
    ) -> "SpeciationPlot":
        """Construct a :class:`SpeciationPlot` from an
        :class:`EquilibriumSweepResult`.

        Parameters
        ----------
        sweep_result:
            A 1-D ``EquilibriumSweepResult`` (returned by
            ``EquilibriumSweepSolver.sweep()``).
        species:
            List of species names to include.
        xlabel:
            Variable name for the x-axis.
        xvalues:
            Explicit x-axis values; if *None* the index is used.
        palette:
            Colormap name.
        """
        if xvalues is None:
            n = len(sweep_result.states)
            xvalues = np.arange(n, dtype=float)
        else:
            xvalues = np.asarray(xvalues, dtype=float)

        sp = cls(xvalues, xlabel=xlabel, palette=palette)

        for name in species:
            log_a = np.array(
                [
                    float(
                        __import__("reaktoro")
                        .ChemicalProps(state)
                        .speciesActivityLg(name)
                    )
                    for state in sweep_result.states
                ]
            )
            sp.add_log_activity(name, log_a)

        for elem in elements:
            mol = sweep_result.elementMolalityArray(elem)
            sp.add_element_molality(elem, mol, log10=True)

        return sp


# ===========================================================================
# T3-1 – SolubilityPlot
# ===========================================================================


class SolubilityPlot:
    """Filled-contour diagram of total dissolved element concentration.

    Visualises log₁₀ of the total dissolved molality of an element (sum over
    all aqueous species containing it) across a 2-D sweep grid.  Iso-solubility
    contour lines are optionally overlaid.

    Typical workflow::

        grid = solver.sweepPHEhGrid(state, pH_vals, Eh_vals, "V")
        sp = SolubilityPlot.from_grid_result(grid, element="Fe")
        fig, ax = sp.plot(levels=20)
        plt.show()

    Parameters
    ----------
    xvalues:
        1-D array of x-axis values (pH).
    yvalues:
        1-D array of y-axis values (Eh, V).
    molality_grid:
        2-D array of shape ``(nx, ny)`` with total element molality (mol/kg).
    element:
        Element symbol used for the colorbar label.
    xlabel, ylabel:
        Axis variable names; passed through :func:`axis_label`.
    """

    def __init__(
        self,
        xvalues,
        yvalues,
        molality_grid,
        element: str = "",
        xlabel: str = "pH",
        ylabel: str = "Eh",
    ):
        if not _HAS_MPL:
            raise ImportError("matplotlib is required for SolubilityPlot")
        self.xvalues = np.asarray(xvalues, dtype=float)
        self.yvalues = np.asarray(yvalues, dtype=float)
        self.molality_grid = np.asarray(molality_grid, dtype=float)
        self.element = element
        self.xlabel = axis_label(xlabel)
        self.ylabel = axis_label(ylabel)

    @classmethod
    def from_grid_result(
        cls,
        grid_result,
        element: str,
        xlabel: str = "pH",
        ylabel: str = "Eh",
    ) -> "SolubilityPlot":
        """Build from an :class:`EquilibriumSweepGridResult`.

        Parameters
        ----------
        grid_result:
            Result of ``solver.sweepPHEhGrid(...)`` or ``sweepLgActivityGrid(...)``
        element:
            Element symbol (e.g. ``"Fe"``, ``"Ca"``, ``"S"``).
        """
        molality = np.asarray(grid_result.elementMolalityGrid(element), dtype=float)
        return cls(
            xvalues=np.asarray(grid_result.xvalues, dtype=float),
            yvalues=np.asarray(grid_result.yvalues, dtype=float),
            molality_grid=molality,
            element=element,
            xlabel=xlabel,
            ylabel=ylabel,
        )

    def plot(
        self,
        ax: Optional["Axes"] = None,
        figsize: tuple = (7, 5),
        levels: int = 20,
        iso_levels: Optional[Sequence[float]] = None,
        cmap: str = "viridis",
        clabel_fmt: str = "%.1f",
        clabel_inline: bool = True,
    ) -> tuple:
        """Render the solubility diagram.

        Parameters
        ----------
        ax:
            Existing axes; a new figure is created if *None*.
        figsize:
            Figure size when creating a new figure.
        levels:
            Number of filled contour levels for the log₁₀ molality field.
        iso_levels:
            Explicit log₁₀ molality values at which to draw labelled iso-solubility
            lines.  If *None*, no extra contour lines are drawn.
        cmap:
            Colormap name for the filled contour plot.
        clabel_fmt:
            Format string for iso-level labels.
        clabel_inline:
            Draw contour labels inline.

        Returns
        -------
        fig, ax
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        # Log-scale the molality, guard against zeros/negatives.
        with np.errstate(divide="ignore", invalid="ignore"):
            log_mol = np.where(
                self.molality_grid > 0,
                np.log10(self.molality_grid),
                np.nan,
            )

        # Shape is (nx, ny); contourf expects (ny, nx) = (rows, cols).
        X, Y = np.meshgrid(self.xvalues, self.yvalues)
        cf = ax.contourf(X, Y, log_mol.T, levels=levels, cmap=cmap)

        cb = fig.colorbar(cf, ax=ax)
        elem_label = self.element if self.element else "element"
        cb.set_label(
            rf"$\log_{{10}}$ molality of {elem_label} (mol/kg)",
            fontsize=9,
        )

        if iso_levels is not None and len(iso_levels) > 0:
            cs = ax.contour(
                X,
                Y,
                log_mol.T,
                levels=sorted(iso_levels),
                colors="white",
                linewidths=1.0,
                linestyles="--",
            )
            ax.clabel(cs, inline=clabel_inline, fmt=clabel_fmt, fontsize=8)

        x_min, x_max = self.xvalues[0], self.xvalues[-1]
        y_min, y_max = self.yvalues[0], self.yvalues[-1]
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.tick_params(which="both", direction="in", top=True, right=True)
        try:
            ax.minorticks_on()
        except Exception:
            pass

        return fig, ax


# ===========================================================================
# T3-2 – ActivityDiagram
# ===========================================================================


class ActivityDiagram(PredominancePlot):
    """Stability-field diagram with log(aX) vs log(aY) axes.

    Wraps :class:`PredominancePlot` with axis labels formatted as
    log activity of the two sweep species.  Data is obtained from
    ``solver.sweepLgActivityGrid(...).predominantSpeciesGrid(...)``.

    Typical workflow::

        # specs must have lgActivity(speciesX) and lgActivity(speciesY) declared
        grid = solver.sweepLgActivityGrid(state,
            "Ca+2", lga_ca_vals, "CO3-2", lga_co3_vals)
        species = ["Calcite", "Dolomite", "Ca+2"]
        predominance = grid.predominantSpeciesGrid(species)

        ad = ActivityDiagram(grid, predominance, species,
                             speciesX="Ca+2", speciesY="CO3-2")
        fig, ax = ad.plot()
        plt.show()

    Parameters
    ----------
    grid_result:
        :class:`EquilibriumSweepGridResult` from ``sweepLgActivityGrid``.
    predominance:
        2-D integer array (nx, ny) from ``grid_result.predominantSpeciesGrid``.
    species:
        List of species names.
    speciesX, speciesY:
        Names of the x and y sweep species (used for axis labels).
    palette:
        Colormap name.
    """

    def __init__(
        self,
        grid_result,
        predominance,
        species: Sequence[str],
        speciesX: str = "",
        speciesY: str = "",
        palette: str = "tab20",
    ):
        xlabel = (
            rf"$\log_{{10}}\,a_{{\mathrm{{{format_species_label(speciesX)}}}}}$"
            if speciesX
            else r"$\log_{10}\,a_X$"
        )
        ylabel = (
            rf"$\log_{{10}}\,a_{{\mathrm{{{format_species_label(speciesY)}}}}}$"
            if speciesY
            else r"$\log_{10}\,a_Y$"
        )
        # Bypass axis_label lookup — labels are already formatted.
        super().__init__(
            xvalues=np.asarray(grid_result.xvalues, dtype=float),
            yvalues=np.asarray(grid_result.yvalues, dtype=float),
            predominance=predominance,
            species=species,
            xlabel="pH",  # placeholder; overwritten below
            ylabel="Eh",
            palette=palette,
        )
        # Overwrite with the properly formatted log-activity labels.
        self.xlabel = xlabel
        self.ylabel = ylabel

    @classmethod
    def from_grid_result(
        cls,
        grid_result,
        species: Sequence[str],
        speciesX: str = "",
        speciesY: str = "",
        palette: str = "tab20",
    ) -> "ActivityDiagram":
        """Convenience constructor: compute predominance from *grid_result*."""
        predominance = np.asarray(
            grid_result.predominantSpeciesGrid(list(species)), dtype=float
        )
        return cls(grid_result, predominance, species, speciesX, speciesY, palette)


# ===========================================================================
# T3-3 – MosaicPlot
# ===========================================================================


class MosaicPlot:
    """Layered predominance diagram for multiple species groups.

    In CHNOSZ, ``mosaic()`` draws separate stability fields for overlapping
    groups (e.g., Fe-aqueous ions on top of Fe-oxide minerals) using different
    colours and alpha blending.  This class replicates that pattern.

    Each *layer* is a dict with keys:

    * ``"species"`` – list of species names for this group
    * ``"predominance"`` – 2-D array (nx, ny) from ``predominantSpeciesGrid``
    * ``"palette"`` – (optional) colormap name, defaults to ``"tab20"``
    * ``"alpha"`` – (optional) float 0–1, defaults to 1.0 for the first layer
      and 0.7 for subsequent layers

    Typical workflow::

        aq_species  = ["Fe+2", "Fe+3", "FeOH+"]
        min_species = ["Hematite", "Magnetite", "Siderite"]

        layers = [
            {"species": min_species,
             "predominance": grid.predominantSpeciesGrid(min_species),
             "palette": "Pastel1", "alpha": 1.0},
            {"species": aq_species,
             "predominance": grid.predominantSpeciesGrid(aq_species),
             "palette": "tab10", "alpha": 0.65},
        ]
        mp = MosaicPlot(grid.xvalues, grid.yvalues, layers)
        fig, ax = mp.plot()
        mp.add_water_lines(ax)
        plt.show()

    Parameters
    ----------
    xvalues, yvalues:
        1-D axis arrays.
    layers:
        List of layer dicts (see above), rendered in order (first = bottom).
    xlabel, ylabel:
        Axis variable names passed through :func:`axis_label`.
    """

    def __init__(
        self,
        xvalues,
        yvalues,
        layers: list,
        xlabel: str = "pH",
        ylabel: str = "Eh",
    ):
        if not _HAS_MPL:
            raise ImportError("matplotlib is required for MosaicPlot")
        self.xvalues = np.asarray(xvalues, dtype=float)
        self.yvalues = np.asarray(yvalues, dtype=float)
        self.layers = layers
        self.xlabel = axis_label(xlabel)
        self.ylabel = axis_label(ylabel)

    def plot(
        self,
        ax: Optional["Axes"] = None,
        figsize: tuple = (7, 5),
        label_min_fraction: float = 0.01,
        label_fontsize: int = 9,
        boundary_color: str = "black",
        boundary_linewidth: float = 0.8,
        show_legend: bool = True,
    ) -> tuple:
        """Render all layers onto a single axes.

        Parameters
        ----------
        ax:
            Existing axes; a new figure is created if *None*.
        figsize:
            Figure size when creating a new figure.
        label_min_fraction:
            Minimum fractional area for field labels.
        label_fontsize:
            Font size for labels.
        boundary_color:
            Colour for stability-field boundaries within each layer.
        boundary_linewidth:
            Width for boundary contour lines.
        show_legend:
            Whether to draw a combined legend for all layers.

        Returns
        -------
        fig, ax
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        nx = len(self.xvalues)
        ny = len(self.yvalues)
        x_min, x_max = self.xvalues[0], self.xvalues[-1]
        y_min, y_max = self.yvalues[0], self.yvalues[-1]
        X, Y = np.meshgrid(self.xvalues, self.yvalues)

        # all_patches stores (group_label, patch) tuples so the final
        # legend can be split into titled groups per layer.
        all_patches: list = []  # list of (group_label, mpatches.Patch)

        for layer_idx, layer in enumerate(self.layers):
            species = list(layer["species"])
            predominance = np.asarray(layer["predominance"], dtype=float)
            palette = layer.get("palette", "tab20")
            default_alpha = 1.0 if layer_idx == 0 else 0.65
            alpha = float(layer.get("alpha", default_alpha))
            layer_boundary_color = layer.get("boundary_color", boundary_color)
            layer_boundary_linewidth = float(
                layer.get("boundary_linewidth", boundary_linewidth)
            )
            layer_boundary_linestyle = layer.get("boundary_linestyle", "-")
            n_species = len(species)
            colours = _make_colormap(max(n_species, 2), palette)

            Z = predominance.T  # (ny, nx)
            Zint = np.ma.masked_invalid(Z)
            Zint = np.ma.masked_less(Zint, 0)

            cmap = mcolors.ListedColormap(colours[:n_species])
            norm = mcolors.BoundaryNorm(
                boundaries=np.arange(-0.5, n_species + 0.5, 1),
                ncolors=n_species,
            )

            ax.imshow(
                Zint,
                origin="lower",
                extent=[x_min, x_max, y_min, y_max],
                aspect="auto",
                cmap=cmap,
                norm=norm,
                interpolation="nearest",
                alpha=alpha,
            )

            # Boundary lines for this layer — one line per pair of
            # adjacent species, drawn on a binary per-pair mask.
            Zfloat = np.ma.filled(Zint.astype(float), np.nan)
            if not np.all(np.isnan(Zfloat)):
                try:
                    from scipy.ndimage import gaussian_filter as _gf

                    _has_scipy = True
                except ImportError:
                    _has_scipy = False
                seen_pairs: set = set()
                for drow, dcol in [(0, 1), (1, 0)]:
                    a = Zfloat[: -drow or None, : -dcol or None]
                    b = Zfloat[drow:, dcol:]
                    pair_mask = np.isfinite(a) & np.isfinite(b) & (a != b)
                    seen_pairs |= {
                        tuple(sorted(p))
                        for p in zip(a[pair_mask].astype(int), b[pair_mask].astype(int))
                    }
                for lo, hi in seen_pairs:
                    binary = np.where(
                        Zfloat == hi, 1.0, np.where(Zfloat == lo, 0.0, np.nan)
                    )
                    if np.all(np.isnan(binary)):
                        continue
                    if _has_scipy:
                        binary = _gf(
                            np.nan_to_num(binary, nan=0.5), sigma=0.6, mode="nearest"
                        )
                    try:
                        ax.contour(
                            X,
                            Y,
                            binary,
                            levels=[0.5],
                            colors=layer_boundary_color,
                            linewidths=layer_boundary_linewidth,
                            linestyles=layer_boundary_linestyle,
                        )
                    except Exception:
                        pass

            # Optional: Draw thin edges around each region in this layer
            # to visually frame the fields and improve readability under
            # transparency. Useful for the base (mineral) layer.
            if layer.get("draw_edges", False):
                try:
                    # Find all region boundaries in the integer predominance field
                    # and draw thin dark outlines to visually frame them.
                    from scipy.ndimage import binary_erosion as _be

                    Zfloat_filled = np.nan_to_num(Zfloat, nan=-1.0)
                    valid_mask = Zfloat_filled >= 0
                    for idx in range(n_species):
                        region = Zfloat_filled == idx
                        if not np.any(region):
                            continue
                        # Erode slightly so only the boundary remains.
                        eroded = _be(region, iterations=1)
                        boundary = region & ~eroded
                        if np.any(boundary):
                            edge_contour = np.where(boundary, 1.0, np.nan)
                            ax.contour(
                                X,
                                Y,
                                edge_contour,
                                levels=[0.5],
                                colors="black",
                                linewidths=0.4,
                                alpha=0.3,
                            )
                except Exception:
                    pass

            # Species labels.
            total_valid = np.sum(~Zint.mask) if np.ma.is_masked(Zint) else Zint.size
            for idx, name in enumerate(species):
                mask = Zint == idx
                count = np.sum(mask)
                if total_valid == 0 or count == 0:
                    continue
                if count / total_valid < label_min_fraction:
                    continue
                rows, cols = np.where(mask)
                cy, cx = np.mean(rows), np.mean(cols)
                px = x_min + cx / max(nx - 1, 1) * (x_max - x_min)
                py = y_min + cy / max(ny - 1, 1) * (y_max - y_min)
                px, py, ha, va = _anchored_label_position(
                    px, py, x_min, x_max, y_min, y_max
                )
                ax.text(
                    px,
                    py,
                    format_species_label(name),
                    ha=ha,
                    va=va,
                    fontsize=label_fontsize,
                    fontweight="bold",
                    color="black",
                    clip_on=False,
                )

            if show_legend:
                group_label = layer.get("label", None)
                for idx, name in enumerate(species):
                    if np.sum(Zint == idx) == 0:
                        continue
                    all_patches.append(
                        (
                            group_label,
                            mpatches.Patch(
                                facecolor=colours[idx],
                                alpha=alpha,
                                label=format_species_label(name),
                            ),
                        )
                    )

        if show_legend and all_patches:
            # Check whether any layer supplied a group label.
            has_groups = any(g is not None for g, _ in all_patches)
            if not has_groups:
                # Simple flat legend — original behaviour.
                ax.legend(
                    handles=[p for _, p in all_patches],
                    loc="upper left",
                    fontsize=8,
                    framealpha=0.7,
                )
            else:
                # Grouped legend: one titled sub-block per named group.
                from matplotlib.lines import Line2D as _L2D

                legend_handles = []
                seen_groups: list = []
                for group_label, patch in all_patches:
                    if group_label not in seen_groups:
                        seen_groups.append(group_label)
                if None in seen_groups:
                    seen_groups.remove(None)
                    seen_groups.append(None)  # unnamed group last

                for group_label in seen_groups:
                    # Insert a transparent title handle; text is styled bold
                    # via the legend handler map rather than mathtext so that
                    # spaces in the group label are preserved.
                    title_text = group_label if group_label is not None else "Other"
                    title_patch = mpatches.Patch(
                        color="none",
                        label=title_text,
                    )
                    legend_handles.append(title_patch)
                    for g, patch in all_patches:
                        if g == group_label:
                            legend_handles.append(patch)

                leg = ax.legend(
                    handles=legend_handles,
                    loc="upper left",
                    fontsize=8,
                    framealpha=0.7,
                    handlelength=1.4,
                    handleheight=0.85,
                )
                # Bold the group-title entries (transparent patches).
                title_labels = set(seen_groups)
                for text in leg.get_texts():
                    if text.get_text() in title_labels:
                        text.set_fontweight("bold")
                        text.set_fontsize(8.5)

        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.tick_params(which="both", direction="in", top=True, right=True)
        try:
            ax.minorticks_on()
        except Exception:
            pass

        return fig, ax

    def add_water_lines(
        self, ax: "Axes", T_K: float = 298.15, P_Pa: float = 1e5, **line_kw
    ) -> None:
        """Overlay O₂ / H₂ water-stability lines. Delegates to :func:`water_lines`."""
        water_lines(ax, T_K=T_K, P_Pa=P_Pa, **line_kw)


# ===========================================================================
# T3-4 – saturation_curve (standalone helper)
# ===========================================================================


def saturation_curve(
    ax: "Axes",
    xvalues,
    yvalues,
    si_grid,
    label: Optional[str] = None,
    si_value: float = 0.0,
    label_inline: bool = True,
    **line_kw,
) -> None:
    """Draw an iso-saturation contour line (SI = *si_value*) on *ax*.

    This is a convenience wrapper around :func:`matplotlib.axes.Axes.contour`
    that extracts the SI = 0 boundary of a mineral from a saturation-index
    grid produced by :meth:`EquilibriumSweepGridResult.saturationIndexGrid`.

    Parameters
    ----------
    ax:
        Matplotlib axes to draw on.
    xvalues:
        1-D array of x-axis values (pH or log aX).
    yvalues:
        1-D array of y-axis values (Eh or log aY).
    si_grid:
        2-D array of shape ``(nx, ny)`` from ``saturationIndexGrid(mineral)``.
    label:
        Legend / inline label for the contour line.  If *None*, no label is
        drawn.
    si_value:
        The saturation-index value to contour.  Default is 0 (saturation).
    label_inline:
        Draw the label inline on the contour line.
    **line_kw:
        Extra keyword arguments forwarded to ``ax.contour()`` (e.g.
        ``colors``, ``linewidths``, ``linestyles``).
    """
    if not _HAS_MPL:
        raise ImportError("matplotlib is required for saturation_curve()")

    xvalues = np.asarray(xvalues, dtype=float)
    yvalues = np.asarray(yvalues, dtype=float)
    si = np.asarray(si_grid, dtype=float)

    X, Y = np.meshgrid(xvalues, yvalues)  # (ny, nx)
    kw = dict(colors="black", linewidths=1.2, linestyles="-")
    kw.update(line_kw)

    cs = ax.contour(X, Y, si.T, levels=[si_value], **kw)

    if label is not None:
        if label_inline:
            ax.clabel(cs, inline=True, fmt=lambda _v: label, fontsize=8)
        else:
            # Attach label to the first path segment for the legend.
            for coll in cs.collections:
                coll.set_label(label)
            ax.legend(fontsize=8, framealpha=0.7)


# ===========================================================================
# T4-2 -- TPDiagram
# ===========================================================================


class TPDiagram(PredominancePlot):
    """Stability-field diagram with Temperature (x) vs Pressure (y) axes.

    Wraps :class:`PredominancePlot` with axis labels for T and P.
    Data is obtained from ``solver.sweepTPGrid(...).predominantSpeciesGrid(...)``.

    Typical workflow::

        grid = solver.sweepTPGrid(state,
            np.linspace(273.15, 573.15, 60), "K",
            np.linspace(1e5, 1e8, 60), "Pa")
        species = ["Calcite", "Aragonite", "Quartz"]
        td = TPDiagram.from_grid_result(grid, species, T_unit="K", P_unit="Pa")
        fig, ax = td.plot()
        plt.show()

    Parameters
    ----------
    grid_result:
        :class:`EquilibriumSweepGridResult` from ``sweepTPGrid``.
    predominance:
        2-D integer array from ``grid_result.predominantSpeciesGrid(species)``.
    species:
        Ordered list of species names matching the *predominance* grid.
    T_unit:
        Unit for the x-axis label (``"K"`` or ``"C"``).
    P_unit:
        Unit for the y-axis label (``"Pa"`` or ``"bar"``).
    """

    _T_LABELS = {"K": "T (K)", "C": "T (\xb0C)", "celsius": "T (\xb0C)"}
    _P_LABELS = {
        "Pa": "P (Pa)",
        "bar": "P (bar)",
        "MPa": "P (MPa)",
        "GPa": "P (GPa)",
    }

    def __init__(
        self,
        grid_result,
        predominance,
        species: Sequence[str],
        T_unit: str = "K",
        P_unit: str = "Pa",
        **kwargs,
    ):
        xlabel = self._T_LABELS.get(T_unit, "T ({})".format(T_unit))
        ylabel = self._P_LABELS.get(P_unit, "P ({})".format(P_unit))
        super().__init__(
            grid_result.xvalues,
            grid_result.yvalues,
            predominance,
            species,
            xlabel=xlabel,
            ylabel=ylabel,
            **kwargs,
        )
        # Override so axis_label() does not double-format the string.
        self.xlabel = xlabel
        self.ylabel = ylabel

    @classmethod
    def from_grid_result(
        cls,
        grid_result,
        species: Sequence[str],
        T_unit: str = "K",
        P_unit: str = "Pa",
        **kwargs,
    ) -> "TPDiagram":
        """Convenience constructor from a ``sweepTPGrid`` result."""
        predominance = grid_result.predominantSpeciesGrid(species)
        return cls(grid_result, predominance, species, T_unit, P_unit, **kwargs)


# ===========================================================================
# T4-1 -- LogfO2pHDiagram
# ===========================================================================


class LogfO2pHDiagram(PredominancePlot):
    r"""Stability-field diagram with log10(fO2) (x) vs pH (y).

    Wraps :class:`PredominancePlot` with axis labels for oxygen fugacity and
    pH.  Data is obtained from
    ``solver.sweepLogfO2pHGrid(...).predominantSpeciesGrid(...)``.

    EquilibriumSpecs must have ``fugacity("O2")`` and ``pH`` declared before
    constructing the solver.

    Typical workflow::

        grid = solver.sweepLogfO2pHGrid(state,
            np.linspace(-80, 0, 80), "bar",
            np.linspace(0, 14, 70))
        species = ["Fe+2", "Fe+3", "Hematite", "Magnetite", "Pyrite"]
        fd = LogfO2pHDiagram.from_grid_result(grid, species)
        fig, ax = fd.plot()
        plt.show()

    Parameters
    ----------
    grid_result:
        :class:`EquilibriumSweepGridResult` from ``sweepLogfO2pHGrid``.
    predominance:
        2-D integer array from ``grid_result.predominantSpeciesGrid(species)``.
    species:
        Ordered list of species names matching the *predominance* grid.
    """

    def __init__(
        self,
        grid_result,
        predominance,
        species: Sequence[str],
        **kwargs,
    ):
        xlabel = r"$\log_{10} f_{\mathrm{O_2}}$"
        ylabel = "pH"
        super().__init__(
            grid_result.xvalues,
            grid_result.yvalues,
            predominance,
            species,
            xlabel=xlabel,
            ylabel=ylabel,
            **kwargs,
        )
        self.xlabel = xlabel
        self.ylabel = ylabel

    @classmethod
    def from_grid_result(
        cls,
        grid_result,
        species: Sequence[str],
        **kwargs,
    ) -> "LogfO2pHDiagram":
        """Convenience constructor from a ``sweepLogfO2pHGrid`` result."""
        predominance = grid_result.predominantSpeciesGrid(species)
        return cls(grid_result, predominance, species, **kwargs)
