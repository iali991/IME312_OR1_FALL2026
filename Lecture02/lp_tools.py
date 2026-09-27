"""
lp_tools.py -- a small, teaching-oriented toolkit for two-variable Linear Programs.
Parses LP files (CPLEX LP-format subset), solves them by the corner-point method,
and draws the graphical solution.
"""
from __future__ import annotations

import re
import itertools
from dataclasses import dataclass, field

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

TOL = 1e-7

# ----------------------------------------------------------------------------
# 1. Model containers
# ----------------------------------------------------------------------------

@dataclass
class Constraint:
    name: str
    coeffs: dict          # {var: coefficient}
    op: str               # '<=', '>=', '='
    rhs: float

    def as_text(self, order):
        return f"{_expr_text(self.coeffs, order)} {self.op} {_num(self.rhs)}"


@dataclass
class LPModel:
    sense: str = "max"                  # 'max' or 'min'
    obj_name: str = "obj"
    obj: dict = field(default_factory=dict)
    obj_const: float = 0.0
    constraints: list = field(default_factory=list)
    lower: dict = field(default_factory=dict)   # var -> lower bound (-inf allowed)
    upper: dict = field(default_factory=dict)   # var -> upper bound (+inf allowed)
    var_order: list = field(default_factory=list)

    # -- convenience ---------------------------------------------------------
    @property
    def n_vars(self):
        return len(self.var_order)

    def __str__(self):
        o = self.var_order
        head = "Maximize" if self.sense == "max" else "Minimize"
        lines = [f"{head}  Z = {_expr_text(self.obj, o)}"
                 + (f" + {_num(self.obj_const)}" if self.obj_const else ""),
                 "subject to:"]
        for c in self.constraints:
            lines.append(f"    {c.name+':':<14}{c.as_text(o)}")
        for v in o:
            lo, up = self.lower.get(v, 0.0), self.upper.get(v, np.inf)
            bits = []
            if lo == -np.inf and up == np.inf:
                bits.append(f"{v} free")
            else:
                if lo != -np.inf:
                    bits.append(f"{v} >= {_num(lo)}")
                if up != np.inf:
                    bits.append(f"{v} <= {_num(up)}")
            lines.append("    " + ", ".join(bits))
        return "\n".join(lines)


def _num(x):
    if x == int(x):
        return str(int(x))
    return f"{x:g}"


def _expr_text(coeffs, order):
    parts = []
    for v in order:
        c = coeffs.get(v, 0.0)
        if c == 0:
            continue
        sign = "-" if c < 0 else ("+" if parts else "")
        a = abs(c)
        mag = "" if a == 1 else _num(a)
        parts.append(f"{sign} {mag}{v}".strip() if parts else f"{sign}{mag}{v}")
    return " ".join(parts) if parts else "0"
# ----------------------------------------------------------------------------
# 2. LP-file parser  (CPLEX LP format, the subset we need)
# ----------------------------------------------------------------------------

_SECTION = {
    "maximize": "max", "maximise": "max", "max": "max",
    "minimize": "min", "minimise": "min", "min": "min",
    "subject to": "st", "such that": "st", "st": "st", "s.t.": "st",
    "bounds": "bounds", "bound": "bounds",
    "free": "free",
    "end": "end",
}

_TERM = re.compile(r"""([+-]?)\s*           # sign
                       (\d*\.?\d*(?:[eE][+-]?\d+)?)\s*  # optional number
                       \*?\s*
                       ([A-Za-z_][A-Za-z0-9_.\[\]]*)     # variable name
                    """, re.X)
_OPS = re.compile(r"(<=|>=|=<|=>|<|>|=)")


def _parse_expr(s, model):
    """Turn '3 x1 - 2.5x2 + 7' into ({'x1':3,'x2':-2.5}, 7.0)."""
    coeffs, const = {}, 0.0
    s = s.strip()
    pos = 0
    token = re.compile(r"\s*([+-]?)\s*([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)?\s*"
                       r"\*?\s*([A-Za-z_][A-Za-z0-9_.\[\]]*)?")
    while pos < len(s):
        m = token.match(s, pos)
        if not m or m.end() == pos:
            break
        sign, num, var = m.group(1), m.group(2), m.group(3)
        if num is None and var is None:
            break
        mult = -1.0 if sign == "-" else 1.0
        val = float(num) if num is not None else 1.0
        if var:
            if var not in model.var_order:
                model.var_order.append(var)
            coeffs[var] = coeffs.get(var, 0.0) + mult * val
        else:
            const += mult * val
        pos = m.end()
    return coeffs, const


def parse_lp(text: str) -> LPModel:
    """Parse an LP file (string) in CPLEX LP format into an LPModel."""
    model = LPModel()
    # strip comments ( \ ... to end of line ) and blank lines
    raw = []
    for line in text.splitlines():
        line = line.split("\\")[0]
        if line.strip():
            raw.append(line.rstrip())

    section, buf, unnamed = None, [], itertools.count(1)

    def flush():
        """Process the statement accumulated in `buf`."""
        stmt = " ".join(buf).strip()
        buf.clear()
        if not stmt:
            return
        if section == "obj":
            name = model.obj_name
            if ":" in stmt.split("<")[0] and re.match(r"^\s*[A-Za-z_][\w.]*\s*:", stmt):
                name, stmt = stmt.split(":", 1)
                name = name.strip()
            model.obj_name = name
            c, k = _parse_expr(stmt, model)
            model.obj, model.obj_const = c, k
        elif section == "st":
            name = None
            if re.match(r"^\s*[A-Za-z_][\w.]*\s*:", stmt):
                name, stmt = stmt.split(":", 1)
                name = name.strip()
            m = _OPS.search(stmt)
            if not m:
                raise ValueError(f"Constraint without a relational operator: {stmt!r}")
            lhs, op, rhs = stmt[:m.start()], m.group(1), stmt[m.end():]
            op = {"=<": "<=", "=>": ">=", "<": "<=", ">": ">="}.get(op, op)
            lc, lk = _parse_expr(lhs, model)
            rc, rk = _parse_expr(rhs, model)
            coeffs = {v: lc.get(v, 0.0) - rc.get(v, 0.0)
                      for v in set(lc) | set(rc)}
            coeffs = {v: c for v, c in coeffs.items() if abs(c) > TOL}
            model.constraints.append(
                Constraint(name or f"c{next(unnamed)}", coeffs, op, rk - lk))
        elif section == "bounds":
            _parse_bound(stmt, model)
        elif section == "free":
            for v in re.findall(r"[A-Za-z_][A-Za-z0-9_.\[\]]*", stmt):
                if v not in model.var_order:
                    model.var_order.append(v)
                model.lower[v], model.upper[v] = -np.inf, np.inf

    for line in raw:
        key = line.strip().lower().rstrip(":")
        if key in _SECTION:
            flush()
            tag = _SECTION[key]
            if tag in ("max", "min"):
                model.sense, section = tag, "obj"
            elif tag == "end":
                flush()
                section = None
                break
            else:
                section = tag
            continue
        if section is None:
            continue
        # a new statement starts when the previous one is complete
        if section in ("st", "bounds", "free") and buf:
            flush()
        buf.append(line)
    flush()

    # default bounds: x >= 0
    for v in model.var_order:
        model.lower.setdefault(v, 0.0)
        model.upper.setdefault(v, np.inf)
    return model


def _parse_bound(stmt, model):
    s = stmt.strip()
    if re.match(r"^[A-Za-z_][\w.\[\]]*\s+free$", s, re.I):
        v = s.split()[0]
        if v not in model.var_order:
            model.var_order.append(v)
        model.lower[v], model.upper[v] = -np.inf, np.inf
        return
    parts = _OPS.split(s)
    parts = [p.strip() for p in parts]
    def norm(o):
        return {"=<": "<=", "=>": ">=", "<": "<=", ">": ">="}.get(o, o)

    if len(parts) == 5:                      # lo <= x <= up
        lo, o1, v, o2, up = parts
        if v not in model.var_order:
            model.var_order.append(v)
        model.lower[v] = _val(lo)
        model.upper[v] = _val(up)
    elif len(parts) == 3:
        a, op, b = parts[0], norm(parts[1]), parts[2]
        if re.match(r"^[A-Za-z_]", a):
            v, val = a, _val(b)
        else:
            v, val = b, _val(a)
            op = {"<=": ">=", ">=": "<=", "=": "="}[op]
        if v not in model.var_order:
            model.var_order.append(v)
        if op == "<=":
            model.upper[v] = val
        elif op == ">=":
            model.lower[v] = val
        else:
            model.lower[v] = model.upper[v] = val


def _val(tok):
    tok = tok.strip().lower().replace(" ", "")
    if tok in ("inf", "+inf", "infinity", "+infinity"):
        return np.inf
    if tok in ("-inf", "-infinity"):
        return -np.inf
    return float(tok)


def load_lp(path) -> LPModel:
    with open(path) as f:
        return parse_lp(f.read())
# ----------------------------------------------------------------------------
# 3. Half-plane representation + corner-point solver
# ----------------------------------------------------------------------------

@dataclass
class HalfPlane:
    a: float
    b: float
    op: str      # '<=', '>=', '='
    c: float
    label: str

    def slack(self, x, y):
        return self.c - (self.a * x + self.b * y)

    def holds(self, x, y, tol=1e-6):
        lhs = self.a * x + self.b * y
        if self.op == "<=":
            return lhs <= self.c + tol
        if self.op == ">=":
            return lhs >= self.c - tol
        return abs(lhs - self.c) <= tol


def to_halfplanes(model: LPModel):
    """All constraints + bounds of a 2-variable model as half-planes."""
    if model.n_vars != 2:
        raise ValueError(
            f"The graphical method needs exactly 2 decision variables, "
            f"this model has {model.n_vars}: {model.var_order}")
    v1, v2 = model.var_order
    hps = [HalfPlane(c.coeffs.get(v1, 0.0), c.coeffs.get(v2, 0.0),
                     c.op, c.rhs, c.name) for c in model.constraints]
    for i, v in enumerate((v1, v2)):
        a, b = (1.0, 0.0) if i == 0 else (0.0, 1.0)
        lo, up = model.lower[v], model.upper[v]
        if lo != -np.inf:
            hps.append(HalfPlane(a, b, ">=", lo, f"{v} >= {_num(lo)}"))
        if up != np.inf:
            hps.append(HalfPlane(a, b, "<=", up, f"{v} <= {_num(up)}"))
    return hps


def _intersect(h1, h2):
    A = np.array([[h1.a, h1.b], [h2.a, h2.b]], float)
    if abs(np.linalg.det(A)) < 1e-12:
        return None
    return np.linalg.solve(A, np.array([h1.c, h2.c], float))


def corner_points(model: LPModel, tol=1e-6):
    """Every vertex of the feasible region, with the constraints binding there."""
    hps = to_halfplanes(model)
    pts = []
    for h1, h2 in itertools.combinations(hps, 2):
        p = _intersect(h1, h2)
        if p is None:
            continue
        if all(h.holds(p[0], p[1], tol) for h in hps):
            for q, _ in pts:
                if abs(q[0] - p[0]) < 1e-6 and abs(q[1] - p[1]) < 1e-6:
                    break
            else:
                binding = [h.label for h in hps
                           if abs(h.a * p[0] + h.b * p[1] - h.c) < 1e-6]
                pts.append((p, binding))
    return pts


def objective_value(model, x, y):
    v1, v2 = model.var_order
    return (model.obj.get(v1, 0.0) * x + model.obj.get(v2, 0.0) * y
            + model.obj_const)


@dataclass
class LPResult:
    status: str                       # 'optimal', 'infeasible', 'unbounded'
    model: LPModel = None
    x: tuple = None
    z: float = None
    corners: list = field(default_factory=list)     # [(point, binding, z)]
    optimal_corners: list = field(default_factory=list)
    multiple: bool = False
    message: str = ""

    def __str__(self):
        if self.status != "optimal":
            return f"Status: {self.status.upper()} -- {self.message}"
        v1, v2 = self.model.var_order
        s = [f"Status : OPTIMAL",
             f"{v1:<6}= {self.x[0]:.6g}",
             f"{v2:<6}= {self.x[1]:.6g}",
             f"Z      = {self.z:.6g}"]
        if self.multiple:
            alt = ", ".join(f"({p[0]:.4g}, {p[1]:.4g})"
                            for p, _, _ in self.optimal_corners)
            s.append(f"\nMultiple optimal solutions: every point on the segment "
                     f"between {alt} is optimal.")
        return "\n".join(s)

    def corner_table(self):
        try:
            import pandas as pd
        except ImportError:
            return None
        v1, v2 = self.model.var_order
        rows = [{v1: round(p[0], 6), v2: round(p[1], 6), "Z": round(z, 6),
                 "binding constraints": ", ".join(b), "optimal": "  <== " if
                 any(p is q for q, _, _ in self.optimal_corners) else ""}
                for p, b, z in self.corners]
        return pd.DataFrame(rows).sort_values("Z",
                ascending=(self.model.sense == "min")).reset_index(drop=True)


def solve_lp(model: LPModel, tol=1e-6) -> LPResult:
    """Corner-point (extreme-point) solution of a 2-variable LP."""
    hps = to_halfplanes(model)
    raw = corner_points(model, tol)
    corners = [(p, b, objective_value(model, p[0], p[1])) for p, b in raw]

    feasible = _feasible_somewhere(model, hps, corners)
    if not feasible:
        return LPResult("infeasible", model,
                        message="every point violates at least one constraint, "
                                "the feasible region is empty.")

    unbounded, ray = _unbounded_direction(model, hps)
    if unbounded:
        return LPResult("unbounded", model, corners=corners,
                        message=f"the objective improves without limit along the "
                                f"direction {np.round(ray, 4).tolist()}.")

    if not corners:
        return LPResult("unbounded", model,
                        message="the feasible region has no extreme point "
                                "(unbounded in both directions).")

    best = max(c[2] for c in corners) if model.sense == "max" \
        else min(c[2] for c in corners)
    opt = [c for c in corners if abs(c[2] - best) < 1e-7]
    return LPResult("optimal", model, x=tuple(opt[0][0]), z=best,
                    corners=corners, optimal_corners=opt,
                    multiple=len(opt) > 1)


def _feasible_somewhere(model, hps, corners):
    if corners:
        return True
    # fall back on a linear-programming feasibility check
    from scipy.optimize import linprog
    A_ub, b_ub, A_eq, b_eq = [], [], [], []
    for h in hps:
        if h.op == "<=":
            A_ub.append([h.a, h.b]); b_ub.append(h.c)
        elif h.op == ">=":
            A_ub.append([-h.a, -h.b]); b_ub.append(-h.c)
        else:
            A_eq.append([h.a, h.b]); b_eq.append(h.c)
    r = linprog(c=[0, 0],
                A_ub=A_ub or None, b_ub=b_ub or None,
                A_eq=A_eq or None, b_eq=b_eq or None,
                bounds=[(None, None), (None, None)], method="highs")
    return r.status == 0


def _unbounded_direction(model, hps):
    """Is there a feasible ray along which the objective keeps improving?"""
    from scipy.optimize import linprog
    v1, v2 = model.var_order
    c = np.array([model.obj.get(v1, 0.0), model.obj.get(v2, 0.0)], float)
    if not np.any(c):
        return False, None
    # maximise c.d over the recession cone, |d| bounded by a box
    A_ub, b_ub, A_eq, b_eq = [], [], [], []
    for h in hps:
        if h.op == "<=":
            A_ub.append([h.a, h.b]); b_ub.append(0.0)
        elif h.op == ">=":
            A_ub.append([-h.a, -h.b]); b_ub.append(0.0)
        else:
            A_eq.append([h.a, h.b]); b_eq.append(0.0)
    obj = -c if model.sense == "max" else c
    r = linprog(c=obj, A_ub=A_ub or None, b_ub=b_ub or None,
                A_eq=A_eq or None, b_eq=b_eq or None,
                bounds=[(-1, 1), (-1, 1)], method="highs")
    if r.status == 0:
        gain = c @ r.x
        improves = gain > 1e-7 if model.sense == "max" else gain < -1e-7
        if improves:
            return True, r.x
    return False, None
# ----------------------------------------------------------------------------
# 4. Plotting
# ----------------------------------------------------------------------------

def _plot_window(model, corners, pad=0.25, unbounded=False):
    xs, ys = [0.0], [0.0]
    for p, _, _ in corners:
        xs.append(p[0]); ys.append(p[1])
    for h in to_halfplanes(model):
        if abs(h.b) > TOL:
            ys.append(h.c / h.b)
        if abs(h.a) > TOL:
            xs.append(h.c / h.a)
    xs = [v for v in xs if np.isfinite(v)]
    ys = [v for v in ys if np.isfinite(v)]
    xlo, xhi = min(xs), max(xs)
    ylo, yhi = min(ys), max(ys)
    if unbounded:
        pad = max(pad, 0.9)
    dx = max(xhi - xlo, 1.0) * pad
    dy = max(yhi - ylo, 1.0) * pad
    return (xlo - dx, xhi + dx, ylo - dy, yhi + dy)


def plot_lp(model: LPModel, result: LPResult = None, iso_values=None,
            ax=None, show_corners=True, shade=True, title=None,
            window=None, figsize=(8, 6.5), grid_n=600):
    """Draw constraint lines, the feasible region, corner points and iso-lines."""
    if result is None:
        result = solve_lp(model)
    v1, v2 = model.var_order
    hps = to_halfplanes(model)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    x0, x1, y0, y1 = window or _plot_window(
        model, result.corners, unbounded=(result.status == "unbounded"))
    x0, y0 = min(x0, 0), min(y0, 0)

    # --- feasible region as a pixel mask (works for bounded and unbounded) ---
    if shade:
        gx = np.linspace(x0, x1, grid_n)
        gy = np.linspace(y0, y1, grid_n)
        X, Y = np.meshgrid(gx, gy)
        mask = np.ones_like(X, dtype=bool)
        scale = max(x1 - x0, y1 - y0)
        for h in hps:
            lhs = h.a * X + h.b * Y
            eps = 1e-9 * max(1.0, abs(h.c))
            if h.op == "<=":
                mask &= lhs <= h.c + eps
            elif h.op == ">=":
                mask &= lhs >= h.c - eps
            else:
                mask &= np.abs(lhs - h.c) <= 5e-3 * scale
        if mask.any():
            ax.contourf(X, Y, mask.astype(float), levels=[0.5, 1.5],
                        colors=["#7ba7d7"], alpha=0.35, zorder=0)

        # a feasible set that is only a segment or a point is hard to see
        degenerate = (mask.sum() < 0.01 * mask.size
                      or any(h.op == "=" for h in hps))
        if degenerate and len(result.corners) == 2:
            (p, q) = [c[0] for c in result.corners]
            ax.plot([p[0], q[0]], [p[1], q[1]], color="#1f4e79", lw=6,
                    alpha=0.45, solid_capstyle="round", zorder=1)
        elif len(result.corners) == 1 and not mask.any():
            ax.plot(*result.corners[0][0], "o", ms=14, color="#7ba7d7",
                    alpha=0.6, zorder=1)

    # --- constraint lines ----------------------------------------------------
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    for i, c in enumerate(model.constraints):
        a, b = c.coeffs.get(v1, 0.0), c.coeffs.get(v2, 0.0)
        col = colors[i % 10]
        lbl = f"{c.name}: {c.as_text(model.var_order)}"
        ls = "-" if c.op != "=" else "-"
        if abs(b) > TOL:
            xx = np.array([x0, x1])
            ax.plot(xx, (c.rhs - a * xx) / b, color=col, lw=2, ls=ls, label=lbl)
        elif abs(a) > TOL:
            ax.axvline(c.rhs / a, color=col, lw=2, ls=ls, label=lbl)

    # --- axes / non-negativity ----------------------------------------------
    ax.axhline(0, color="0.35", lw=1.2)
    ax.axvline(0, color="0.35", lw=1.2)

    # --- iso-objective lines -------------------------------------------------
    c1, c2 = model.obj.get(v1, 0.0), model.obj.get(v2, 0.0)
    if iso_values is None and (abs(c1) > TOL or abs(c2) > TOL):
        # space the iso-lines so that they are visibly apart on this window
        step = 0.26 * np.hypot(c1, c2) * max(x1 - x0, y1 - y0)
        base = result.z if result.status == "optimal" else \
            objective_value(model, (x0 + x1) / 2, (y0 + y1) / 2)
        sgn = -1.0 if model.sense == "max" else 1.0
        iso_values = [base + sgn * k * step for k in range(4)]
    if iso_values and (abs(c1) > TOL or abs(c2) > TOL):
        for k, z in enumerate(iso_values):
            zz = z - model.obj_const
            lab = "iso-objective lines" if k == 0 else None
            if abs(c2) > TOL:
                xx = np.array([x0, x1])
                ax.plot(xx, (zz - c1 * xx) / c2, color="0.45", lw=1,
                        ls="--", alpha=0.8, label=lab)
            elif abs(c1) > TOL:
                ax.axvline(zz / c1, color="0.45", lw=1, ls="--",
                           alpha=0.8, label=lab)

    # --- corner points -------------------------------------------------------
    if show_corners:
        for p, _, z in result.corners:
            ax.plot(*p, "o", ms=7, mfc="white", mec="#12365c", mew=1.8, zorder=5)
            ax.annotate(f"({p[0]:.4g}, {p[1]:.4g})\nZ={z:.5g}", p,
                        textcoords="offset points", xytext=(7, 7),
                        fontsize=8, color="#12365c")
    if result.status == "optimal":
        for p, _, _ in result.optimal_corners:
            ax.plot(*p, "*", ms=20, color="#d1495b", zorder=6)

    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_xlabel(v1); ax.set_ylabel(v2)
    head = "Maximize" if model.sense == "max" else "Minimize"
    ax.set_title(title or f"{head} Z = {_expr_text(model.obj, model.var_order)}"
                          f"   -- {result.status.upper()}")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(alpha=0.25, ls=":")
    return ax


def solve_and_plot(source, **kw):
    """`source` may be a path to an .lp file, LP text, or an LPModel."""
    if isinstance(source, LPModel):
        model = source
    elif isinstance(source, str) and ("\n" in source or not source.endswith(".lp")):
        model = parse_lp(source)
    else:
        model = load_lp(source)
    res = solve_lp(model)
    print(model, end="\n\n")
    print(res)
    tbl = res.corner_table()
    if tbl is not None and len(tbl):
        print("\nCorner points:")
        print(tbl.to_string(index=False))
    plot_lp(model, res, **kw)
    plt.show()
    return model, res
