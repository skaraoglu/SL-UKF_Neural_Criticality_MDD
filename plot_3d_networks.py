"""
plot_3d_networks.py
───────────────────────────────────────────────────────────────────────────────
3D interactive brain network visualization using Plotly.
Called from main_analysis.ipynb (R) via system2() or reticulate.

Reads two CSVs written by Section 7 of main_analysis.ipynb:
  results/network_edges.csv  — all edges with weights and MNI endpoints
  results/atlas_coords.csv   — node MNI (x,y,z) and lobe labels

Writes:
  results/figures/network_3d_w_k.html
  results/figures/network_3d_w_norm.html
  results/figures/network_3d_w_sync.html
  results/figures/network_3d_w_frac.html
  results/figures/network_3d_all.html   (4-panel combined)
"""

import argparse, os, sys
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── CLI ───────────────────────────────────────────────────────────────────────
par = argparse.ArgumentParser()
par.add_argument("--data_dir", default="results")
par.add_argument("--out_dir",  default="results/figures")
par.add_argument("--top_n",    type=int, default=200)
args, _ = par.parse_known_args()

DATA_DIR, OUT_DIR, TOP_N = args.data_dir, args.out_dir, args.top_n
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
for fname in ["network_edges.csv", "atlas_coords.csv"]:
    fpath = os.path.join(DATA_DIR, fname)
    if not os.path.exists(fpath):
        print(f"ERROR: {fpath} not found. Run Section 7 of main_analysis.ipynb first.")
        sys.exit(1)

edges = pd.read_csv(os.path.join(DATA_DIR, "network_edges.csv"))
nodes = pd.read_csv(os.path.join(DATA_DIR, "atlas_coords.csv"))
print(f"Loaded {len(edges)} edges, {len(nodes)} nodes.")

# ── Lobe colour palette — mirrors LOBE_COLOURS in R/atlas.R ──────────────────
LOBE_COL = {
    "Frontal":   "#E74C3C",
    "Temporal":  "#3498DB",
    "Parietal":  "#2ECC71",
    "Occipital": "#F39C12",
    "Cingulate": "#9B59B6",
    "Insula":    "#1ABC9C",
    "SCGM":      "#95A5A6",
}

# ── Metric descriptors ────────────────────────────────────────────────────────
METRICS = [
    dict(col="w_k",    log_col="w_k_log",
         title="Raw k",
         subtitle="Absolute spring coupling constant",
         colour="#2C3E50"),
    dict(col="w_norm", log_col="w_norm_log",
         title="k / sqrt(ab)  (Normalized)",
         subtitle="Coupling relative to intrinsic oscillator dynamics",
         colour="#2980B9"),
    dict(col="w_sync", log_col="w_sync_log",
         title="k / |a-b|  (Synchronizability)",
         subtitle="Proximity to synchronization -- Arnold tongue measure",
         colour="#27AE60"),
    dict(col="w_frac", log_col="w_frac_log",
         title="2k/(a+b+2k)  (Coupling Fraction)",
         subtitle="Fraction of total system stiffness due to coupling; bounded (0,1)",
         colour="#8E44AD"),
]


def hex_rgba(h, a):
    r, g, b = int(h[1:3],16), int(h[3:5],16), int(h[5:7],16)
    return f"rgba({r},{g},{b},{a:.3f})"


def build_traces(edf, ndf, metric, top_n):
    """Scatter3d traces: binned edge segments + lobe-grouped node markers."""
    col, lcol, colour = metric["col"], metric["log_col"], metric["colour"]

    df      = edf.nlargest(top_n, col).copy()
    lw      = df[lcol].values
    lw_norm = (lw - lw.min()) / (lw.max() - lw.min() + 1e-12)

    traces  = []
    N_BINS  = 12
    bins    = np.floor(lw_norm * N_BINS).astype(int).clip(0, N_BINS-1)

    for b in range(N_BINS):
        mask = bins == b
        if not mask.any():
            continue
        frac    = b / (N_BINS - 1)
        alpha   = 0.04 + 0.88 * frac
        width   = 0.4  + 3.5  * frac
        sub     = df[mask]
        xs, ys, zs = [], [], []
        for _, row in sub.iterrows():
            xs += [row.xA, row.xB, None]
            ys += [row.yA, row.yB, None]
            zs += [row.zA, row.zB, None]
        traces.append(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="lines",
            line=dict(color=hex_rgba(colour, alpha), width=width),
            hoverinfo="skip", showlegend=False))

    # ── Node strength from full edge list (before top_n filter) ─────────────
    str_A = edf.groupby("ROI_A")[col].sum().reset_index().rename(
        columns={"ROI_A": "name", col: "strength"})
    str_B = edf.groupby("ROI_B")[col].sum().reset_index().rename(
        columns={"ROI_B": "name", col: "strength"})
    strength_df = (pd.concat([str_A, str_B])
                     .groupby("name", as_index=False)["strength"].sum())

    active    = pd.concat([df["ROI_A"], df["ROI_B"]]).unique()
    sub_nodes = ndf[ndf["name"].isin(active)].copy().reset_index(drop=True)
    sub_nodes = sub_nodes.merge(strength_df, on="name", how="left")
    sub_nodes["strength"] = sub_nodes["strength"].fillna(0)

    s_min = sub_nodes["strength"].min()
    s_max = sub_nodes["strength"].max()
    # Map strength → marker size [5, 18]
    sub_nodes["marker_size"] = (5 + 13 * (sub_nodes["strength"] - s_min)
                                  / (s_max - s_min + 1e-12))

    for lobe, grp in sub_nodes.groupby("lobe"):
        c = LOBE_COL.get(lobe, "#aaa")
        traces.append(go.Scatter3d(
            x=grp.x.tolist(), y=grp.y.tolist(), z=grp.z.tolist(),
            mode="markers+text",
            text=grp.name.tolist(),
            textposition="top center",
            textfont=dict(size=7, color="rgba(30,30,30,0.75)"),
            marker=dict(
                size=grp["marker_size"].tolist(),
                color="white",
                line=dict(color=c, width=2.2)
            ),
            name=lobe, legendgroup=lobe,
            hovertemplate=(
                "<b>%{text}</b><br>"
                "x=%{x:.1f}  y=%{y:.1f}  z=%{z:.1f}<br>"
                "strength=%{customdata:.3g}<extra></extra>"
            ),
            customdata=grp["strength"].tolist()
        ))
    return traces


SCENE = dict(
    xaxis=dict(title="x MNI  (L<->R)", showgrid=False, zeroline=False,
               showbackground=True, backgroundcolor="rgba(240,242,248,0.9)"),
    yaxis=dict(title="y MNI  (P<->A)", showgrid=False, zeroline=False,
               showbackground=True, backgroundcolor="rgba(240,242,248,0.9)"),
    zaxis=dict(title="z MNI  (I<->S)", showgrid=False, zeroline=False,
               showbackground=True, backgroundcolor="rgba(240,242,248,0.9)"),
    camera=dict(eye=dict(x=1.45, y=1.45, z=0.8)),
    aspectmode="data")

# ── Individual plots ──────────────────────────────────────────────────────────
print(f"\nBuilding individual 3D network plots (top {TOP_N} edges each)...")
for m in METRICS:
    traces = build_traces(edges, nodes, m, TOP_N)
    fig    = go.Figure(data=traces)
    fig.update_layout(
        title=dict(
            text=(f"<b>Brain Network -- {m['title']}</b><br>"
                  f"<span style='font-size:11px'>{m['subtitle']}"
                  f" | Top {TOP_N} edges | HOA 110 ROIs | MNI 3D</span>"),
            x=0.5, xanchor="center", font=dict(size=14)),
        scene=SCENE,
        legend=dict(title="Lobe", x=1.0, y=0.92,
                    bgcolor="rgba(255,255,255,0.88)",
                    bordercolor="#ccc", borderwidth=1),
        margin=dict(l=0, r=0, t=85, b=0),
        paper_bgcolor="white", width=950, height=720)
    out = os.path.join(OUT_DIR, f"network_3d_{m['col']}.html")
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"  Saved: {out}")

# ── 4-panel combined ──────────────────────────────────────────────────────────
print("\nBuilding 4-panel combined figure...")
fig4 = make_subplots(
    rows=2, cols=2,
    specs=[[{"type":"scene"},{"type":"scene"}],
           [{"type":"scene"},{"type":"scene"}]],
    subplot_titles=[m["title"] for m in METRICS],
    horizontal_spacing=0.01, vertical_spacing=0.05)

scenes   = ["scene","scene2","scene3","scene4"]
seen_lob = set()

for idx, m in enumerate(METRICS):
    row = idx // 2 + 1
    col = idx  % 2 + 1
    for tr in build_traces(edges, nodes, m, TOP_N):
        is_node = hasattr(tr, "name") and tr.name in LOBE_COL
        if is_node:
            tr.showlegend  = tr.name not in seen_lob
            tr.legendgroup = tr.name
            if tr.showlegend:
                seen_lob.add(tr.name)
        else:
            tr.showlegend = False
        fig4.add_trace(tr, row=row, col=col)

    panel_scene = dict(
        xaxis=dict(title="", showgrid=False, zeroline=False,
                   showticklabels=False, showbackground=True,
                   backgroundcolor="rgba(240,242,248,0.9)"),
        yaxis=dict(title="", showgrid=False, zeroline=False,
                   showticklabels=False, showbackground=True,
                   backgroundcolor="rgba(240,242,248,0.9)"),
        zaxis=dict(title="", showgrid=False, zeroline=False,
                   showticklabels=False, showbackground=True,
                   backgroundcolor="rgba(240,242,248,0.9)"),
        camera=dict(eye=dict(x=1.5, y=1.5, z=0.75)),
        aspectmode="data")
    fig4.update_layout(**{scenes[idx]: panel_scene})

fig4.update_layout(
    title=dict(
        text=("<b>3D Brain Coupling Networks -- All Four Metrics</b><br>"
              f"<span style='font-size:11px'>Top {TOP_N} edges per metric | "
              "Edge opacity and width proportional to log1p(weight) | HOA 110 ROIs</span>"),
        x=0.5, xanchor="center", font=dict(size=15)),
    legend=dict(title="Lobe", x=1.0, y=0.95,
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="#ccc", borderwidth=1),
    margin=dict(l=0, r=130, t=90, b=0),
    paper_bgcolor="white", width=1300, height=960)

out4 = os.path.join(OUT_DIR, "network_3d_all.html")
fig4.write_html(out4, include_plotlyjs="cdn")
print(f"  Saved: {out4}")
print("\nAll 3D figures written. Open any .html file in a browser for interactive viewing.")
