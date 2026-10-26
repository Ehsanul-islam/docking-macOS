"""2D protein-ligand interaction image (RDKit ligand depiction, residue discs placed by ray-casting then relaxed,
coloured dashed lines by interaction type, legend). The `reference` style reproduces the supplied manuscript figures;
`mono` (print-friendly greyscale) and `dark` use the same layout. Input: a profile dict from interactions.profile()."""
import io
import json
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

_REF = {  # type: (disc face, line colour, line style, priority)
    "Salt bridge":            ("#F5A623", "#E07B00", (0, (4, 2)), 0),
    "Hydrogen bond":          ("#6FCF6F", "#1E8C1E", (0, (4, 2)), 1),
    "Water bridge":           ("#8FD6F0", "#1C8FC4", (0, (2, 2)), 2),
    "π-stacking (P)":         ("#E27FD8", "#B0169E", (0, (5, 2)), 3),
    "π-stacking (T)":         ("#E27FD8", "#B0169E", (0, (5, 2)), 3),
    "π-cation":               ("#F2B870", "#C46A00", (0, (5, 2)), 4),
    "Hydrophobic":            ("#F4B6C8", "#C2587A", (0, (1, 1.6)), 5),
    "van der Waals":          ("#DDF2D4", None, None, 6),
}
_GREY = {  # print-friendly: shades of grey, line styles tell the types apart
    "Salt bridge":            ("#BDBDBD", "#000000", (0, (4, 1.5)), 0),
    "Hydrogen bond":          ("#8C8C8C", "#000000", (0, (4, 2)), 1),
    "Water bridge":           ("#D9D9D9", "#333333", (0, (2, 2)), 2),
    "π-stacking (P)":         ("#6E6E6E", "#000000", (0, (7, 2, 1, 2)), 3),
    "π-stacking (T)":         ("#6E6E6E", "#000000", (0, (7, 2, 1, 2)), 3),
    "π-cation":               ("#A6A6A6", "#222222", (0, (5, 2)), 4),
    "Hydrophobic":            ("#E6E6E6", "#555555", (0, (1, 1.6)), 5),
    "van der Waals":          ("#F6F6F6", None, None, 6),
}
_DARK = {
    "Salt bridge":            ("#C98A2E", "#FFB84D", (0, (4, 2)), 0),
    "Hydrogen bond":          ("#3F8F3F", "#6FE06F", (0, (4, 2)), 1),
    "Water bridge":           ("#3B7C94", "#6CC6EA", (0, (2, 2)), 2),
    "π-stacking (P)":         ("#9A3F94", "#F077E6", (0, (5, 2)), 3),
    "π-stacking (T)":         ("#9A3F94", "#F077E6", (0, (5, 2)), 3),
    "π-cation":               ("#A86B2A", "#FFB066", (0, (5, 2)), 4),
    "Hydrophobic":            ("#8E4A62", "#FF8FB0", (0, (1, 1.6)), 5),
    "van der Waals":          ("#3A4A3A", None, None, 6),
}
THEMES = {
    "reference": dict(styles=_REF, bg="white", fg="black", edge="#555555", legend_bg="white", palette=None),
    "mono": dict(styles=_GREY, bg="white", fg="black", edge="#444444", legend_bg="white", palette="bw"),
    "dark": dict(styles=_DARK, bg="#161a21", fg="#F2F2F2", edge="#9AA0AA", legend_bg="#1f242d", palette="dark"),
}
LEGEND_NAME = {"π-stacking (P)": "π–π stacking (parallel)", "π-stacking (T)": "π–π stacking (T-shaped)",
               "Hydrophobic": "Hydrophobic (alkyl / π-alkyl)", "van der Waals": "van der Waals contact (≤ 4 Å)",
               "Water bridge": "Water-mediated H-bond"}

W, H = 1500, 1150
S = 4  # supersample the RDKit depiction so the molecule stays sharp at 600 dpi


def draw(profile, png, style="reference", dpi=600):
    """Render `profile` (dict from interactions.profile, or a path to its JSON) to `png`. Returns the number of residues."""
    d = json.load(open(profile)) if isinstance(profile, (str, pathlib.Path)) else profile
    th = THEMES[style]
    ST, BG, FG, EDGE = th["styles"], th["bg"], th["fg"], th["edge"]
    mol3d = (Chem.MolFromMolBlock(d["ligand_molblock"], removeHs=False) if "ligand_molblock"
             in d else Chem.MolFromMolFile(d["ligand_sdf"], removeHs=False))
    heavy = [a.GetIdx() for a in mol3d.GetAtoms() if a.GetAtomicNum() > 1]
    xyz = np.array([list(mol3d.GetConformer().GetAtomPosition(i)) for i in heavy])
    mol = Chem.RemoveHs(mol3d)                        # keeps heavy-atom order
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(W * S, H * S)
    o = drawer.drawOptions()
    o.padding = 0.30 if mol.GetNumAtoms() < 30 else 0.27
    o.bondLineWidth = 3 * S; o.minFontSize = 20 * S; o.maxFontSize = 26 * S; o.additionalAtomLabelPadding = 0.1
    o.clearBackground = True
    bgc = matplotlib.colors.to_rgba(BG)
    o.setBackgroundColour(bgc)
    if th["palette"] == "bw":
        o.useBWAtomPalette()
    elif th["palette"] == "dark":
        o.setAtomPalette({6: (0.93, 0.93, 0.93), 7: (0.45, 0.62, 1.0), 8: (1.0, 0.42, 0.42), 16: (1.0, 0.85, 0.3),
                          1: (0.93, 0.93, 0.93), 15: (1.0, 0.6, 0.2), 9: (0.5, 1, 0.7), 17: (0.5, 1, 0.5)})
    drawer.DrawMolecule(mol); drawer.FinishDrawing()
    img = Image.open(io.BytesIO(drawer.GetDrawingText()))
    pix = np.array([[drawer.GetDrawCoords(i).x, drawer.GetDrawCoords(i).y] for i in range(mol.GetNumAtoms())]) / S

    def anchor(p3):
        dd = np.linalg.norm(xyz - np.array(p3), axis=1)
        return pix[int(np.argmin(dd))]

    def ring_anchor(ring):
        idx = [int(np.argmin(np.linalg.norm(xyz - np.array(a), axis=1))) for a in ring]
        return pix[idx].mean(0)

    res = {}
    S_ = S
    for it in d["interactions"]:
        r = res.setdefault(it["residue"], dict(lines=[], types=set()))
        a = ring_anchor(it["lig_ring"]) if "lig_ring" in it else anchor(it["lig"])
        r["types"].add(it["type"]); r["lines"].append((it["type"], a, it["dist"], it.get("water_name")))
    for c in d["vdw_contacts"]:
        if c["residue"] not in res:
            res[c["residue"]] = dict(lines=[], types={"van der Waals"}, vdw_anchor=anchor(c["lig"]))
    for r in res.values():            # closest hydrophobic line per residue only
        hyd = [l for l in r["lines"] if l[0] == "Hydrophobic"]
        if len(hyd) > 1:
            best = min(hyd, key=lambda l: l[2])
            r["lines"] = [l for l in r["lines"] if l[0] != "Hydrophobic"] + [best]
        uniq = {}
        for l in r["lines"]:
            k = (l[0], round(l[1][0]), round(l[1][1]))
            if k not in uniq or l[2] < uniq[k][2]: uniq[k] = l
        r["lines"] = list(uniq.values())
        r["main"] = min(r["types"], key=lambda t: ST[t][3])

    R_DISC, GAP = 48, 30
    cen = pix.mean(0)
    obst = [pix] + [np.array([(pix[b.GetBeginAtomIdx()] + pix[b.GetEndAtomIdx()]) / 2 for b in mol.GetBonds()])]
    rings = [pix[list(r)].mean(0) for r in mol.GetRingInfo().AtomRings()]
    if rings: obst.append(np.array(rings))
    OB = np.vstack(obst)
    names = list(res)
    target = {}
    for n in names:
        r = res[n]
        anc = [l[1] for l in r["lines"]] or [r["vdw_anchor"]]
        a = np.mean(anc, axis=0)
        u = a - cen; u = u / (np.linalg.norm(u) + 1e-6)
        if np.linalg.norm(a - cen) < 1: u = np.array([0.0, -1.0])
        t = 0.0
        extra = 0 if r["lines"] else 25
        while np.linalg.norm(OB - (a + u * t), axis=1).min() < R_DISC + GAP + extra and t < 1200:
            t += 4
        target[n] = a + u * t
        res[n]["anchor"] = a
    rng = np.random.default_rng(7)
    P = {n: target[n] + rng.normal(0, 25, 2) for n in names}
    MIN = 2 * R_DISC + 14
    for it in range(2500):
        F = {n: 0.03 * (target[n] - P[n]) for n in names}
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                v = P[a] - P[b]; dd = np.linalg.norm(v)
                if dd < 1e-3: v = rng.normal(0, 1, 2); dd = np.linalg.norm(v)
                if dd < MIN:
                    push = 0.5 * (MIN - dd) * v / dd
                    F[a] += push; F[b] -= push
            dv = P[a] - OB; dist = np.linalg.norm(dv, axis=1) + 1e-6
            close = dist < R_DISC + GAP
            if close.any():
                F[a] += ((R_DISC + GAP - dist[close])[:, None] * dv[close] / dist[close][:, None]).sum(0) * 0.6
        for n in names: P[n] += np.clip(F[n], -20, 20)
    for it in range(400):
        moved = False
        for a in names:
            dv = P[a] - OB; dist = np.linalg.norm(dv, axis=1); j = int(np.argmin(dist))
            if dist[j] < R_DISC + GAP - 4:
                P[a] = P[a] + dv[j] / (dist[j] + 1e-6) * (R_DISC + GAP - dist[j]); moved = True
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                v = P[a] - P[b]; dd = np.linalg.norm(v) + 1e-6
                if dd < MIN - 2:
                    P[a] = P[a] + 0.5 * (MIN - dd) * v / dd; P[b] = P[b] - 0.5 * (MIN - dd) * v / dd; moved = True
        for a in names:
            for b in names:
                if a == b: continue
                for l in res[b]["lines"]:
                    p0, p1 = P[b], l[1]
                    seg = p1 - p0; L2 = seg @ seg
                    if L2 < 1: continue
                    tt = np.clip((P[a] - p0) @ seg / L2, 0, 1); q = p0 + tt * seg
                    dv = P[a] - q; dd = np.linalg.norm(dv)
                    if dd < R_DISC + 6 and 0.05 < tt < 0.95:
                        nrm = dv / dd if dd > 1e-3 else np.array([-seg[1], seg[0]]) / np.sqrt(L2)
                        P[a] = P[a] + nrm * (R_DISC + 8 - dd) * 0.5; moved = True
        if not moved: break
    for n in names: res[n]["xy"] = P[n]

    allxy = np.vstack([pix] + [P[n] for n in names]) if names else pix
    x0, y0 = allxy.min(0) - R_DISC - 40; x1, y1 = allxy.max(0) + R_DISC + 40
    y0 -= 90
    bw, bh = x1 - x0, y1 - y0
    LEG = 1.5
    fig = plt.figure(figsize=(bw / 100, bh / 100 + LEG), dpi=150, facecolor=BG)
    ax = fig.add_axes([0, LEG / (bh / 100 + LEG), 1, (bh / 100) / (bh / 100 + LEG)])
    ax.imshow(img, extent=[0, W, H, 0], zorder=1)
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis("off")
    placed = []
    for n, r in res.items():
        for t, a, dist, wname in r["lines"]:
            face, lc, ls, _ = ST[t]
            p = r["xy"]
            v = a - p; Ln = np.linalg.norm(v); start = p + v / Ln * R_DISC
            ax.plot([start[0], a[0]], [start[1], a[1]], color=lc, lw=2.4, ls=ls, zorder=2, solid_capstyle="round")
            if t != "Hydrophobic":
                for frac in (0.28, 0.45, 0.62, 0.78, 0.2):
                    m = start + (a - start) * frac
                    if all(np.linalg.norm(m - q) > 42 for q in placed) and \
                       all(np.linalg.norm(m - rr["xy"]) > R_DISC + 12 for rr in res.values()):
                        break
                seg = a - start; L = np.linalg.norm(seg) + 1e-6; nrm = np.array([-seg[1], seg[0]]) / L
                for _ in range(12):          # keep the distance label off ligand atoms and bonds
                    if np.linalg.norm(pix - m, axis=1).min() > 34: break
                    m = m + nrm * 8
                placed.append(m)
                lab = f"{dist:.2f} Å" if t != "Water bridge" else f"via {wname}\n{dist:.2f} Å"
                ax.text(m[0], m[1], lab, fontsize=8, color=lc, ha="center", va="center", zorder=4, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.15", fc=BG, ec="none", alpha=0.9))
        face = ST[r["main"]][0]
        ax.add_patch(Circle(r["xy"], R_DISC, fc=face, ec=EDGE, lw=1.2, zorder=3))
        label = n.replace(" (", "\n(")
        ax.text(*r["xy"], label, ha="center", va="center", fontsize=10 if "\n" not in label else 9, fontweight="bold", zorder=5, color=FG if style == "dark" and False else "black")

    title = d["compound"] + (f"  —  {d['target_name']}" if d.get("target_name") else "")
    ax.text((x0 + x1) / 2, y0 + 18, title, ha="center", va="top", fontsize=15, fontweight="bold", color=FG)
    if d.get("affinity") is not None:
        ax.text((x0 + x1) / 2, y0 + 58, f"{d.get('score_label', 'AutoDock Vina binding affinity')}: {d['affinity']:.2f} kcal/mol",
                ha="center", va="top", fontsize=11.5, color=FG)
    present = sorted({r["main"] for r in res.values()} | {l[0] for r in res.values() for l in r["lines"]}, key=lambda t: ST[t][3])
    handles = []
    for t in present:
        face, lc, ls, _ = ST[t]
        kw = dict(color=lc, lw=2.4, ls=ls) if lc else dict(lw=0)
        handles.append(Line2D([0], [0], marker="o", markersize=12, markerfacecolor=face, markeredgecolor=EDGE,
                              label=LEGEND_NAME.get(t, t), **kw))
    lax = fig.add_axes([0, 0, 1, LEG / (bh / 100 + LEG)]); lax.axis("off")
    if handles:
      lax.legend(handles=handles, loc="center", ncol=max(1, min(3 if bw < 1300 else 4, len(handles))), fontsize=10, frameon=True,
                 title="Interactions", title_fontsize=11, handlelength=3.2, columnspacing=1.4,
                 facecolor=th["legend_bg"], edgecolor=EDGE, labelcolor=FG)
      lax.get_legend().get_title().set_color(FG)
    fig.savefig(png, dpi=dpi, facecolor=BG, pil_kwargs={"optimize": True}); plt.close(fig)
    return len(res)


