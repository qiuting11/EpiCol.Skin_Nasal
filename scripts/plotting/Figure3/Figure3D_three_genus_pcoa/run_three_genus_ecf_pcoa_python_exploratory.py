from pathlib import Path
from html import escape
import sys

import numpy as np
import pandas as pd

DEFAULT_INPUT = Path("data/external/figure3D_three_genus_pcoa/mag_cf_long_with_core_phenotypes.tsv")
INPUT = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_INPUT
OUT_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("outputs/Figure3/Figure3D_three_genus_pcoa_python_exploratory")

TARGET_GENERA = {
    "Staphylococcus": "Staphylococcus",
    "Corynebacterium": "Corynebacterium",
    "Cutibacterium_or_Propionibacterium": "Cutibacterium",
}

COLORS = {
    "Staphylococcus": "#3B7EA1",
    "Corynebacterium": "#D97742",
    "Cutibacterium": "#4E9F6D",
}


def jaccard_distance_binary(x):
    x = x.astype(bool).astype(float)
    intersection = x @ x.T
    row_sum = x.sum(axis=1)
    union = row_sum[:, None] + row_sum[None, :] - intersection
    similarity = np.ones_like(union, dtype=float)
    np.divide(intersection, union, out=similarity, where=union != 0)
    distance = 1.0 - similarity
    np.fill_diagonal(distance, 0.0)
    return distance


def braycurtis_distance(x):
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    distance = np.zeros((n, n), dtype=float)
    row_sum = x.sum(axis=1)
    for i in range(n):
        numerator = np.abs(x[i + 1 :] - x[i]).sum(axis=1)
        denominator = row_sum[i + 1 :] + row_sum[i]
        values = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator != 0)
        distance[i, i + 1 :] = values
        distance[i + 1 :, i] = values
    return distance


def classical_pcoa(distance_matrix):
    n = distance_matrix.shape[0]
    d2 = distance_matrix**2
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ d2 @ centering
    eigvals, eigvecs = np.linalg.eigh(gram)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > 1e-12
    coords = eigvecs[:, positive] * np.sqrt(eigvals[positive])
    explained = eigvals[positive] / eigvals[positive].sum()
    return coords, explained


def permanova(distance_matrix, groups, n_perm=999, seed=20260827):
    rng = np.random.default_rng(seed)
    groups = np.asarray(groups)
    n = len(groups)
    k = len(np.unique(groups))
    d2 = distance_matrix**2

    def pseudo_f(labels):
        total_ss = d2.sum() / n
        within_ss = 0.0
        for group in np.unique(labels):
            idx = np.where(labels == group)[0]
            within_ss += d2[np.ix_(idx, idx)].sum() / len(idx)
        between_ss = total_ss - within_ss
        f_stat = (between_ss / (k - 1)) / (within_ss / (n - k))
        r2 = between_ss / total_ss
        return f_stat, r2

    observed_f, observed_r2 = pseudo_f(groups)
    perm_f = np.empty(n_perm)
    for i in range(n_perm):
        perm_f[i], _ = pseudo_f(rng.permutation(groups))
    p_value = (np.sum(perm_f >= observed_f) + 1) / (n_perm + 1)
    return observed_f, observed_r2, p_value


def anova_f(values, groups):
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    overall = values.mean()
    unique_groups = np.unique(groups)
    ss_between = sum(
        np.sum(groups == group) * (values[groups == group].mean() - overall) ** 2
        for group in unique_groups
    )
    ss_within = sum(
        ((values[groups == group] - values[groups == group].mean()) ** 2).sum()
        for group in unique_groups
    )
    return (ss_between / (len(unique_groups) - 1)) / (ss_within / (len(values) - len(unique_groups)))


def betadisper(coords, groups, n_perm=999, seed=20260828):
    rng = np.random.default_rng(seed)
    groups = np.asarray(groups)

    def distances_to_centroid(labels):
        distances = np.zeros(len(labels))
        for group in np.unique(labels):
            idx = np.where(labels == group)[0]
            centroid = coords[idx].mean(axis=0)
            distances[idx] = np.linalg.norm(coords[idx] - centroid, axis=1)
        return distances

    observed_distances = distances_to_centroid(groups)
    observed_f = anova_f(observed_distances, groups)
    perm_f = np.empty(n_perm)
    for i in range(n_perm):
        labels = rng.permutation(groups)
        distances = distances_to_centroid(labels)
        perm_f[i] = anova_f(distances, labels)
    p_value = (np.sum(perm_f >= observed_f) + 1) / (n_perm + 1)
    return observed_f, p_value, observed_distances


def plot_svg(meta, coords, explained, title, stats_text, output_path):
    width, height = 720, 620
    left, right, top, bottom = 82, 36, 62, 82
    plot_w, plot_h = width - left - right, height - top - bottom
    x, y = coords[:, 0], coords[:, 1]
    xpad = (x.max() - x.min()) * 0.08 or 1
    ypad = (y.max() - y.min()) * 0.08 or 1
    xmin, xmax = x.min() - xpad, x.max() + xpad
    ymin, ymax = y.min() - ypad, y.max() + ypad

    def sx(value):
        return left + (value - xmin) / (xmax - xmin) * plot_w

    def sy(value):
        return top + (ymax - value) / (ymax - ymin) * plot_h

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{left}" y="32" font-family="Arial" font-size="18" font-weight="700">{escape(title)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="white" stroke="#222" stroke-width="1"/>',
    ]
    if ymin <= 0 <= ymax:
        zy = sy(0)
        parts.append(f'<line x1="{left}" x2="{left + plot_w}" y1="{zy:.2f}" y2="{zy:.2f}" stroke="#c7c7c7" stroke-width="1"/>')
    if xmin <= 0 <= xmax:
        zx = sx(0)
        parts.append(f'<line x1="{zx:.2f}" x2="{zx:.2f}" y1="{top}" y2="{top + plot_h}" stroke="#c7c7c7" stroke-width="1"/>')

    for genus, sub in meta.groupby("genus_display"):
        for i in sub.index.to_numpy():
            parts.append(
                f'<circle cx="{sx(coords[i, 0]):.2f}" cy="{sy(coords[i, 1]):.2f}" r="3.2" '
                f'fill="{COLORS[genus]}" fill-opacity="0.55" stroke="none"/>'
            )

    parts.append(f'<text x="{left + plot_w / 2}" y="{height - 28}" text-anchor="middle" font-family="Arial" font-size="15">PCoA1 ({explained[0] * 100:.1f}%)</text>')
    parts.append(f'<text x="24" y="{top + plot_h / 2}" text-anchor="middle" font-family="Arial" font-size="15" transform="rotate(-90 24 {top + plot_h / 2})">PCoA2 ({explained[1] * 100:.1f}%)</text>')

    legend_x, legend_y = left + plot_w - 190, top + 18
    for j, (genus, sub) in enumerate(meta.groupby("genus_display")):
        y0 = legend_y + j * 22
        parts.append(f'<circle cx="{legend_x}" cy="{y0}" r="5" fill="{COLORS[genus]}" fill-opacity="0.75"/>')
        parts.append(f'<text x="{legend_x + 12}" y="{y0 + 4}" font-family="Arial" font-size="12">{escape(genus)} (n={len(sub)})</text>')

    box_x, box_y = left + 14, top + plot_h - 58
    parts.append(f'<rect x="{box_x}" y="{box_y}" width="295" height="44" rx="4" fill="white" fill-opacity="0.9" stroke="#dddddd"/>')
    for j, line in enumerate(stats_text.split("\n")):
        parts.append(f'<text x="{box_x + 10}" y="{box_y + 18 + j * 17}" font-family="Arial" font-size="12">{escape(line)}</text>')
    parts.append("</svg>")
    output_path.write_text("\n".join(parts), encoding="utf-8")


def run_analysis(matrix, meta, metric, label, prefix):
    if metric == "jaccard":
        distance = jaccard_distance_binary(matrix.values)
    elif metric == "braycurtis":
        distance = braycurtis_distance(np.log1p(matrix.values))
    else:
        raise ValueError(f"Unknown metric: {metric}")

    coords, explained = classical_pcoa(distance)
    permanova_f, permanova_r2, permanova_p = permanova(distance, meta["genus_display"].values)
    betadisper_f, betadisper_p, centroid_distances = betadisper(coords, meta["genus_display"].values)

    coord_df = meta.copy()
    coord_df["PCoA1"] = coords[:, 0]
    coord_df["PCoA2"] = coords[:, 1]
    coord_df["distance_to_genus_centroid"] = centroid_distances
    coord_df.to_csv(OUT_DIR / f"{prefix}_coordinates.tsv", sep="\t", index=False)

    stats_text = (
        f"PERMANOVA R2={permanova_r2:.3f}, P={permanova_p:.3g}\n"
        f"betadisper P={betadisper_p:.3g}; 999 permutations"
    )
    plot_svg(meta, coords, explained, label, stats_text, OUT_DIR / f"{prefix}.svg")

    return {
        "analysis": label,
        "metric": metric,
        "n_mags": len(meta),
        "n_features": matrix.shape[1],
        "pcoa1_explained": explained[0],
        "pcoa2_explained": explained[1],
        "permanova_F": permanova_f,
        "permanova_R2": permanova_r2,
        "permanova_p": permanova_p,
        "permutations": 999,
        "betadisper_F": betadisper_f,
        "betadisper_p": betadisper_p,
    }


def maybe_write_pngs():
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return

    stats = pd.read_csv(OUT_DIR / "three_genus_pcoa_statistics.tsv", sep="\t")
    coord_files = [
        ("three_genus_ecf_presence_jaccard_pcoa_coordinates.tsv", "A. eCF presence/absence profiles"),
        ("three_genus_ecf_copy_braycurtis_pcoa_coordinates.tsv", "B. eCF copy-number profiles"),
    ]
    try:
        font = ImageFont.truetype("arial.ttf", 18)
        font_small = ImageFont.truetype("arial.ttf", 13)
        font_title = ImageFont.truetype("arialbd.ttf", 20)
    except Exception:
        font = font_small = font_title = ImageFont.load_default()

    panels = []
    width, height = 760, 620
    for coord_file, title in coord_files:
        df = pd.read_csv(OUT_DIR / coord_file, sep="\t")
        row = stats[stats["analysis"].eq(title)].iloc[0]
        img = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(img, "RGBA")
        left, right, top, bottom = 86, 40, 64, 82
        plot_w, plot_h = width - left - right, height - top - bottom
        x, y = df["PCoA1"], df["PCoA2"]
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        xpad = (xmax - xmin) * 0.08 or 1
        ypad = (ymax - ymin) * 0.08 or 1
        xmin, xmax = xmin - xpad, xmax + xpad
        ymin, ymax = ymin - ypad, ymax + ypad

        def sx(value):
            return left + (value - xmin) / (xmax - xmin) * plot_w

        def sy(value):
            return top + (ymax - value) / (ymax - ymin) * plot_h

        draw.text((left, 26), title, fill=(0, 0, 0), font=font_title)
        draw.rectangle((left, top, left + plot_w, top + plot_h), outline=(35, 35, 35), width=1)
        if ymin <= 0 <= ymax:
            zy = sy(0)
            draw.line((left, zy, left + plot_w, zy), fill=(190, 190, 190), width=1)
        if xmin <= 0 <= xmax:
            zx = sx(0)
            draw.line((zx, top, zx, top + plot_h), fill=(190, 190, 190), width=1)
        for genus, sub in df.groupby("genus_display"):
            color = tuple(int(COLORS[genus].lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
            for _, point in sub.iterrows():
                cx, cy = sx(point["PCoA1"]), sy(point["PCoA2"])
                draw.ellipse((cx - 3, cy - 3, cx + 3, cy + 3), fill=color + (135,))
        draw.text((left + plot_w / 2 - 95, height - 42), f"PCoA1 ({row['pcoa1_explained'] * 100:.1f}%)", fill=(0, 0, 0), font=font)
        draw.text((16, top + plot_h / 2 - 10), f"PCoA2 ({row['pcoa2_explained'] * 100:.1f}%)", fill=(0, 0, 0), font=font)
        lx, ly = left + plot_w - 205, top + 15
        for i, (genus, sub) in enumerate(df.groupby("genus_display")):
            yy = ly + i * 24
            color = tuple(int(COLORS[genus].lstrip("#")[j : j + 2], 16) for j in (0, 2, 4))
            draw.ellipse((lx - 5, yy - 5, lx + 5, yy + 5), fill=color + (190,))
            draw.text((lx + 14, yy - 8), f"{genus} (n={len(sub)})", fill=(0, 0, 0), font=font_small)
        bx, by = left + 14, top + plot_h - 60
        draw.rounded_rectangle((bx, by, bx + 320, by + 48), radius=5, fill=(255, 255, 255, 230), outline=(215, 215, 215))
        draw.text((bx + 10, by + 8), f"PERMANOVA R2={row['permanova_R2']:.3f}, P={row['permanova_p']:.3g}", fill=(0, 0, 0), font=font_small)
        draw.text((bx + 10, by + 27), f"betadisper P={row['betadisper_p']:.3g}; 999 permutations", fill=(0, 0, 0), font=font_small)
        panels.append(img)

    combined = Image.new("RGB", (width * 2, height), "white")
    combined.paste(panels[0], (0, 0))
    combined.paste(panels[1], (width, 0))
    panels[0].save(OUT_DIR / "three_genus_ecf_presence_jaccard_pcoa.png", dpi=(300, 300))
    panels[1].save(OUT_DIR / "three_genus_ecf_copy_braycurtis_pcoa.png", dpi=(300, 300))
    combined.save(OUT_DIR / "three_genus_ecf_pcoa_combined.png", dpi=(300, 300))


def main():
    if not INPUT.exists():
        raise FileNotFoundError(f"Missing input file: {INPUT}. See docs/data_availability.md for the expected external data layout.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    long = pd.read_csv(INPUT, sep="\\t")
    long = long[long["target_genus_group"].isin(TARGET_GENERA)].copy()
    long["genus_display"] = long["target_genus_group"].map(TARGET_GENERA)
    long["presence_binary"] = (long["copy_number"] > 0).astype(int)

    presence = long.pivot_table(index="MAG", columns="CF_id", values="presence_binary", aggfunc="max", fill_value=0)
    copy_number = long.pivot_table(index="MAG", columns="CF_id", values="copy_number", aggfunc="sum", fill_value=0)

    meta_cols = ["MAG", "genus_display", "target_genus_group", "Species"]
    if "source_group" in long.columns:
        meta_cols.append("source_group")
    meta = long[meta_cols].drop_duplicates("MAG").set_index("MAG").loc[presence.index].reset_index()

    stats = [
        run_analysis(presence, meta, "jaccard", "A. eCF presence/absence profiles", "three_genus_ecf_presence_jaccard_pcoa"),
        run_analysis(copy_number, meta, "braycurtis", "B. eCF copy-number profiles", "three_genus_ecf_copy_braycurtis_pcoa"),
    ]
    pd.DataFrame(stats).to_csv(OUT_DIR / "three_genus_pcoa_statistics.tsv", sep="\t", index=False)
    meta.groupby("genus_display").size().reset_index(name="n_mags").to_csv(OUT_DIR / "three_genus_pcoa_sample_counts.tsv", sep="\t", index=False)
    maybe_write_pngs()
    print(f"Wrote outputs to {OUT_DIR.resolve()}")
    print(pd.DataFrame(stats).to_string(index=False))


if __name__ == "__main__":
    main()


