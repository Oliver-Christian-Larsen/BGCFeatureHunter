import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import pyopenms as poms
import glob
import os

MZML_DIR      = "./data/"
CSV_FILE_PATH = "./output/unique_features.csv"
OUTPUT_PDF    = "waterfall.pdf"

PPM_TOLERANCE = 10
RT_BUFFER_SEC = 120.0

# Values are substrings matched against mzML *filenames* (case-sensitive).
#Plotting name : filesubstring
GENOTYPE_SUBSTRINGS = {
    "KO": "_KO_",
    "WT": "_WT_",
    "OE": "_OE_",
}
#Can be used if more than one medium is present
#plotting name : medium substring in mzml filename
MEDIA_SUBSTRINGS = {
    "PDA": "_P_",
    "YES": "_Y_",
}

# Color family encodes genotype; linestyle encodes media.
#Use the two previously established plotting names
GENOTYPE_COLORMAPS = {
    "WT": cm.Blues,
    "KO": cm.Reds,
    "OE": cm.Purples
}
MEDIA_LINESTYLES = {
    "YES": "solid",
    "PDA": "dashed",
}

# Extra vertical spacing inserted between genotype blocks in the waterfall.
GENOTYPE_GAP = 2.0

def _sci_formatter(x, pos):
    if x == 0:
        return "0"
    if abs(x) >= 1000:
        exp   = int(np.floor(np.log10(abs(x))))
        coeff = x / 10**exp
        return fr"${coeff:.1f}\times10^{{{exp}}}$"
    return f"{x:g}"

SCI_FORMATTER = ticker.FuncFormatter(_sci_formatter)


def classify_files(all_mzml_files):
    classified = {g: {m: [] for m in MEDIA_SUBSTRINGS} for g in GENOTYPE_SUBSTRINGS}
    unmatched  = []

    for fp in all_mzml_files:
        fname    = os.path.basename(fp)
        genotype = next((g for g, sub in GENOTYPE_SUBSTRINGS.items() if sub in fname), None)
        media    = next((m for m, sub in MEDIA_SUBSTRINGS.items()    if sub in fname), None)

        if genotype and media:
            classified[genotype][media].append(fp)
        else:
            unmatched.append(fname)

    if unmatched:
        print(f"[Warning] {len(unmatched)} file(s) could not be classified "
              f"(no matching genotype+media substring): {unmatched}")

    return classified


def build_ordered_samples(classified):
    samples = []

    for genotype in GENOTYPE_SUBSTRINGS:
        cmap = GENOTYPE_COLORMAPS[genotype]

        files_for_genotype = []
        for media in MEDIA_SUBSTRINGS:
            for fp in sorted(classified[genotype][media]):
                files_for_genotype.append((fp, media))

        n = len(files_for_genotype)
        shade_values = np.linspace(0.40, 0.90, max(n, 1))

        for i, (fp, media) in enumerate(files_for_genotype):
            samples.append({
                "filepath":  fp,
                "label":     os.path.basename(fp),
                "genotype":  genotype,
                "media":     media,
                "color":     cmap(shade_values[i]),
                "linestyle": MEDIA_LINESTYLES[media],
            })

    y = 0.0
    current_geno = None
    for s in samples:
        if current_geno is not None and s["genotype"] != current_geno:
            y += GENOTYPE_GAP
        s["y_pos"] = y
        y += 1.0
        current_geno = s["genotype"]

    return samples

class MS1Cache:
    def __init__(self, filepath):
        self.filepath = filepath
        self.filename = os.path.basename(filepath)
        self.rts = []
        self.mz_arrays = []
        self.int_arrays = []
        self.loaded = False

        try:
            exp = poms.OnDiscMSExperiment()
            exp.openFile(filepath)

            for i in range(exp.getNrSpectra()):
                spec = exp.getSpectrum(i)
                if spec.getMSLevel() == 1:
                    self.rts.append(spec.getRT())
                    mz, inten = spec.get_peaks()
                    self.mz_arrays.append(mz)
                    self.int_arrays.append(inten)

            self.rts    = np.array(self.rts)
            self.loaded = True

        except Exception as e:
            print(f"[Error] Could not load {self.filename}: {e}")

    def extract_trace(self, target_mz, ppm, rt_min, rt_max):
        if not self.loaded:
            return np.array([]), np.array([])

        idx_start = np.searchsorted(self.rts, rt_min)
        idx_end   = np.searchsorted(self.rts, rt_max)
        if idx_start >= idx_end:
            return np.array([]), np.array([])

        delta   = target_mz * ppm * 1e-6
        mz_low  = target_mz - delta
        mz_high = target_mz + delta

        xic_rt, xic_int = [], []
        for k in range(idx_start, idx_end):
            mzs = self.mz_arrays[k]
            ints = self.int_arrays[k]
            l = np.searchsorted(mzs, mz_low)
            r = np.searchsorted(mzs, mz_high)
            xic_rt.append(self.rts[k])
            xic_int.append(float(np.sum(ints[l:r])) if l < r else 0.0)

        return np.array(xic_rt), np.array(xic_int)

class WaterfallPipeline:
    def __init__(self):
        if not os.path.exists(MZML_DIR):
            raise FileNotFoundError(f"mzML directory not found: {MZML_DIR}")
        if not os.path.exists(CSV_FILE_PATH):
            raise FileNotFoundError(f"Feature CSV not found: {CSV_FILE_PATH}")

        self.df = pd.read_csv(CSV_FILE_PATH)
        if self.df.empty:
            raise ValueError("Feature CSV is empty.")

        all_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
        if not all_files:
            raise FileNotFoundError(f"No .mzML files found in {MZML_DIR}")

        classified = classify_files(all_files)
        self.samples = build_ordered_samples(classified)

        if not self.samples:
            raise ValueError(
                "No mzML files could be classified. "
                "Check GENOTYPE_SUBSTRINGS and MEDIA_SUBSTRINGS."
            )

        print(f"Classified {len(self.samples)} file(s):")
        for s in self.samples:
            print(f"[genotype={s['genotype']:8s}  media={s['media']:8s}]  {s['label']}")

        self.caches = {}

    def _get_cache(self, filepath):
        if filepath not in self.caches:
            print(f"Loading spectra: {os.path.basename(filepath)}")
            self.caches[filepath] = MS1Cache(filepath)
        return self.caches[filepath]


    def _legend_handles(self):
        handles = []

        handles.append(
            Line2D([0], [0], color='none', label='Genotype')
        )
        for geno, cmap in GENOTYPE_COLORMAPS.items():
            handles.append(
                Line2D([0], [0], color=cmap(0.70), lw=2.5, label=geno)
            )

        handles.append(
            Line2D([0], [0], color='none', label='Media')
        )
        for media, ls in MEDIA_LINESTYLES.items():
            handles.append(
                Line2D([0], [0], color='dimgray', lw=2, linestyle=ls, label=media)
            )

        return handles

    def run(self):
        features = [(float(rt), float(mz)) for rt, mz in self.df[['rt', 'mz']].itertuples(index=False, name=None)]
        total    = len(features)
        print(f"Processing {total} feature(s)")

        with PdfPages(OUTPUT_PDF) as pdf:
            for idx, (target_rt, target_mz) in enumerate(features, 1):
                print(f"[{idx:>{len(str(total))}}/{total}]  "
                      f"mz={target_mz:.4f}  rt={target_rt:.1f} s")

                rt_min = max(0.0, target_rt - RT_BUFFER_SEC)
                rt_max = target_rt + RT_BUFFER_SEC
                common_time = np.linspace(rt_min, rt_max, 500)
                t_min_axis = common_time / 60.0

                fig = plt.figure(figsize=(12, 8))
                ax  = fig.add_subplot(111, projection='3d')

                for s in sorted(self.samples, key=lambda x: x['y_pos'], reverse=True):
                    cache = self._get_cache(s['filepath'])
                    if not cache.loaded:
                        continue

                    rt_raw, int_raw = cache.extract_trace(
                        target_mz, PPM_TOLERANCE, rt_min, rt_max
                    )

                    ax.plot(
                        t_min_axis,
                        np.full_like(t_min_axis, s['y_pos']),
                        int_raw,
                        color     = s['color'],
                        linestyle = s['linestyle'],
                        lw        = 1.8,
                        alpha     = 0.90,
                    )

                ax.view_init(elev=25, azim=-60)
                ax.zaxis.set_major_formatter(SCI_FORMATTER)

                ax.set_xlabel("Retention Time (min)", labelpad=10, fontweight='bold')
                ax.set_zlabel("Intensity", labelpad=15, fontweight='bold')
                ax.set_yticks([])

                ax.set_title(
                    f"m/z {target_mz:.4f}   |   RT {target_rt / 60:.2f} min",
                    fontweight='bold', pad=20,
                )

                ax.legend(
                    handles       = self._legend_handles(),
                    loc           = 'center left',
                    bbox_to_anchor= (1.05, 0.5),
                    frameon       = False,
                    title         = "Condition",
                    title_fontsize= 10,
                )

                plt.tight_layout()
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

        print(f"Output saved to {OUTPUT_PDF}")


if __name__ == "__main__":
    pipeline = WaterfallPipeline()
    pipeline.run()
