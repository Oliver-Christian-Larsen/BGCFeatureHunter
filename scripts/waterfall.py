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

def _sci_formatter(x, pos):
    if x == 0:
        return "0"
    if abs(x) >= 1000:
        exp   = int(np.floor(np.log10(abs(x))))
        coeff = x / 10**exp
        return fr"${coeff:.1f}\times10^{{{exp}}}$"
    return f"{x:g}"

SCI_FORMATTER = ticker.FuncFormatter(_sci_formatter)

GENOTYPE_COLORMAPS = {"WT": cm.Blues, "KO": cm.Reds, "OE": cm.Purples}
LINESTYLE_CYCLE    = ["solid", "dashed", "dotted", "dashdot"]
NO_MEDIA = "default"

def build_ordered_samples(classified_files, requested_media):
    plottable = [
        f for f in classified_files
        if f.get("group") not in (None, "blank")
    ]

    if requested_media:
        found = {f["media"] for f in plottable if f["media"] is not None}

        if not found:
            raise ValueError(
                f"None of the requested media {requested_media} were found "
                "in any filename. Check the --media substrings."
            )

        for m in requested_media:
            if m not in found:
                print(f"WARNING: medium '{m}' was requested but matched no files.")

        media_types = [m for m in requested_media if m in found]

        for f in plottable:
            if f["media"] is None:
                print(f"WARNING: {f['filename']} matched no requested medium and will not be plotted.")
        plottable = [f for f in plottable if f["media"] is not None]
    else:
        media_types = [NO_MEDIA]

    linestyles = {
        medium: LINESTYLE_CYCLE[i % len(LINESTYLE_CYCLE)]
        for i, medium in enumerate(media_types)
    }
    media_indices = {medium: i for i, medium in enumerate(media_types)}

    shade_values = np.linspace(0.40, 0.90, max(len(plottable), 1))

    samples = []
    for i, file_info in enumerate(plottable):
        gtype = file_info["group"]
        medium = file_info["media"] if requested_media else NO_MEDIA

        samples.append({
            **file_info,
            "genotype": gtype,
            "media": medium,
            "media_index": media_indices[medium],
            "color": GENOTYPE_COLORMAPS[gtype](shade_values[i]),
            "linestyle": linestyles[medium],
        })

    # Vertical positions, with a gap between genotype blocks
    y = 0.0
    current_genotype = None
    for sample in samples:
        if current_genotype is not None and sample["genotype"] != current_genotype:
            y += 2
        sample["y_pos"] = y
        y += 1.0
        current_genotype = sample["genotype"]

    return samples



class MS1Cache:
    def __init__(self, filepath):
    
        self.filepath  = filepath
        self.filename  = os.path.basename(filepath)
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
    def __init__(self,arg_dic,classified_files):
        self.output_pdf = os.path.join(arg_dic["paths"]["out"], "waterfall.pdf")
        if not os.path.exists(f"{arg_dic['paths']['out']}/unique_features.csv"):
            raise FileNotFoundError(f"Unique Feature CSV not found: {f"{arg_dic['paths']['out']}/unique_features.csv"}")

        self.df = pd.read_csv(f"{arg_dic['paths']['out']}/unique_features.csv")

        if self.df.empty:
            raise ValueError("Feature CSV is empty.")

        all_files = sorted(glob.glob(os.path.join(arg_dic['paths']['mzml'], "*.mzML")))
        if not all_files:
            raise FileNotFoundError(f"No .mzML files found in {arg_dic['paths']['mzml']}")


        self.samples = build_ordered_samples(classified_files, arg_dic["media"])
        self.genotypes = list(dict.fromkeys(s["genotype"] for s in self.samples))
        self.linestyles = {s["media"]: s["linestyle"] for s in self.samples}
        
        if not self.samples:
            raise ValueError(
                "No mzML files could be classified. "
                "Check GENOTYPE_SUBSTRINGS and MEDIA_SUBSTRINGS."
            )

        self.caches = {}

    def _get_cache(self, filepath):
        if filepath not in self.caches:
            self.caches[filepath] = MS1Cache(filepath)
        return self.caches[filepath]


    def _legend_handles(self,arg_dic):
        handles = []

        handles.append(
            Line2D([0], [0], color='none', label='Genotype')
        )

        for geno in self.genotypes:
            show_name = arg_dic["display_names"][geno]
            handles.append(Line2D([0], [0], color=GENOTYPE_COLORMAPS[geno](0.70), lw=2.5, label=show_name))


        handles.append(
            Line2D([0], [0], color='none', label='Media')
        )

        for media, ls in self.linestyles.items():
            handles.append(Line2D([0], [0], color='dimgray', lw=2, linestyle=ls, label=media))
            

        return handles

    def run(self,arg_dic,config):
        print("Hold on, the waterfalls are coming!")
        features = [(float(rt), float(mz)) for rt, mz in self.df[['rt', 'mz']].itertuples(index=False, name=None)]
        total    = len(features)
        print(f"Processing {total} feature(s)")

        with PdfPages(self.output_pdf) as pdf:
            for idx, (target_rt, target_mz) in enumerate(features, 1):
                if idx % 20 == 0:
                    print(f"[{idx:>{len(str(total))}}/{total}]")

                rt_min = max(0.0, target_rt - config["visual_inspection"]["RT_window"])
                rt_max = target_rt + config["visual_inspection"]["RT_window"]
                

                fig = plt.figure(figsize=(12, 8))
                ax  = fig.add_subplot(111, projection='3d')

                for s in sorted(self.samples, key=lambda x: x['y_pos'], reverse=True):
                    cache = self._get_cache(s['filepath'])
                    if not cache.loaded:
                        continue

                    rt_raw, int_raw = cache.extract_trace(
                        target_mz, config["visual_inspection"]["ppm_error"], rt_min, rt_max
                    )

                    t_min_axis = rt_raw / 60.0

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
                    handles       = self._legend_handles(arg_dic),
                    loc           = 'center left',
                    bbox_to_anchor= (1.05, 0.5),
                    frameon       = False,
                    title         = "Condition",
                    title_fontsize= 10,
                )

                plt.tight_layout()
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

        print(f"Output saved to {self.output_pdf}")


if __name__ == "__main__":
    pipeline = WaterfallPipeline(arg_dic)
    pipeline.run(arg_dic)
