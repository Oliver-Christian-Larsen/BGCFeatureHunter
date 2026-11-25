import pyopenms as oms
from collections import defaultdict

def make_merger():
    merger = oms.SpectraMerger()
    param = merger.getParameters()
    param.setValue("precursor_method:mz_tolerance", 0.1)  
    param.setValue("precursor_method:rt_tolerance", 10.0)
    merger.setParameters(param)

    return merger

def get_ce(scan):
    p = scan.getPrecursors()[0]
    ce = p.getMetaValue("collision energy")
    if ce and ce > 0:
        return round(float(ce), 1)
    return None

def get_and_merge_spectrum(merger, mzml):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(mzml, exp)

    by_ce = defaultdict(list)
    for s in exp:
        if s.getMSLevel() == 2 and s.getPrecursors():
            by_ce[get_ce(s)].append(s)

    merged = []
    for ce, scans in by_ce.items():
        sub = oms.MSExperiment()
        for s in scans:
            sub.addSpectrum(s)
        sub.sortSpectra(True)

        n_before = sub.size()
        merger.mergeSpectraPrecursors(sub)
        out = [s for s in sub if s.getMSLevel() == 2]

        print(f"CE {ce}: {n_before} scans -> {len(out)} merged")
        merged.extend((ce, s) for s in out)

    return merged




def main():
    mzml = "/home/user/BGCFeatureHunter/data/KW_TES_OCL_ox_01_ms2p_111_01_18621_calibrated_centroided.mzML"

    merger = make_merger()

    ms2_blocks = get_and_merge_spectrum(merger,mzml)

main()