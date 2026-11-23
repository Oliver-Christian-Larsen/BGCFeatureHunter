from find_mfs import IsotopeMatchConfig, FormulaFinder, FormulaPrior
import numpy as np
import os

def make_np_array(entry):
    mz_m0 = entry["mz_m0"]
    mz_m1 = entry["mz_m1"]
    mz_m2 = entry["mz_m2"]
    mz_m3 = entry["mz_m3"]

    use_envelope = False

    if mz_m3 > 0.0:
        use_envelope = True
        array = np.array([
            [mz_m0,entry["ratio_m0"]],
            [mz_m1,entry["ratio_m1"]],
            [mz_m2,entry["ratio_m2"]],
            [mz_m3,entry["ratio_m3"]]]
            )
        return array, use_envelope, mz_m0

    elif mz_m2 > 0.0:
        use_envelope = True
        array = np.array([
            [mz_m0,entry["ratio_m0"]],
            [mz_m1,entry["ratio_m1"]],
            [mz_m2,entry["ratio_m2"]]],
            )
        return array, use_envelope, mz_m0

    elif mz_m1 > 0.0:
        use_envelope = True
        array = np.array([
            [mz_m0,entry["ratio_m0"]],
            [mz_m1,entry["ratio_m1"]]],
            )
        return array, use_envelope, mz_m0

    else:
        array = np.array([
            [mz_m0,entry["ratio_m0"]]],
            )
        return array, use_envelope, mz_m0

def m0_only(ffinder,charge,mz_m0): 
    results = ffinder.find_formulae(
        mass=mz_m0,
        charge=charge,
        error_ppm=10,
        check_octet=True,
        filter_rdbe=(0,20),
        max_counts={
        'P': 0,
        'S': 2,
    }, 
    )


def envelope(ffinder,charge,mz_m0,array):
    iso_config = IsotopeMatchConfig(
        envelope=array,
        mz_tolerance_da=0.1,
        minimum_rmse=0.05,
    )
    results = ffinder.find_formulae(
    mass=mz_m0,
    charge=charge,
    error_ppm=10.0,
    check_octet=True,
    filter_rdbe=(0, 20),
    max_counts={
        'P': 0,
        'S': 2,
    },
    isotope_match=iso_config,
    )
    prior = FormulaPrior.default()
    prior.score_results(
      results,
      mass_sigma_ppm=10.0,
      isotope_sigma=0.1,
    )
    ranked = results.sort_by_prior()
    return ranked



def save_data(ranked,cid,mf_table):
    with open(mf_table, "a") as fh:
        fh.write(f"{cid}\n{ranked.to_table()}\n\n")

def main(data,arg_dic):
    ffinder = FormulaFinder()
    mf_table = os.path.join(arg_dic["paths"]["out"],"molecular_formulas.txt")

    with open(mf_table, "w"):
        pass

    for entry in data:
        cid = str(entry["consensus_id"])
        charge = entry["charge"]

        if round(charge,0) == 0:
            charge = 1

        array, use_envelope, mz_m0 = make_np_array(entry)
        if use_envelope == False and arg_dic["elemental_single"] == True:
            m0_only(ffinder,charge,mz_m0)

        if use_envelope == True:

            ranked = envelope(ffinder,charge,mz_m0,array)


            save_data(ranked,cid,mf_table)

if __name__ == "__main__":
    main(data,arg_dic)

    

