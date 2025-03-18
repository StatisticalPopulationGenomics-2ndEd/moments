## Compute uncertainties using the FIM method

import demes
import moments 


def main():
    # define paths and parameters
    graph_file = "MSL_GBR_Vindija_round2.misid_fit.yaml"
    options_file = "options_MSL_GBR_Vindija_round2.yaml"
    data_file = "../../spectra/MSL_GBR_Vindija"
    u = 1.5e-8
    L = 960914001
    p_misid = 0.0223

    data = moments.Spectrum.from_file(data_file)

    # load parameter names and best-fit values using moments.Demes functions
    builder = demes.load(graph_file).asdict()
    options = moments.Demes.Inference._get_params_dict(options_file)
    param_names, fit_params, _, __ = \
        moments.Demes.Inference._set_up_params_and_bounds(options, builder)

    std_errs = moments.Demes.Inference.uncerts(
        graph_file,
        options_file,
        data,
        uL=u * L,
        fit_ancestral_misid=True,
        misid_fit=p_misid,
    )

    # print results
    print(r"95% confidence intervals:")
    print("param\t2.5%\t97.5%")
    for name, val, err in zip(param_names, fit_params, std_errs):
        print(f"{name}\t{val - 1.96 * err:.3}\t{val + 1.96 * err:.3}")

    return 


main() 
