## Refit the MSL, GBR, Vindija model to fine-tune parameters

import moments 


def main():
    # define paths and parameters
    graph_file = "MSL_GBR_Vindija_round2.yaml"
    options_file = "options_MSL_GBR_Vindija_round2.yaml"
    data_file = "../../spectra/MSL_GBR_Vindija"
    output = "MSL_GBR_Vindija_round2.misid_fit.yaml"
    u = 1.5e-8
    L = 960914001
    misid_guess = 0.024

    data = moments.Spectrum.from_file(data_file)

    param_names, misid_popt, misid_ll = moments.Demes.Inference.optimize(
        graph_file, 
        options_file, 
        data, 
        maxiter=10000,
        fit_ancestral_misid=True, 
        misid_guess=misid_guess,
        uL=u * L, 
        verbose=1,
        output=output,
        method="lbfgsb"
    )

    # print results
    print("Best-fit log-likelihood:", -misid_ll)
    print("Best-fit parameters:")
    for name, popt in zip(param_names, misid_popt):
        print(f"{name}\t{popt:.3}")

    return 


main() 
