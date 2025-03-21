## Fit a model of MSL, GBR, and Vindija populations.

import moments 


def main():
    # define paths and parameters
    graph_file = "MSL_GBR_Vindija_model.yaml"
    options_file = "options_MSL_GBR_Vindija.yaml"
    data_file = "../../spectra/MSL_GBR_Vindija"
    output = "MSL_GBR_Vindija_model.misid_fit.yaml"
    u = 1.5e-8
    L = 960914001
    misid_guess = 0.03

    data = moments.Spectrum.from_file(data_file)

    # fit using the LBFGSB algorithm
    param_names, fit_params, ll = moments.Demes.Inference.optimize(
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
    print("Best-fit log-likelihood:", -ll)
    print("Best-fit parameters:")
    for name, val in zip(param_names, fit_params):
        print(f"{name}\t{val:.3}")

    return 


main() 
