## Fit a model of MSL, GBR, and Vindija populations.

import moments 
import sys


def main():
    # define paths and parameters
    tag = sys.argv[1]
    graph_file = "MSL_GBR_Vindija_model3.yaml"
    options_file = "options_MSL_GBR_Vindija3.yaml"
    u = 1.5e-8
    L = 960914001
    data = moments.Spectrum.from_file("../../data/data_MSL_GBR_Vindija")
    output = f"MSL_GBR_Vindija_round3_misid_best_fit_{tag}.yaml"
    param_names, misid_popt, misid_ll = moments.Demes.Inference.optimize(
        graph_file, 
        options_file, 
        data, 
        maxiter=10000,
        fit_ancestral_misid=True, 
        misid_guess=0.03,
        uL=u * L, 
        verbose=1,
        output=output, 
        perturb=0.25,
        method="powell"
    )
    print("Best-fit log-likelihood:", -misid_ll)
    print("Best-fit parameters:")
    for name, popt in zip(param_names, misid_popt):
        print(f"{name}\t{popt:.3}")

    return 


main() 
