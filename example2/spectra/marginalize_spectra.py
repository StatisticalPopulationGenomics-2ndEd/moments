
import moments 


full_fs = moments.Spectrum.from_file("MSL_GBR_Vindija.fs")

fs_MSL = full_fs.marginalize([1, 2])
fs_MSL.to_file("MSL.fs")

fs_GBR = full_fs.marginalize([0, 2])
fs_GBR.to_file("GBR.fs")

fs_MSL_GBR = full_fs.marginalize([2])
fs_MSL_GBR.to_file("MSL_GBR.fs")
