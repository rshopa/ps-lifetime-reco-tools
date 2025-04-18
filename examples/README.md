# Examples

Here, a brief description how to execute the exemplary fitting is given. Please refer to the description of the input parameters (in JSON format) in ```/params_json``` directory.

Four executives are available. The basic usage for the single-histogram minimisation is:

```Rscript [--vanilla] fitOneHistogramLMA2.R -p <json_params> -s <lifetime_spectrum_ascii> -o <output_prefix> [-v] [> <output_log>]```

```Rscript [--vanilla] fitOneHistogramLMA3.R -p <json_params> <lifetime_spectrum_ascii> -o <output_prefix> [-fixed-ints-ratio] [-three-stage] [-d <material_density>] [-v] [> <output_log>]```

And for the multi-histogram (multi-voxel) fitting:

```Rscript [--vanilla] runMultiVoxelFitLMA2.R -p <json_params> -dt <min,max,n_bins> -s <lifetime_spectra> -vox-ids <voxel_coordinates> -o <output_prefix> [-vox-size <x_mm,y_mm,z_mm>] [-preserve-na] [-na-to-zero] [-v] [> <output_log>]```

```Rscript [--vanilla] runMultiVoxelFitLMA3.R -p <json_params> -dt <min,max,n_bins> -s <lifetime_spectra> -vox-ids <voxel_coordinates> -o <output_prefix> [-vox-size <x_mm,y_mm,z_mm>] [-fixed-ints-ratio] [-preserve-na] [-na-to-zero] [-three-stage] [-v] [> <output_log>]```

The ```--vanilla``` argument prevents some user-specific R settings to activate and makes sure no saving or restoring of workspaces is made.

The arguments for each executable can be checked via help using keys ```-h```, ```--help```, ```-?``` or without any argument given, for instance: 

```Rscript [--vanilla] fitOneHistogramLMA2.R -?```

The examples put in this directory are for both separate histograms (in ```/single_histograms```, ASCII format) and multi-voxel data -- for a [NEMA IEC phantom](https://www.spect.com/our-products/nema-iec-pet-body-phantom) (in ```\multi_voxel```, in both binary formats and RDS -- see main README). For the latter, only a fragment of data is given -- for the voxels around the largest phantom sphere.

Note that all scripts require a compilation of ```/cpp/kdeCPP.cpp``` functions each time. A proper R package creation is an option, though not implemented in this repository.
