Tomas van Lieshout
Bachelor student - Natuur- en Sterrenkunde, Radboud Universiteit. 

This repository includes the code I used for my bachelor thesis 'Fast X-ray Transients, a search for extragalactic super fast X-ray transients in the Chandra Source Catalog". 

During the thesis I conducted a search in the CSC 2.1 looking for super fast X-ray transient candidates with the usage of Jonathan Quirola Vázques his search algorithm. Where he used 20 ks search windows, I used a 0.2 ks search window to look for faster transients. 



search
    This pipeline automates the end-to-end processing of Chandra ObsIDs for my thesis.
    It downloads raw data, applies CIAO preprocessing, runs wavdetect, executes my custom transient-search algorithm, and outputs detection catalogs     across multiple time windows.
    The code is modular and can be adapted for other Chandra data mining tasks.
    
    Global configuration:
    DEFAULT_DATA_PATH: where downloaded ObsIDs are stored
    WINDOW_SIZES: list of search window durations in ks
    FILENAMES: input ObsID list files
    VERBOSE: 0 = shows minimal output
             1 = shows informative progress messages
             2 = shows detailed debugging information

filter
    This notebook handles catalog cross-matching for the detected transient candidates.
    It loads filtered detection lists, ensures catalog columns exist, merges new detections, and evaluates them using user-supplied catalog              functions (Gaia, NED, SIMBAD, etc.).
    It also aggregates reference matches across all window sizes and identifies strict and relaxed FXRT candidates.

    Global configuration:
    WINDOW: specific search window in ks for filtering 
    VERBOSE: 0 = only essential output
             1 = shows informative progress messages
             2 = shows catalog matching actions
             3 = shows detailed debugging information 

lightcurves
    This notebook performs PSF-adaptive photometry for all FXRT candidates using Chandra ACIS event files. It estimates the PSF R90 radius from off      axis angle, selects source and background regions dynamically, performs adaptive time-binning, computes global T90, and produces a complete          suite of diagnostic plots (light curves, energy-time diagrams, spectra, and region maps). Outputs are saved in structured folders for scientific     analysis and publication.

    Global configuration:
    data_root: where downloaded ObsIDs are populated from
    candidates_file: where final candidates are listed

counterparts_check
    This notebook checks all FXRT candidates against major optical and infrared catalogs. It loads the cleaned FXRT detection list and queries Pan-      STARRS DR2, GAIA DR3, WISE/IRSA, and relevant VIZIER catalogs within a typical 2″ radius. For each candidate, the script records all possible        positional matches, including photometric properties and offsets. The output allows the identification of isolated FXRTs versus those associated     with known stars or galaxies.

    Global configuration:
    fxrt: which ObsIDs are checked against the databases

variability_timescale
    This notebook implements Bayesian Blocks analysis to measure minimum variability timescales for FXRT candidates. Using unbinned photon arrival       times from Chandra event files, the function mvt_bayesian_blocks() computes the Bayesian Blocks segmentation, extracts block widths, and reports     the shortest significant variability time. This provides a quantitative measure of the rapidity of each transient’s temporal evolution.

    Global configuration:
    obsid: Chandra ObsID for which the MVT is computed. 

time_and_location
    This notebook performs detailed time-domain and positional analysis for FXRT candidates. It loads photon arrival times, computes variability         timescales using Bayesian Blocks, and evaluates the positions of Solar System bodies at the event time to rule out planetary or asteroid             contamination. The notebook also transforms coordinates between geocentric and topocentric frames and generates diagnostics for temporal and         spatial photon distributions.

    Global configuration:
    data_root: where downloaded ObsIDs are populated from
    candidate_file: which ObsIDs are checked
    evtfile: observation time and date of ObsID
    date: date and time of ObsIDs observation
    ra_input: RA of candidate (2000/ICRS coordinates, use dot decimal)
    dec_input: DEC of candidate (2000/ICRS coordinates, use dot decimal)

proper_motion_check
    This notebook loads multi-epoch RA/DEC measurements for potential FXRT counterparts and performs a weighted linear fit to measure proper motion.     It computes uncertainties and covariance, visualizes the proper-motion error ellipse, and evaluates whether the motion is statistically              significant. This determines whether the counterpart is a foreground star (moving) or a stationary extragalactic source.

    Global configuration:
    df: input file with RA and DEC of potential counterpart

OAA_check
    This module converts PSF R90 radii into off-axis angles and computes the 95% positional error radius for Chandra/ACIS detections using empirical     relations from Kim et al. (2007) and Wang et al. (2016). These estimates are used throughout the analysis pipeline to evaluate X-ray position        accuracy, compute catalog-matching radii, and assess the reliability of potential optical/IR counterparts.

    Global configuration:
    R90: the PSF R_90 radius in arcsec
    counts: net counts inside 95% radius

flux_spectra
    This script processes each FXRT flare candidate by extracting a flare-filtered event file, computing bootstrap timing (T90), deriving Bayesian       hardness ratios, extracting spectra via CIAO, and fitting three physical spectral models (power law, blackbody, bremsstrahlung) in Sherpa. It        computes absorbed/unabsorbed fluxes, fluences, peak flux, and adaptive light curves. Final outputs include full spectral diagnostics, HR             evolution, model residuals, and a machine-readable summary file with all derived physical properties and model selection metrics (W-stat, AIC,       BIC).

    Global configuration:
    nh_cm2_by_obsid: hydrogen column density for specific ObsID / location

