# Analysis of inclusive jets in Au+Au collisions at $\sqrt{s_{\textrm{NN}}} = 200$ GeV

## Tree production
`data/` folder contains the jet algorithm and tree production.


### Usage and setup
- The project should be compiled with `64-bit` architecture using command:
```bash
setup 64b
```
- For now the code needs to be run on old rcas60XX machines, one exception is unfolding, which needs to be executed in starsub0X machines

#### FastJet
One does not have to use own `fastjet` installation. There already exists compiled version -> one may check it using
```bash
echo $LD_LIBRARY_PATH
```
and see something like
`/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86/gcc-4.8.5/fastjet-3.3.4-2ro35ixrxr4b5jn4dprn46h3t37n64od/lib` along the libraries.
Another possibility is to install `fastjet` using `data/fastjet_install.sh` script, which will download and compile the `fastjet` library in users `gpfs` directory



For local use - `data/runLocal.sh` to create a test tree
And mass production - `data/run.sh` which accepts the list of input files as an argument.

`./run.sh` runs both data and embedding (same logic for Local running)
`./run.sh embedding` runs only embedding (same logic for Local running)
`./run.sh data` runs only data (same logic for Local running)


### Main components
---
- The main compiled source code is situated in `data/StRoot/StPicoHFJetMaker/StPicoHFJetMaker.{cxx,h}`


- Macro used for running is `data/macros/runPicoHFJetMaker.C`

### Tree structure
The output tree dependent on the input picoDST file (embedding or real data) contains the following branches:

#### Event information:
- `runId`: Run number
- `centrality`: Centrality bin
- `centralityWeight`: Centrality weight

#### Reconstructed jet:
- `reco_pt`: Jet transverse momentum
- `reco_pt_corr`: Jet transverse momentum corrected for background
- `reco_eta`: Jet pseudorapidity
- `reco_phi`: Jet azimuthal angle
- `reco_area`: Jet area
- `reco_rho`: Background density of event
- `reco_pt_lead`: Leading constituent transverse momentum
- `reco_n_constituents`: Number of constituents in the jet
- `reco_neutral_fraction`: Fraction of neutral constituents in the jet
- `reco_trigger_match`: Flag indicating if the jet matches a trigger

#### Monte Carlo jet (if available):
- `xsecWeight`: Cross section weight
- `mc_pt`: Jet transverse momentum
- `mc_eta`: Jet pseudorapidity
- `mc_phi`: Jet azimuthal angle
- `mc_area`: Jet area
- `mc_pt_lead`: Leading constituent transverse momentum
- `mc_n_constituents`: Number of constituents in the jet
- `mc_neutral_fraction`: Fraction of neutral constituents in the jet
- `deltaR`: Delta R between reconstructed and Monte Carlo jet

## Tree merging
- Output from jobs can be merged in `trees/` folder, all of the scripts select the newest job output and merge contents from production
- To merge data, use `./merge_data_all.sh`
- To merge all embedding into one file, use `./merge_embedding_all.sh`
- To separately merge p_T^hat bins, use `./merge_pThatbins.sh`
- To merge different job output, use `./merge_data_all.sh ../data/submit/2025-11-28`, similarly for embedding
- Merged file will be generated in this folder

## Analysis
`analysis/` contains the analysis of the produced trees like filling histograms, unfolding, drawing

### Histogram filling
- Some of the QA histograms can be created by using `histograms/` folder
- Execute `make_hists.C` using `./run_hists.sh`

### Efficiencies

### Unfolding
The unfolding is performed using the [`RooUnfold`](https://gitlab.cern.ch/RooUnfold/RooUnfold) package. The unfolding procedure - `analysis/unfolding/` 