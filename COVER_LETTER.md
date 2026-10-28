*(Draft. Complete the bracketed fields before sending.)*

[Date]

The Editors
*Molecular Ecology Resources*

Dear Editors,

I submit for your consideration **"Minimum-evidence thresholds decide which taxa a metabarcoding survey
can measure, not how specialised they look: why per-taxon habitat preferences do not transfer between
datasets"** as a Resource Article.

**How the work fits the Aims and Scope.** The paper is about a tool-level property of a method the
molecular ecology community uses constantly, not about the ecology of any taxon. Every analysis that tests
taxa for habitat preference must first decide which taxa carry enough evidence to test. We show that this
minimum-evidence rule is not a neutral filter: it admits taxa whose expected read count is near one at
whatever depth is used, so sequencing depth changes *which* taxa a survey can measure rather than how
specialised they appear. The deliverables are therefore methodological: two prerequisites for
cross-dataset comparison (matched read depth and matched taxonomic rank, each with a diagnostic we report
in full), a demonstration that single-draw rarefaction does not control its stated error rate, and a
benchmark showing that no estimator dominates because the ranking flips with how preference is encoded. The
aquatic-fungal dataset is the vehicle, not the subject.

**What is new relative to the closest prior work.** That sequencing depth affects diversity metrics is
established (Shirazi et al., 2021). The closest result to ours is García-Machado et al. (2023), who found
that depth strongly affected species occupancy *before* a minimum-read threshold was applied and had
negligible effect after it, and read the threshold as curing depth sensitivity. We read the same phenomenon
the other way: the threshold conceals depth sensitivity by changing which taxa are measured, so the
stability it buys is stability of the marginal distribution and not of any named taxon. A depth ladder
inside one dataset, with ecosystem, samples, reference database and pairing rule fixed and only depth
varied, separates those readings.

**Pre-registration and reproducibility.** Eleven analysis plans were committed to version control before
the corresponding outcomes were computed; their SHA-256 hashes are recorded outside the plan files and the
analysis scripts assert them before computing. Predictions that failed are reported as failures, including
two that carried an explanation we subsequently withdrew, and one estimator benchmark that was discarded
and replaced because its generator had handed one estimator its own estimand. Code, all frozen plans and
all result tables are at github.com/kgurgen97/aquatic-fungi-habitat; the FunAqua data are archived at
Zenodo under CC-BY and the NEON products are public.

**Declarations.** The manuscript is original, is not under consideration elsewhere, and has not been
published previously. I am the sole author and declare no competing interests. [Funding statement.]
[Suggested reviewers, if the journal invites them.]

Thank you for considering this manuscript.

Yours sincerely,

Gurgen Karapetyan
[Affiliation]
karagurgen23@gmail.com
[ORCID]
