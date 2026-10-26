# Internal notes, not part of the submission

Split out of `MANUSCRIPT.md` so the submission copy carries only the paper. These record decisions, their reasons and what each pre-registered analysis found, including the ones that failed.

## Author's notes on what still has to happen before submission

1. **References — done, 2026-10-22.** All 14 entries verified against Crossref or the publisher's record;
   author lists, volumes, article numbers and DOIs come from those records. Four corrections came out of
   the check and are now in the text: (a) the Malfertheiner paper is *Community conservatism is widespread
   across microbial phyla and environments*, not an aquatic-specific paper, and the draft's gloss
   "leaving the fungal question open" was replaced by their own sentence, since they never mention fungi;
   (b) § 3.6 had attributed Levins' *B* to both niche-breadth studies — Malard et al. use
   environmental-gradient breadth, and Hou et al. use Levins only to split generalists from specialists;
   (c) Hou et al. is an indoor glass-tank experiment, not a field survey, and is now described as such;
   (d) Reñé et al. was in the list but cited nowhere, and reports the *opposite* vertical pattern, so it
   is now cited in § 3.1 as the disagreement it is. Madin et al. (2020) was removed as belonging to the
   companion study. The Open Tree synthesis version was never recorded, so the query was re-run: synthesis
   `opentree16.1`, all 77 OTT ids identical to the stored file, 70 placed in the tree — the number the
   text reports. NEON was pulled without pinning a release; the data availability statement now says so
   and shows that no provisional month entered the analysis.
2. **Framing — decided, 2026-10-22: mechanism-led.** The biology-led alternative was rejected on the
   paper's own evidence. Its headline would have been that fungi are generalists among eukaryotes, but
   § 3.5 shows fungi (0.150) are indistinguishable from Alveolata (0.157) and that the pooled gap against
   non-fungal eukaryotes is driven by Metazoa (0.294), Straminipila (0.271) and Viridiplantae (0.227); on
   Levins' *B* pooled against protists, p = 0.163. The title would have been contradicted by the paper's
   own table. § 3.5 also has no external validation available — NEON carries only 16S and ITS, and
   metaPR2 failed rank-comparability — whereas the mechanism is derived under known truth and confirmed in
   two independent datasets. The depth objection to § 3.5 does *not* hold, and that is worth keeping in
   mind if a reviewer raises it: the comparison uses the same 204 pairs at the same per-pair depth, the
   ≥ 20 discordant-pair rule applies identically to both groups, and robustness check (ii) leaves the
   difference intact in 4 of 4 discordant-count bins. The mechanism in § 3.4 reinforces this rather than
   undermining it, since both groups are filtered to the same position on the detection curve.
   Consequently: the title, abstract, introduction and § 4.1 now lead with the mechanism; the mechanism
   section was promoted from the interim label 3.3b to § 3.4, shifting the two sections after it; and the
   phrase "fungi sit at the generalist end of the eukaryotic range" was replaced in the abstract, the
   § 3.5 heading and the Figure 3 caption by what the table actually shows.
2a. **External review, 2026-10-23 — what it changed and what it leaves open.** A reviewer read the draft,
   the results files and the figures. Four criticisms were correct and are now fixed in the text:
   (a) § 3.5's heading read "What does transfer", which § 4.4 contradicts — no external validation of that
   section was possible, so nothing showed it transfers; it is now "A group-level contrast within
   FunAqua". (b) § 3.3's heading asserted non-transfer as established, and § 4.1 called the agreement
   "at chance"; with 14 of 27 and a 95 % CI of 34.0–69.3 %, the interval contains both chance and moderate
   agreement, so both now say replication was *not demonstrated* rather than *absent*. (c) The Jaccard
   0.51 turnover between the two datasets is **confounded** — they differ in geography, water-body type
   and in what "sediment" means, so orders testable in one and not the other may be ecologically absent
   rather than filtered out by depth. § 3.4 now says the simulation shows the mechanism is *sufficient* to
   produce such turnover, not that it produced this one. (d) The λ ≈ 1 selection was derived at a single
   effect size (3×) and one abundance distribution; a taxon with a large effect accumulates discordant
   pairs at λ ≫ 1, so the statement holds in the small-to-moderate effect regime, which § 3.4 now says.
   The abstract was rewritten to 241 words, inside Molecular Ecology Resources' 250-word cap.

   Still open, in the order I would do them: (i) the controlled within-NEON depth experiment —
   `rarefied_depth` runs 232 to 425 195, so a ladder of 100/300/1 000/3 000/10 000 keeps 102–149 pairs and
   holds ecosystem, samples and reference database fixed while varying only depth. This is the test (c)
   requires and the one that would turn the mechanism from sufficient to demonstrated. (ii) the parameter
   sweep for (d). (iii) a benchmark of alternative estimators — single and repeated rarefaction,
   abundance-based paired models, hierarchical occupancy, beta-binomial — on simulated truth, comparing
   bias, RMSE, FDR, power and coverage. (iv) sensitivity to site clustering (149 NEON pairs come from
   23 sites) and to taxonomic non-independence among orders, neither of which the agreement test accounts
   for. Also: the reviewer's literature check and ours agree that depth and rarefaction effects on
   diversity metrics are known (Shirazi et al., 2021, *Ecol Evol* 11(22)), and that the closest result to
   ours is García-Machado et al. (2023, *Environmental DNA*), who found depth strongly affected occupancy
   *before* thresholding and not after. They read that as the threshold curing depth sensitivity; we read
   the same phenomenon as the threshold concealing it by changing which taxa are measured. That contrast
   is the paper's clearest novelty claim and should be stated explicitly in the introduction.
3. **Figures — drawn.** `make_figures.py` produces all four as vector PDF and 600 dpi PNG in `figures/`,
   sized to 89 mm (single column) and 183 mm (double). Palette is categorical slots 1–2 of a validated
   default (blue `#2a78d6`, orange `#eb6834`); the validator passes all checks on all pairs, and the
   aqua slot was dropped because it warns below 3:1 against the light surface. Two layout bugs and one
   numerical inconsistency were found by looking at the rendered output and are fixed:
   Fig 3's mirrored layout clipped "Straminipila" to "Straminip", Fig 1 printed a blue annotation on top
   of blue bars, and Fig 1 counted stability against the draws in which an order was *testable* rather
   than against all 200, which made it say 36 where the text says 35. An assertion in the script now
   fails if the figure and `draw_stability.json` ever disagree again.

   | figure | panel | file |
   |---|---|---|
   | 1 | per-order significance frequency across 200 draws, 35-order stable set marked | `fig1_draw_stability` |
   | 2 | FunAqua against NEON preference index, 27 orders, agreement quadrants | `fig2_transfer_failure` |
   | 3 | median effect magnitude and Levins *B* by kingdom, two panels | `fig3_kingdoms` |
   | 5 | the algebra and the idealised simulation, three panels | `fig4_simulation` |
   | 6 | the mechanism in real data: ladder, sweep, four abundance variants | `fig5_mechanism_in_data` |

3a. **Why § 3.4a–3.4c keep letters when § 3.3b did not.** These look like the interim label that was
   removed earlier, and they are not the same case. "3.3b" was an independent section wedged after § 3.3
   while it waited for a number. § 3.4a–3.4c are three parts of the single argument § 3.4 opens — the
   algebra, then the same question in real data, then its generality and the model that reproduces it —
   which is an ordinary convention. They are also referenced by `SWEEP_PLAN.md`, `PREVALENCE_PLAN.md` and
   `BENCHMARK_PLAN.md`, which are frozen pre-registrations and must not be edited; renumbering the
   manuscript would leave those references pointing at sections that no longer exist and weaken the
   provenance trail they exist to provide. A journal copy-editor may renumber at proof stage, which is
   the right moment for it.

   Still to decide: whether Fig 4 belongs in the main text or the supplement, given that it is post hoc
   and reaches significance at no cap.
4. **Abstract word count** is 248; check the target journal's limit.
5. **Statement on the corrected numbers.** § 3.2 reports that our own two implementations disagreed. Decide
   whether that belongs in the Results as written, or in a methods note. It is honest either way and a
   reviewer will respect it, but it does draw attention.
