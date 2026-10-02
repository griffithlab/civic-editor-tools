# CIViC Evidence Item Basics

Summarized from CIViC's public documentation (see Sources at the bottom). This is reference
material for an LLM assisting with editorial review — it is deliberately concise, not a full
restatement of the curation guide. When in doubt, the linked pages are authoritative.

## What an evidence item is

An Evidence Item is a manually curated clinical statement, drawn from a single publication, that
describes how a **Molecular Profile** (one or more variants — a single variant in the simple case,
or a more complex combination) affects protein function, oncogenicity, cancer predisposition,
diagnosis, prognosis, or treatment response. Example: "Patients with BRAF V600 mutations respond
well to the drug dabrafenib."

Every evidence item has these required components: Molecular Profile, Source, Variant Origin,
Disease, Evidence Statement, Evidence Type, Evidence Level, Evidence Direction, Significance, and
an Evidence (trust) Rating. Associated Phenotypes and Therapies are populated when relevant.

**A single evidence item should describe one disease and one evidence type.** If a source supports
multiple diseases, or multiple distinct types of claim, CIViC curators create separate evidence
items rather than combining them into one.

## Evidence type

Evidence Type is chosen first, since it determines which Significance values are valid (see
below).

| Type | Meaning |
|---|---|
| **Predictive** | The variant's effect on therapeutic response. |
| **Diagnostic** | The variant's impact on patient diagnosis (cancer subtype). |
| **Prognostic** | The variant's impact on disease progression, severity, or patient survival. |
| **Predisposing** | A *germline* molecular profile's role in conferring susceptibility to disease (including pathogenicity evaluations). |
| **Oncogenic** | A *somatic* variant's involvement in tumor pathogenesis (per the Hallmarks of Cancer). |
| **Functional** | The variant alters biological function relative to the reference state. |

## Evidence direction

Direction states whether the evidence **supports** or **does not support** the claim implied by
the Evidence Type + Significance combination — it does not stand alone. Values: `Supports`,
`Does Not Support`, and (Oncogenic evidence only) `N/A`, for cases where a directional assessment
doesn't apply.

What "Supports"/"Does Not Support" mean shifts with Evidence Type:
- **Predictive** — supports/refutes a drug–variant response relationship.
- **Diagnostic** — supports/refutes the variant's diagnostic impact.
- **Prognostic** — supports/refutes the stated outcome association.
- **Predisposing** — supports/refutes a pathogenic (or protective) germline role.
- **Oncogenic** — supports/refutes an oncogenic (or protective) somatic role.
- **Functional** — supports/refutes the stated functional alteration.

A few specific distinctions the curation guide calls out, since they're easy to get wrong:
- **"Resistance"** (Significance) applies when a variant *actively induces* resistance; if a
  variant merely *fails to sensitize* a patient to a drug, that's `Does Not Support` +
  `Sensitivity/Response`, not `Resistance`.
- **"Reduced Sensitivity"** is a comparison against an already-established *sensitizing* profile
  for that drug — not a comparison between two different drugs.
- For Predisposing/Oncogenic evidence, `Supports` + the positive significance value (Predisposition
  / Oncogenicity) and `Does Not Support` + the same value are used to express the same clinical
  concept in opposite directions (i.e. likely pathogenic vs. likely benign) — don't confuse
  `Does Not Support` here with "no data."

## Evidence level

Ranks how directly the study supports the claim — roughly, how close the evidence is to
"actionable in a real patient today."

| Level | Name | What it means |
|---|---|---|
| **A** | Validated association | Proven/consensus association in human medicine — typically Phase III trials or an FDA-approved companion diagnostic already in routine clinical use. |
| **B** | Clinical evidence | Clinical trial or other primary patient data (generally >5 patients) supports the association; may be a smaller/less definitive Phase I–III trial. |
| **C** | Case study | Individual case report(s) from the clinical literature — fewer than 5 patients. |
| **D** | Preclinical evidence | In vivo or in vitro models (cell lines, animal models, molecular assays) — not direct human clinical data. |
| **E** | Inferential association | Indirect evidence, at least one step removed from a direct variant↔clinical-relevance link (e.g. in silico prediction, inferred mechanism). |

Levels A and B are the highest curation priority since they're the most immediately clinically
applicable.

## Significance

The valid Significance values depend on the Evidence Type chosen:

**Predictive** — `Sensitivity/Response`, `Reduced Sensitivity`, `Resistance`, `Adverse Response`,
or `N/A` (variant lacks clinical interpretive value).

**Diagnostic** — `Positive` (associated with diagnosis of the disease/subtype) or `Negative`
(associated with lack of it).

**Prognostic** — `Better Outcome`, `Poor Outcome`, or `N/A`.

**Predisposing** — `Predisposition` (germline variant with cancer-risk potential, potentially
meeting ACMG/AMP criteria) or `Protectiveness` (germline variant protective against cancer).

**Oncogenic** — `Oncogenicity` (somatic variant with cancer-driving potential, per
ClinGen/CGC/VICC criteria) or `Protectiveness` (somatic variant protective against cancer). Note
oncogenic effects can be cellular-context-dependent, so the Disease field matters here too.

**Functional** — `Gain of Function`, `Loss of Function`, `Unaltered Function`, `Neomorphic`
(novel function), `Dominant Negative` (abrogates the wild-type allele's product), or `Unknown`.

Note the current ACMG/AMP and ClinGen/CGC/VICC scoring criteria referenced above apply only to
*simple, single-variant* molecular profiles, not complex multi-variant ones.

## Variant origin

Whether the variant was inherited (germline) or acquired (somatic) in the study being cited:

- **Somatic** — found only in tumor/diseased cells; not expected to be heritable.
- **Rare Germline** — present in every cell; population frequency <1% for the relevant population.
- **Common Germline** — present in every cell; population frequency ≥1%.
- **Combined** — a complex molecular profile whose component variants have heterogeneous origins.
- **Mixed** — the patient population described mixes somatic and germline cases.
- **Unknown** — origin can't be determined from the available evidence.
- **N/A** — the variant type (e.g. an expression pattern) doesn't fit the origin concept at all.

CIViC prioritizes somatic events and clinically-relevant germline variants; common germline
polymorphisms (≥1% frequency) are lower curation priority unless independently clinically
significant.

## Disease

The cancer type/subtype associated with the molecular profile, **as described in the specific
source being curated** — not the curator's general knowledge of where a variant is relevant.
Drawn from the Disease Ontology (DO). Curation guidance:
- Use the **most specific** DO subtype available, not a broad parent category.
- **One disease per evidence item.** If a source supports multiple diseases, split into multiple
  evidence items.
- If the right term doesn't exist in DO yet, a curator can add it, or submit it to the DO Term
  Tracker.

## Therapy and therapy interaction type

Therapy is required for **Predictive** evidence. Curators use the NCI Thesaurus (NCIT) and prefer
the most specific standardized name (not a trade name or broad drug class) — a novel therapy can
be added if it's genuinely missing from NCIT.

When an evidence item involves more than one therapy, the relationship between them must be one
of, and must be **explicitly stated in the source** (not inferred):
- **Combination** — administered together.
- **Sequential** — given one after another.
- **Substitutes** — used interchangeably.

If a source actually describes more than one of these relationships at once (e.g. some patients
got a combination, others a substitution), consider splitting into separate evidence items rather
than conflating interaction types.

## Associated phenotype

Optional. Drawn from the Human Phenotype Ontology (HPO); standardizes a symptom/abnormality (e.g.
pheochromocytoma) seen in a patient carrying the variant. Populate it only when:
- the cited source **specifically documents** that phenotype in a patient with the variant, and
- it adds information beyond what the Disease field already conveys — most useful for
  Predisposing evidence involving syndromes or non-binary phenotype presentations.

Leave it blank for phenotypes observed across a group of patients with merely *related* (not
identical) variants, or when the association is the curator's general knowledge rather than
something the specific source states.

## Source

Links the evidence item to exactly one publication. CIViC currently accepts:
- **PubMed** (by PMID),
- **ASCO** Meeting Library abstracts (by ASCO Web ID),
- **ASH** (Blood journal) meeting abstracts (by DOI).

Curation should draw from **primary literature, not review articles** — avoid pulling a claim from
a review's introduction/discussion where it's just restating someone else's earlier finding.
Within that constraint, a full published article and an ASCO/ASH meeting abstract are both treated
as legitimate primary sources — the operative distinction is *primary vs. secondary* literature,
not simply *full text vs. abstract*. That said, **when a full article is available, curation should
use it rather than stopping at the abstract**: an abstract alone often omits exactly the detail
(sample size, statistics, disease stage, prior treatment, confounding factors) needed to curate
and review evidence accurately, so treat "abstract only, when a full article exists" as reduced
context, not as equivalent to the full text.

Watch for **overlapping patient populations** — the same patients/outcomes reported in more than
one publication (e.g. an early case series later folded into a larger trial) should be flagged
rather than treated as independent supporting evidence.

## Evidence rating

A 1–5 star rating of curator confidence in *this specific evidence statement*, evaluated in
isolation — not a rating of the source publication's overall quality.

| Stars | Meaning |
|---|---|
| **5** | Strong, well-supported evidence from a lab/journal of respected standing; reproducible, well-controlled, independently confirmed, statistically powerful. |
| **4** | Strong, well-supported evidence; well-controlled with convincing results; any discrepancies are well explained. |
| **3** | Compelling but limited in breadth — e.g. a smaller-scale or novel finding without extensive follow-up validation; discrepancies still adequately explained. |
| **2** | Not well supported — may lack proper controls, adequate sample size, or statistical power; little follow-up data available. |
| **1** | Not well supported — results not reproducible or from a very small sample; no validation of a novel claim. |

## Evidence statement (the `description` field)

The evidence statement is a brief summary of the molecular profile's clinical relevance in the
given disease/evidence context. Good statements typically include: the evidence type, the
molecular profile and gene, the disease, any therapies/comparisons, sample size, key conclusions,
and supporting statistics (p-values, confidence intervals, etc.) — aimed at a general scientific
audience (avoid unexplained field-specific jargon), and written from primary literature rather
than reviews.

Common pitfalls to flag when reviewing a proposed statement or its revision:
- Mixing more than one Evidence Type's worth of claims into a single statement.
- Missing context that changes how the evidence should be interpreted (disease stage, prior
  treatment, other confounders).
- Dropping necessary specifics (sample size, statistics) that let a reader judge evidence strength.
- Including protected health information (PHI).
- Drawing from a secondary/review source instead of the primary literature.
- For clinical trials, omitting the trial name/ID where the source provides one.
- For Functional evidence, omitting the experimental system (cell type/model, expression vector).

## Common revision pitfalls

- A revision that changes `significance` without checking whether `evidence_direction` and
  `evidence_type` still form a coherent, valid combination (see the Direction and Significance
  sections above).
- A revision that changes `disease` or `therapies` beyond what the cited source actually supports.
- A description edit that drifts away from what the source states, or blends in claims that
  belong to a different evidence type.
- Treating an abstract-only source as if it were full-text-equivalent when a fuller article is
  actually available (see Source, above).
- Overlapping-patient-population sources being used as if they were independent corroboration.
- ACMG/AMP or ClinGen/CGC/VICC-style scoring language applied to a complex (multi-variant)
  molecular profile, where those criteria don't apply.

## Sources

- https://docs.civicdb.org/en/latest/model/evidence.html (and its linked sub-pages: type, direction,
  level, significance, origin, evidence_rating, disease, therapy, source, statement,
  molecular_profile, associated_phenotype)
- https://docs.civicdb.org/en/latest/curating/evidence.html
