---
name: point-mutation-design
description: Design and generate point-mutant plasmids and primers, especially SnapGene .dna workflows for site-directed mutagenesis, back-to-back PCR, Gibson assembly, and batch mutation tables. Use when asked to infer mutation rules from an example plasmid, choose codons, design F/R mutagenesis primers, update SnapGene features/primers/colors, or validate generated mutant .dna files.
---

# Point Mutation Design

Use this skill for batch point-mutant plasmid construction, especially when the user provides a wild-type SnapGene `.dna` plasmid, a hand-built example mutant, and a mutation list.

## Workflow

1. Inspect the mutation list and any hand-built example before generating files. Infer numbering, file naming, primer naming, visual marking, and feature update conventions from the example.
2. Parse SnapGene files structurally instead of treating them as plain text. The useful chunks are usually: type `0` sequence, type `10` features XML, type `5` primers XML, type `20` strand colors XML, and type `6` notes XML.
3. Confirm coding coordinates. Do not guess protein numbering when a CDS/start feature or example mutant can reveal it. Validate each requested wild-type residue against the wild-type sequence before editing.
4. Select the mutant codon according to the user rule. If no rule is specified, use the project convention or ask. For human expression projects, default to the human optimal codon table in `scripts/snapgene_point_mutagenesis.py`.
5. Preserve the plasmid length for substitutions. Write the mutant codon in lowercase in the plasmid sequence and in both mutagenesis primers when the project uses lowercase for visual checking.
6. Update relevant SnapGene annotations: translated feature sequence, mutation-site base color, new F/R primers, primer binding sites, and modification timestamp. Preserve unrelated features, notes, traces, alignments, and base primers.
7. Produce a primer summary table and run file-level validation before reporting success.

## Primer Rules

Use the user's project rules first. For the SUPREM-style back-to-back PCR plus Gibson workflow:

- Use `partial homology arm + mutant codon + partial homology arm` as the complete overlap between the two primers.
- Complete primer-overlap homology arm length: 15 to 25 bp.
- Prefer GC percent near 50 percent, balanced sequence composition, and at least about 5 bp on both sides of the mutant codon when possible.
- Tm near 60 C is calculated only on the 3 prime segment downstream of the mutant codon in each primer's own orientation. Do not calculate this Tm on the full primer or on the full overlap.
- Extend the 3 prime segment until its SnapGene Tm is near the target. If SnapGene is not available programmatically, use the bundled nearest-neighbor approximation, then verify in SnapGene when possible.
- Prefer a 3 prime terminal G/C clamp when it does not create a worse design.
- For reverse primers, design in reverse-complement orientation and apply the same rule to the segment after the reverse-complemented mutant codon.

## Bundled Script

Use `scripts/snapgene_point_mutagenesis.py` for deterministic SnapGene batch generation. Read or patch it when project-specific naming, coordinates, codon policy, or primer rules differ.

Typical command:

```powershell
python scripts/snapgene_point_mutagenesis.py `
  --wt "0#SUPREM_HS-pEGFP.dna" `
  --template "1#pEGFP-SUPREM-M1-01-W59C.dna" `
  --mutations mu1.txt `
  --output-dir . `
  --series M1 `
  --primer-prefix SU `
  --start-codon 625 `
  --gene-start 628 `
  --gene-end 1665 `
  --gene-feature-id 17 `
  --skip-number 1
```

The script writes generated `.dna` files, a TSV primer/design summary, and `primer_order.tsv` with two columns: primer name and primer sequence. Its Tm calculation is an approximation calibrated for design ranking, not a substitute for SnapGene verification when exact Tm reporting matters.

## Validation Checklist

Before finishing, verify:

- Every requested wild-type residue matches the template before mutation.
- Only the intended codon changes, ignoring intended lowercase marking.
- The translated protein contains the requested amino acid substitution.
- Mutant codons follow the requested codon policy.
- F/R primers have expected names, lowercase mutant codon, valid binding sites, overlap length 15 to 25 bp, and 3 prime segment Tm near target.
- Green or project-specific color marking covers the mutant codon in both strands if the example uses it.
- Existing user files are not overwritten unless the user explicitly requested regeneration.
