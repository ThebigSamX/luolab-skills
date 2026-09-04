from __future__ import annotations

import argparse
import csv
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET


WT_FILE = Path("0#SUPREM_HS-pEGFP.dna")
TEMPLATE_FILE = Path("1#pEGFP-SUPREM-M1-01-W59C.dna")
MUTATION_FILE = Path("mu1.txt")
OUTPUT_DIR = Path(".")
SUMMARY_FILE = Path("point_mutation_primers.tsv")
PRIMER_ORDER_FILE = Path("primer_order.tsv")
OUTPUT_PATTERN = "{n}#pEGFP-SUPREM-{series}-{nn}-{mutation}.dna"
SERIES = "M1"
PRIMER_PREFIX = "SU"
SKIP_NUMBERS: set[int] = set()

START_CODON_1BASED = 625
SUPREM_FEATURE_START_1BASED = 628
SUPREM_FEATURE_END_1BASED = 1665
GENE_FEATURE_ID = 17

CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

# Human-preferred codons used for the requested substitutions.
HUMAN_OPTIMAL_CODON = {
    "A": "GCC",
    "C": "TGC",
    "D": "GAC",
    "E": "GAG",
    "F": "TTC",
    "G": "GGC",
    "H": "CAC",
    "I": "ATC",
    "K": "AAG",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCC",
    "Q": "CAG",
    "R": "CGC",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAC",
}

AA_MASS = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "Q": 128.1307,
    "E": 129.1155,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
}

NN_PARAMS = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

REVCOMP = str.maketrans("ACGTacgt", "TGCAtgca")


def read_chunks(path: Path) -> list[tuple[int, bytes]]:
    data = path.read_bytes()
    chunks: list[tuple[int, bytes]] = []
    offset = 0
    while offset + 5 <= len(data):
        chunk_type = data[offset]
        length = int.from_bytes(data[offset + 1 : offset + 5], "big")
        payload = data[offset + 5 : offset + 5 + length]
        if len(payload) != length:
            raise ValueError(f"{path}: truncated chunk at offset {offset}")
        chunks.append((chunk_type, payload))
        offset += 5 + length
    if offset != len(data):
        raise ValueError(f"{path}: trailing bytes after chunks")
    return chunks


def write_chunks(path: Path, chunks: Iterable[tuple[int, bytes]]) -> None:
    out = bytearray()
    for chunk_type, payload in chunks:
        out.append(chunk_type)
        out.extend(len(payload).to_bytes(4, "big"))
        out.extend(payload)
    path.write_bytes(bytes(out))


def get_sequence(chunks: list[tuple[int, bytes]]) -> tuple[int, str]:
    for chunk_type, payload in chunks:
        if chunk_type == 0:
            return payload[0], payload[1:].decode("ascii")
    raise ValueError("sequence chunk not found")


def get_xml(chunks: list[tuple[int, bytes]], chunk_type: int) -> str:
    for current_type, payload in chunks:
        if current_type == chunk_type:
            return payload.decode("utf-8")
    raise ValueError(f"XML chunk {chunk_type} not found")


def with_xml_declaration(element: ET.Element) -> bytes:
    return b'<?xml version="1.0"?>' + ET.tostring(element, encoding="utf-8")


def gc_percent(seq: str) -> float:
    return 100.0 * sum(base.upper() in {"G", "C"} for base in seq) / len(seq)


def revcomp(seq: str) -> str:
    return seq.translate(REVCOMP)[::-1]


def translate(seq: str) -> str:
    upper = seq.upper()
    if len(upper) % 3:
        raise ValueError("sequence length is not divisible by 3")
    return "".join(CODON_TABLE[upper[i : i + 3]] for i in range(0, len(upper), 3))


def protein_mw(sequence: str) -> float:
    return sum(AA_MASS[aa] for aa in sequence) + 18.0152


def tm_nearest_neighbor(seq: str, sodium_equiv_molar: float = 0.25, oligo_molar: float = 250e-9) -> float:
    seq = seq.upper()
    if len(seq) < 2:
        return 0.0
    delta_h = 0.2
    delta_s = -5.7
    for terminal in (seq[0], seq[-1]):
        if terminal in {"A", "T"}:
            delta_h += 2.2
            delta_s += 6.9
        else:
            delta_h += 0.1
            delta_s -= 2.8
    for first, second in zip(seq, seq[1:]):
        h, s = NN_PARAMS[first + second]
        delta_h += h
        delta_s += s
    return (
        delta_h * 1000.0 / (delta_s + 1.987 * math.log(oligo_molar / 4.0))
        - 273.15
        + 16.6 * math.log10(sodium_equiv_molar)
    )


def parse_mutations(path: Path) -> list[tuple[int, str, str, int, str]]:
    rows: list[tuple[int, str, str, int, str]] = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line or line.startswith("Id"):
            continue
        identifier, mutation = re.split(r"\s+", line)
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", mutation)
        if not match:
            raise ValueError(f"bad mutation entry: {line}")
        wt_aa, pos_text, new_aa = match.groups()
        m1_number = int(identifier.split("-")[1])
        rows.append((m1_number, wt_aa, int(pos_text), new_aa, mutation))
    return rows


def choose_overlap(seq: str, codon_start0: int) -> tuple[int, int, str]:
    best: tuple[float, int, int, str] | None = None
    for overlap_len in range(15, 26):
        for left_len in range(5, overlap_len - 7):
            right_len = overlap_len - 3 - left_len
            overlap_start0 = codon_start0 - left_len
            overlap_end0 = codon_start0 + 3 + right_len
            if overlap_start0 < 0 or overlap_end0 > len(seq):
                continue
            overlap = seq[overlap_start0:overlap_end0]
            score = (
                abs(gc_percent(overlap) - 50.0)
                + 0.20 * abs(overlap_len - 20)
                + 0.12 * abs(left_len - 10)
            )
            candidate = (score, overlap_start0, overlap_end0, overlap)
            if best is None or candidate < best:
                best = candidate
    if best is None:
        raise ValueError("could not choose overlap")
    _, start0, end0, overlap = best
    return start0, end0, overlap


def choose_three_prime_extension(
    seq: str,
    segment_after_codon_start0: int,
    fixed_segment_end0: int,
    direction: str,
) -> tuple[int, str, float]:
    best: tuple[float, int, str, float] | None = None
    if direction == "forward":
        for extension_len in range(0, 45):
            segment = seq[segment_after_codon_start0 : fixed_segment_end0 + extension_len]
            if len(segment) < 8 or fixed_segment_end0 + extension_len > len(seq):
                continue
            tm = tm_nearest_neighbor(segment)
            clamp_penalty = 0.0 if segment[-1].upper() in {"G", "C"} else 1.25
            len_penalty = 0.05 * max(0, len(segment) - 30)
            score = abs(tm - 60.0) + clamp_penalty + len_penalty
            candidate = (score, extension_len, segment, tm)
            if best is None or candidate < best:
                best = candidate
    elif direction == "reverse":
        overlap_start0 = fixed_segment_end0
        codon_start0 = segment_after_codon_start0
        for extension_len in range(0, 45):
            start0 = overlap_start0 - extension_len
            if start0 < 0:
                continue
            segment = revcomp(seq[start0:codon_start0])
            if len(segment) < 8:
                continue
            tm = tm_nearest_neighbor(segment)
            clamp_penalty = 0.0 if segment[-1].upper() in {"G", "C"} else 1.25
            len_penalty = 0.05 * max(0, len(segment) - 30)
            score = abs(tm - 60.0) + clamp_penalty + len_penalty
            candidate = (score, extension_len, segment, tm)
            if best is None or candidate < best:
                best = candidate
    else:
        raise ValueError(f"unknown direction {direction}")
    if best is None:
        raise ValueError("could not choose 3' extension")
    _, extension_len, segment, tm = best
    return extension_len, segment, tm


def design_primers(seq: str, codon_start0: int) -> dict[str, object]:
    overlap_start0, overlap_end0, overlap = choose_overlap(seq, codon_start0)
    f_extension_len, f_tm_segment, f_tm = choose_three_prime_extension(
        seq,
        codon_start0 + 3,
        overlap_end0,
        "forward",
    )
    r_extension_len, r_tm_segment, r_tm = choose_three_prime_extension(
        seq,
        codon_start0,
        overlap_start0,
        "reverse",
    )
    f_binding_start0 = overlap_start0
    f_binding_end0 = overlap_end0 + f_extension_len
    r_binding_start0 = overlap_start0 - r_extension_len
    r_binding_end0 = overlap_end0
    f_primer = seq[f_binding_start0:f_binding_end0]
    r_primer = revcomp(seq[r_binding_start0:r_binding_end0])
    return {
        "overlap_start0": overlap_start0,
        "overlap_end0": overlap_end0,
        "overlap": overlap,
        "f_primer": f_primer,
        "r_primer": r_primer,
        "f_binding_start0": f_binding_start0,
        "f_binding_end0": f_binding_end0,
        "r_binding_start0": r_binding_start0,
        "r_binding_end0": r_binding_end0,
        "f_tm_segment": f_tm_segment,
        "r_tm_segment": r_tm_segment,
        "f_tm_3prime": f_tm,
        "r_tm_3prime": r_tm,
        "f_tm_full": tm_nearest_neighbor(f_primer),
        "r_tm_full": tm_nearest_neighbor(r_primer),
        "overlap_gc": gc_percent(overlap),
    }


def update_features_xml(features_xml: str, suprem_translation: str) -> bytes:
    root = ET.fromstring(features_xml)
    for feature in root.findall("Feature"):
        if feature.attrib.get("recentID") == str(GENE_FEATURE_ID):
            feature.attrib["translationMW"] = f"{protein_mw(suprem_translation):.2f}"
            feature.attrib["consecutiveNumberingStartsFrom"] = "2"
            for qualifier in feature.findall("Q"):
                if qualifier.attrib.get("name") == "translation":
                    value = qualifier.find("V")
                    if value is None:
                        value = ET.SubElement(qualifier, "V")
                    value.attrib["text"] = suprem_translation
                    value.text = None
                    break
            else:
                qualifier = ET.SubElement(feature, "Q", {"name": "translation"})
                ET.SubElement(qualifier, "V", {"text": suprem_translation})
            break
    else:
        raise ValueError("SUPREM feature recentID 17 not found")
    return with_xml_declaration(root)


def primer_element(
    recent_id: int,
    name: str,
    sequence: str,
    start0: int,
    end0_exclusive: int,
    bound_strand: str,
    tm_full: float,
    date_added: str,
) -> ET.Element:
    location = f"{start0}-{end0_exclusive - 1}"
    tm_value = str(round(tm_full))
    primer = ET.Element(
        "Primer",
        {
            "recentID": str(recent_id),
            "name": name,
            "sequence": sequence,
            "description": "<html><body></body></html>",
            "dateAdded": date_added,
        },
    )
    for simplified in (False, True):
        attrs = {
            "location": location,
            "boundStrand": bound_strand,
            "annealedBases": sequence,
            "meltingTemperature": tm_value,
        }
        if simplified:
            attrs = {"simplified": "1", **attrs}
        binding = ET.SubElement(primer, "BindingSite", attrs)
        ET.SubElement(binding, "Component", {"hybridizedRange": location, "bases": sequence})
    return primer


def update_primers_xml(primers_xml: str, m1_number: int, design: dict[str, object]) -> bytes:
    root = ET.fromstring(primers_xml)
    for primer in list(root.findall("Primer")):
        name = primer.attrib.get("name", "")
        if name.startswith(f"{PRIMER_PREFIX}-{SERIES}-"):
            root.remove(primer)
    root.attrib.pop("recycledIDs", None)
    root.attrib["nextValidID"] = "10"
    date_added = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    root.append(
        primer_element(
            8,
            f"{PRIMER_PREFIX}-{SERIES}-{m1_number:02d}F",
            str(design["f_primer"]),
            int(design["f_binding_start0"]),
            int(design["f_binding_end0"]),
            "0",
            float(design["f_tm_full"]),
            date_added,
        )
    )
    root.append(
        primer_element(
            9,
            f"{PRIMER_PREFIX}-{SERIES}-{m1_number:02d}R",
            str(design["r_primer"]),
            int(design["r_binding_start0"]),
            int(design["r_binding_end0"]),
            "1",
            float(design["r_tm_full"]),
            date_added,
        )
    )
    return with_xml_declaration(root)


def update_strand_colors_xml(colors_xml: str, codon_start0: int) -> bytes:
    root = ET.fromstring(colors_xml)
    mutation_range = f"{codon_start0}..{codon_start0 + 2}"
    for strand_name in ("TopStrand", "BottomStrand"):
        strand = root.find(strand_name)
        if strand is None:
            strand = ET.SubElement(root, strand_name)
        for color_range in list(strand.findall("ColorRange")):
            if color_range.attrib.get("colors") == "green":
                strand.remove(color_range)
        ET.SubElement(strand, "ColorRange", {"range": mutation_range, "colors": "green"})
        ranges = list(strand.findall("ColorRange"))
        ranges.sort(key=lambda item: int(item.attrib["range"].split("..")[0]))
        strand[:] = ranges
    return ET.tostring(root, encoding="utf-8")


def update_notes_xml(notes_xml: str) -> bytes:
    root = ET.fromstring(notes_xml)
    last_modified = root.find("LastModified")
    if last_modified is not None:
        now = datetime.now()
        last_modified.attrib["UTC"] = f"{now.hour}:{now.minute}:{now.second}"
        last_modified.text = f"{now.year}.{now.month}.{now.day}"
    return ET.tostring(root, encoding="utf-8") + b"\n"


def replace_chunks(
    template_chunks: list[tuple[int, bytes]],
    sequence_flag: int,
    sequence: str,
    features_payload: bytes,
    primers_payload: bytes,
    colors_payload: bytes,
    notes_payload: bytes,
) -> list[tuple[int, bytes]]:
    replaced: list[tuple[int, bytes]] = []
    for chunk_type, payload in template_chunks:
        if chunk_type == 0:
            replaced.append((chunk_type, bytes([sequence_flag]) + sequence.encode("ascii")))
        elif chunk_type == 10:
            replaced.append((chunk_type, features_payload))
        elif chunk_type == 5:
            replaced.append((chunk_type, primers_payload))
        elif chunk_type == 20:
            replaced.append((chunk_type, colors_payload))
        elif chunk_type == 6:
            replaced.append((chunk_type, notes_payload))
        else:
            replaced.append((chunk_type, payload))
    return replaced



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate SnapGene point-mutant plasmids and mutagenesis primers."
    )
    parser.add_argument("--wt", default=str(WT_FILE), help="Wild-type SnapGene .dna file")
    parser.add_argument("--template", default=str(TEMPLATE_FILE), help="Example/template SnapGene .dna file to preserve metadata style")
    parser.add_argument("--mutations", default=str(MUTATION_FILE), help="Mutation table with Id and Mut. columns")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Directory for generated .dna files")
    parser.add_argument("--summary", default=str(SUMMARY_FILE), help="TSV summary path or filename")
    parser.add_argument("--order-summary", default=str(PRIMER_ORDER_FILE), help="Two-column primer order TSV path or filename")
    parser.add_argument("--output-pattern", default=OUTPUT_PATTERN, help="Python format pattern using {n}, {nn}, {series}, {mutation}")
    parser.add_argument("--series", default=SERIES, help="Mutation series label, e.g. M1")
    parser.add_argument("--primer-prefix", default=PRIMER_PREFIX, help="Primer prefix, e.g. SU")
    parser.add_argument("--start-codon", type=int, default=START_CODON_1BASED, help="1-based coordinate of the ATG start codon")
    parser.add_argument("--gene-start", type=int, default=SUPREM_FEATURE_START_1BASED, help="1-based start of the translated feature excluding start codon if appropriate")
    parser.add_argument("--gene-end", type=int, default=SUPREM_FEATURE_END_1BASED, help="1-based inclusive end of the translated feature")
    parser.add_argument("--gene-feature-id", type=int, default=GENE_FEATURE_ID, help="SnapGene Feature recentID to update")
    parser.add_argument("--skip-number", type=int, action="append", default=[], help="Mutation number to skip; repeatable")
    return parser.parse_args()
def main() -> None:
    global WT_FILE, TEMPLATE_FILE, MUTATION_FILE, OUTPUT_DIR, SUMMARY_FILE, PRIMER_ORDER_FILE, OUTPUT_PATTERN
    global SERIES, PRIMER_PREFIX, START_CODON_1BASED, SUPREM_FEATURE_START_1BASED
    global SUPREM_FEATURE_END_1BASED, GENE_FEATURE_ID, SKIP_NUMBERS

    args = parse_args()
    WT_FILE = Path(args.wt)
    TEMPLATE_FILE = Path(args.template)
    MUTATION_FILE = Path(args.mutations)
    OUTPUT_DIR = Path(args.output_dir)
    SUMMARY_FILE = Path(args.summary)
    PRIMER_ORDER_FILE = Path(args.order_summary)
    OUTPUT_PATTERN = args.output_pattern
    SERIES = args.series
    PRIMER_PREFIX = args.primer_prefix
    START_CODON_1BASED = args.start_codon
    SUPREM_FEATURE_START_1BASED = args.gene_start
    SUPREM_FEATURE_END_1BASED = args.gene_end
    GENE_FEATURE_ID = args.gene_feature_id
    SKIP_NUMBERS = set(args.skip_number)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    wt_chunks = read_chunks(WT_FILE)
    template_chunks = read_chunks(TEMPLATE_FILE)
    sequence_flag, wt_seq = get_sequence(wt_chunks)
    template_features = get_xml(template_chunks, 10)
    template_primers = get_xml(template_chunks, 5)
    wt_colors = get_xml(wt_chunks, 20)
    template_notes = get_xml(template_chunks, 6)
    mutations = parse_mutations(MUTATION_FILE)
    summary_rows: list[dict[str, object]] = []
    primer_order_rows: list[dict[str, object]] = []

    for m1_number, wt_aa, aa_position, new_aa, mutation_label in mutations:
        if m1_number in SKIP_NUMBERS:
            continue
        codon_start0 = START_CODON_1BASED - 1 + (aa_position - 1) * 3
        wt_codon = wt_seq[codon_start0 : codon_start0 + 3].upper()
        observed_wt_aa = CODON_TABLE[wt_codon]
        if observed_wt_aa != wt_aa:
            raise ValueError(
                f"{mutation_label}: expected {wt_aa} at AA {aa_position}, "
                f"found {observed_wt_aa} from codon {wt_codon}"
            )
        new_codon = HUMAN_OPTIMAL_CODON[new_aa]
        mutated_seq = wt_seq[:codon_start0] + new_codon.lower() + wt_seq[codon_start0 + 3 :]
        full_cds = mutated_seq[START_CODON_1BASED - 1 : SUPREM_FEATURE_END_1BASED]
        full_translation = translate(full_cds)
        if full_translation[aa_position - 1] != new_aa:
            raise ValueError(f"{mutation_label}: translation did not produce requested mutation")
        suprem_translation = translate(mutated_seq[SUPREM_FEATURE_START_1BASED - 1 : SUPREM_FEATURE_END_1BASED])
        design = design_primers(mutated_seq, codon_start0)
        features_payload = update_features_xml(template_features, suprem_translation)
        primers_payload = update_primers_xml(template_primers, m1_number, design)
        colors_payload = update_strand_colors_xml(wt_colors, codon_start0)
        notes_payload = update_notes_xml(template_notes)
        output_name = OUTPUT_PATTERN.format(n=m1_number, nn=f"{m1_number:02d}", series=SERIES, mutation=mutation_label)
        forward_primer_name = f"{PRIMER_PREFIX}-{SERIES}-{m1_number:02d}F"
        reverse_primer_name = f"{PRIMER_PREFIX}-{SERIES}-{m1_number:02d}R"
        output_chunks = replace_chunks(
            template_chunks,
            sequence_flag,
            mutated_seq,
            features_payload,
            primers_payload,
            colors_payload,
            notes_payload,
        )
        write_chunks(OUTPUT_DIR / output_name, output_chunks)
        summary_rows.append(
            {
                "file": output_name,
                "id": f"{SERIES}-{m1_number:02d}",
                "mutation": mutation_label,
                "aa_position": aa_position,
                "codon_1based": f"{codon_start0 + 1}-{codon_start0 + 3}",
                "wt_codon": wt_codon,
                "mut_codon_human_optimal": new_codon,
                "overlap": design["overlap"],
                "overlap_len": len(str(design["overlap"])),
                "overlap_gc": f"{float(design['overlap_gc']):.1f}",
                "forward_primer": design["f_primer"],
                "reverse_primer": design["r_primer"],
                "forward_3prime_tm_segment": design["f_tm_segment"],
                "reverse_3prime_tm_segment": design["r_tm_segment"],
                "forward_3prime_tm": f"{float(design['f_tm_3prime']):.1f}",
                "reverse_3prime_tm": f"{float(design['r_tm_3prime']):.1f}",
                "forward_full_tm_written": round(float(design["f_tm_full"])),
                "reverse_full_tm_written": round(float(design["r_tm_full"])),
            }
        )
        primer_name_col = "\u5f15\u7269\u540d"
        primer_sequence_col = "\u5f15\u7269\u5e8f\u5217"
        primer_order_rows.extend(
            [
                {primer_name_col: forward_primer_name, primer_sequence_col: design["f_primer"]},
                {primer_name_col: reverse_primer_name, primer_sequence_col: design["r_primer"]},
            ]
        )

    summary_path = SUMMARY_FILE if SUMMARY_FILE.is_absolute() else OUTPUT_DIR / SUMMARY_FILE
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    order_path = PRIMER_ORDER_FILE if PRIMER_ORDER_FILE.is_absolute() else OUTPUT_DIR / PRIMER_ORDER_FILE
    primer_name_col = "\u5f15\u7269\u540d"
    primer_sequence_col = "\u5f15\u7269\u5e8f\u5217"
    with order_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=[primer_name_col, primer_sequence_col], delimiter="\t")
        writer.writeheader()
        writer.writerows(primer_order_rows)


if __name__ == "__main__":
    main()




