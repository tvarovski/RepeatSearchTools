from __future__ import annotations
from typing import TypedDict
import regex as re

class RepeatRecord(TypedDict):
    """Position-aware repeat metadata for a single repeat pair."""

    query_seq: str
    match_seq: str
    query_start: int
    query_end: int
    match_start: int
    match_end: int
    is_perfect: bool
    length: int
    inverted: bool


def extractSeqFromFileToList(file_path: str) -> list[list[str]]:
    """Read FASTA records into a nested list.

    Args:
        file_path: Path to a FASTA file.

    Returns:
        A list in the form [[title_0, sequence_0], [title_1, sequence_1], ...].
    """
    with open(file_path, "r", encoding="utf-8") as fasta_file:
        contents = fasta_file.readlines()

    fasta_list: list[list[str]] = []
    for line in contents:
        if ">" in line:
            fasta_list.append([line.strip(">").strip(), ""])
        else:
            fasta_list[-1][1] = fasta_list[-1][1] + line.strip()

    print(f"Extraction of sequence information from {file_path} finished.")
    return fasta_list

def validateDNASquence(sequence: str) -> bool:
    """Check whether a sequence contains only expected nucleotide symbols.

    Args:
        sequence: Candidate DNA sequence.

    Returns:
        True when all characters are canonical DNA symbols (plus gap), else False.
    """
    bases = "ATGCatgc-"
    for base in sequence:
        if base not in bases:
            print("warning, sequence doesn't contain canonical nucleotides")
            return False
    return True

def compl(base: str) -> str:
    """Return the complementary nucleotide for a single base.

    Args:
        base: One nucleotide character.

    Returns:
        Complementary nucleotide character.

    Raises:
        ValueError: If the input base is unsupported.
    """
    complements = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "-": "-",
    }
    if base not in complements:
        raise ValueError(f"Unsupported base '{base}'.")
    return complements[base]

def rev_compl(seq: str) -> str:
    """Return reverse complement of a DNA sequence.

    Args:
        seq: Input DNA sequence.

    Returns:
        Reverse-complemented sequence.
    """
    return "".join(compl(base) for base in reversed(seq))

def searchSequenceForRepeats(
    sequence: str,
    min_query_length: int = 4,
    max_query_length: int = 25,
    min_spacer: int = 0,
    window_size: int = 250,
    imperfect_homology: bool = False,
    min_homology: float = 0.8,
    fixed_errors: int | bool = False,
    inverted: bool = True,
) -> list[RepeatRecord]:
    """Search a full DNA sequence and retain exact coordinate information.

    Args:
        sequence: Full DNA sequence.
        min_query_length: Minimum repeat length to evaluate.
        max_query_length: Maximum repeat length to evaluate.
        min_spacer: Minimum spacing between query and match.
        window_size: Sliding-window size used during scanning.
        imperfect_homology: Whether to allow approximate matching.
        min_homology: Minimum homology for approximate matching.
        fixed_errors: Explicit max edit distance for approximate matching.
        inverted: Whether to search for inverted (True) or direct (False) repeats.

    Returns:
        Deduplicated list of position-aware repeat records.
    """
    if not validateDNASquence(sequence):
        print("Sequence validation failed. Please check your input.")
        return []

    if imperfect_homology:
        print("Search has been set to find imperfect repeats")
        if fixed_errors is not False:
            print(f"    Allowing up to {fixed_errors} errors/mismatches...")
        else:
            print(f"    Searching with a minimum of {min_homology} homology")
    else:
        print("Search has been set to find perfect repeats")

    all_repeats: list[RepeatRecord] = []
    seen: set[tuple[int, int, int, int]] = set()
    seq_len = len(sequence)

    for query_length in range(min_query_length, max_query_length + 1):
        print(f"  Searching {query_length}bp repeats with positions...")
        for i in range(seq_len):
            window = sequence[i : i + window_size]
            for record in findRepeatsWithPositions(
                window,
                i,
                query_length,
                min_spacer=min_spacer,
                imperfect_homology=imperfect_homology,
                min_homology=min_homology,
                fixed_errors=fixed_errors,
                inverted=inverted,
            ):
                key = (
                    record["query_start"],
                    record["query_end"],
                    record["match_start"],
                    record["match_end"],
                )
                if key not in seen:
                    seen.add(key)
                    all_repeats.append(record)

    print(f"Found {len(all_repeats)} unique repeat pairs.")
    return all_repeats

def findRepeatsWithPositions(
    window: str,
    window_start: int,
    query_length: int,
    min_spacer: int = 0,
    imperfect_homology: bool = False,
    min_homology: float = 0.8,
    fixed_errors: int | bool = False,
    inverted: bool = True,
) -> list[RepeatRecord]:
    """Find repeats in a local window and map coordinates to full sequence.

    Args:
        window: Sequence window.
        window_start: Start offset of this window in the full sequence.
        query_length: Length of repeat motif.
        min_spacer: Minimum spacing between query and match.
        imperfect_homology: Whether to allow approximate matching.
        min_homology: Minimum homology fraction for approximate matching.
        fixed_errors: Explicit max edit distance for approximate matching.
        inverted: Whether to search for inverted or direct repeats.

    Returns:
        List of position-aware repeat records.
    """
    query_string = window[:query_length]
    query = rev_compl(query_string) if inverted else query_string
    search_seq = window[query_length + min_spacer :]
    results: list[RepeatRecord] = []

    if imperfect_homology:
        errors = round(len(query) * (1 - min_homology))
        if fixed_errors is not False:
            errors = int(fixed_errors)

    for i in range(len(search_seq) - query_length + 1):
        candidate = search_seq[i : query_length + i]

        if imperfect_homology:
            is_match = bool(re.findall(f"({query}){{e<={errors}}}", candidate))
        else:
            is_match = candidate == query

        if is_match:
            results.append(
                {
                    "query_seq": query_string,
                    "match_seq": candidate,
                    "query_start": window_start,
                    "query_end": window_start + query_length,
                    "match_start": window_start + query_length + min_spacer + i,
                    "match_end": window_start + query_length + min_spacer + i + query_length,
                    "is_perfect": candidate == query,
                    "length": query_length,
                    "inverted": inverted,
                }
            )

    return results

def filterRedundantRepeats(repeats: list[RepeatRecord]) -> list[RepeatRecord]:
    """Remove repeats fully contained by longer already-kept repeats.

    Args:
        repeats: List of position-aware repeat records.

    Returns:
        Filtered repeats sorted by descending repeat length.
    """
    sorted_repeats = sorted(repeats, key=lambda record: record["length"], reverse=True)
    non_redundant: list[RepeatRecord] = []

    for record in sorted_repeats:
        redundant = any(
            kept["query_start"] <= record["query_start"]
            and record["query_end"] <= kept["query_end"]
            and kept["match_start"] <= record["match_start"]
            and record["match_end"] <= kept["match_end"]
            for kept in non_redundant
        )
        if not redundant:
            non_redundant.append(record)

    print(f"Reduced from {len(repeats)} to {len(non_redundant)} non-redundant repeats.")
    return non_redundant
