from __future__ import annotations
from functools import lru_cache
from typing import Literal, TypedDict
import regex as re

RepeatType = Literal["direct", "inverted", "both"]
_COMPLEMENT_TABLE = str.maketrans("ATGCatgc-", "TACGtacg-")

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
    repeat_type: Literal["direct", "inverted"]
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

@lru_cache(maxsize=100_000)
def rev_compl(seq: str) -> str:
    """Return reverse complement of a DNA sequence.

    Args:
        seq: Input DNA sequence.

    Returns:
        Reverse-complemented sequence.
    """
    return seq.translate(_COMPLEMENT_TABLE)[::-1]

def normalizeRepeatType(repeat_type: RepeatType | bool) -> RepeatType:
    """Normalize repeat-type input while accepting the previous boolean form.

    Args:
        repeat_type: `"direct"`, `"inverted"`, `"both"`, or legacy boolean.

    Returns:
        Normalized repeat-type string.

    Raises:
        ValueError: If the repeat type is unsupported.
    """
    if repeat_type is True:
        return "inverted"
    if repeat_type is False:
        return "direct"

    normalized = repeat_type.lower()
    if normalized not in {"direct", "inverted", "both"}:
        raise ValueError('repeat_type must be "direct", "inverted", or "both".')
    return normalized

def findInteriorFuzzyMatch(
    query: str,
    candidate_region: str,
    pattern: re.Pattern,
    min_candidate_length: int,
    max_candidate_length: int,
) -> str | None:
    """Find a fuzzy homology match while rejecting terminal-base errors.

    Args:
        query: Expected direct or reverse-complemented repeat sequence.
        candidate_region: Candidate search region from the search window.
        pattern: Compiled fuzzy regex pattern for the query.
        min_candidate_length: Shortest candidate span to evaluate.
        max_candidate_length: Longest candidate span to evaluate.

    Returns:
        Matched candidate sequence when it is within the configured error
        tolerance and both terminal bases match exactly; otherwise None.
    """
    if not query or not candidate_region:
        return None

    best_match: tuple[int, int, str] | None = None
    max_length = min(max_candidate_length, len(candidate_region))
    for candidate_length in range(min_candidate_length, max_length + 1):
        candidate = candidate_region[:candidate_length]
        if query[0] != candidate[0] or query[-1] != candidate[-1]:
            continue

        match = pattern.fullmatch(candidate)
        if match is None:
            continue

        fuzzy_counts = getattr(match, "fuzzy_counts", (0, 0, 0))
        error_count = sum(fuzzy_counts)
        length_delta = abs(len(query) - candidate_length)
        match_score = (error_count, length_delta, candidate)
        if best_match is None or match_score < best_match:
            best_match = match_score

    if best_match is None:
        return None
    return best_match[2]

def searchSequenceForRepeats(
    sequence: str,
    min_query_length: int = 4,
    max_query_length: int = 25,
    min_spacer: int = 0,
    window_size: int = 250,
    imperfect_homology: bool = False,
    min_homology: float = 0.8,
    fixed_errors: int | bool = False,
    repeat_type: RepeatType | bool = "inverted",
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
        repeat_type: Which repeat orientation to search: `"direct"`, `"inverted"`,
            or `"both"`. Legacy booleans are accepted (`True` = inverted,
            `False` = direct).

    Returns:
        Deduplicated list of position-aware repeat records.
    """
    if not validateDNASquence(sequence):
        print("Sequence validation failed. Please check your input.")
        return []

    repeat_type = normalizeRepeatType(repeat_type)

    if imperfect_homology:
        print("Search has been set to find imperfect repeats")
        if fixed_errors is not False:
            print(f"    Allowing up to {fixed_errors} errors/mismatches...")
        else:
            print(f"    Searching with a minimum of {min_homology} homology")
    else:
        print("Search has been set to find perfect repeats")
    print(f"Repeat type: {repeat_type}")

    all_repeats: list[RepeatRecord] = []
    seen: set[tuple[int, int, int, int, str]] = set()
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
                repeat_type=repeat_type,
            ):
                key = (
                    record["query_start"],
                    record["query_end"],
                    record["match_start"],
                    record["match_end"],
                    record["repeat_type"],
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
    repeat_type: RepeatType | bool = "inverted",
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
        repeat_type: Which repeat orientation to search: `"direct"`, `"inverted"`,
            or `"both"`. Legacy booleans are accepted (`True` = inverted,
            `False` = direct).

    Returns:
        List of position-aware repeat records.
    """
    repeat_type = normalizeRepeatType(repeat_type)
    query_string = window[:query_length]
    search_seq = window[query_length + min_spacer :]
    results: list[RepeatRecord] = []

    if imperfect_homology:
        errors = round(query_length * (1 - min_homology))
        if fixed_errors is not False:
            errors = int(fixed_errors)
    else:
        errors = 0
    min_candidate_length = max(1, query_length - errors)
    max_candidate_length = query_length + errors

    types_to_search = (
        ("direct", "inverted") if repeat_type == "both" else (repeat_type,)
    )
    query_by_type = {
        current_type: (
            rev_compl(query_string) if current_type == "inverted" else query_string
        )
        for current_type in types_to_search
    }
    patterns = {
        current_type: re.compile(
            f"({re.escape(query_by_type[current_type])}){{e<={errors}}}"
        )
        for current_type in types_to_search
        if imperfect_homology
    }

    candidate_starts = len(search_seq) - min_candidate_length + 1
    for i in range(candidate_starts):
        candidate_region = search_seq[i : i + max_candidate_length]

        for current_type in types_to_search:
            query = query_by_type[current_type]
            if imperfect_homology:
                candidate = findInteriorFuzzyMatch(
                    query,
                    candidate_region,
                    patterns[current_type],
                    min_candidate_length,
                    max_candidate_length,
                )
                is_match = candidate is not None
            else:
                candidate = candidate_region[:query_length]
                is_match = candidate == query

            if is_match:
                is_inverted = current_type == "inverted"
                match_start = window_start + query_length + min_spacer + i
                match_end = match_start + len(candidate)
                results.append(
                    {
                        "query_seq": query_string,
                        "match_seq": candidate,
                        "query_start": window_start,
                        "query_end": window_start + query_length,
                        "match_start": match_start,
                        "match_end": match_end,
                        "is_perfect": candidate == query,
                        "length": query_length,
                        "repeat_type": current_type,
                        "inverted": is_inverted,
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
