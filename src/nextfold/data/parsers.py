"""Functions for parsing various file formats."""

import collections
import dataclasses
import string
from typing import Iterable, Sequence

DeletionMatrix = Sequence[Sequence[int]]


@dataclasses.dataclass(frozen=True)
class TemplateHit:
    """Class representing a template hit."""

    index: int
    name: str
    aligned_cols: int
    sum_probs: float
    query: str
    hit_sequence: str
    indices_query: list[int]
    indices_hit: list[int]


def parse_fasta(
    fasta_string: str,
) -> tuple[Sequence[str], Sequence[str]]:
    """Parses FASTA string and returns list of strings with amino-acid
    sequences.

    Args:
        fasta_string (str): The string contents of a FASTA file.

    Returns:
        tuple[Sequence[str], Sequence[str]]:
        A tuple of:
            - A list of sequence.
            - A list of descriptions.

    """
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith(">"):
            index += 1
            descriptions.append(line[1:])
            sequences.append("")
            continue
        elif line.startswith("#"):
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line
    return sequences, descriptions


def parse_a3m(
    a3m_string: str,
) -> tuple[Sequence[str], DeletionMatrix, Sequence[str]]:
    """Parses sequences and deletion matrix from a3m format alignment.

    Args:
        a3m_string (str):
            The string contents of a a3m file.
            The first sequence in the file should be the query sequence.

    Returns:
        tuple[Sequence[str], DeletionMatrix, Sequence[str]]:
        A tuple of:
            - A list of sequences that have been aligned to the query.
                These might contain duplicates.
            - The deletion matrix for the alignment as a list of lists.
                The element at `deletion_matrix[i][j]` is the number of
                residues deleted from the aligned sequence `i` at residue
                position `j`.
            - A list of descriptions.
    """
    sequences, descriptions = parse_fasta(a3m_string)
    deletion_matrix = []
    for msa_sequence in sequences:
        deletion_vec = []
        deletion_count = 0
        for j in msa_sequence:
            if j.islower():
                deletion_count += 1
            else:
                deletion_vec.append(deletion_count)
                deletion_count = 0
        deletion_matrix.append(deletion_vec)

    # Make the MSA matrix out of aligned (deletion-free) sequences.
    deletion_table = str.maketrans("", "", string.ascii_lowercase)
    aligned_sequences = [s.translate(deletion_table) for s in sequences]
    return aligned_sequences, deletion_matrix, descriptions


def parse_stockholm(
    stockholm_string: str,
) -> tuple[Sequence[str], DeletionMatrix, Sequence[str]]:
    """Parses sequences and deletion matrix from stockholm format alignment.

    Args:
        stockholm_string (str):
            The string contents of a stockholm file.
            The first sequence in the file should be the query sequence.


    Returns:
        tuple[Sequence[str], DeletionMatrix, Sequence[str]]:
        A tuple of:
            - A list of sequences that have been aligned to the query.
                These might contain duplicates.
            - The deletion matrix for the alignment as a list of lists.
                The element at `deletion_matrix[i][j]` is the number of
                residues deleted from the aligned sequence `i` at residue
                position `j`.
            - The names of the targets matched, including the jackhmmer
                subsequence suffix.

    """
    name_to_sequence = collections.OrderedDict()
    for line in stockholm_string.splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "//")):
            continue
        name, sequence = line.split()
        if name not in name_to_sequence:
            name_to_sequence[name] = ""
        name_to_sequence[name] += sequence

    msa = []
    deletion_matrix = []

    query = ""
    keep_columns = []
    for seq_index, sequence in enumerate(name_to_sequence.values()):
        if seq_index == 0:
            # Gather the columns with gaps from the query
            query = sequence
            keep_columns = [i for i, res in enumerate(query) if res != "-"]

        # Remove the columns with gaps in the query from all sequences.
        aligned_sequence = "".join([sequence[c] for c in keep_columns])

        msa.append(aligned_sequence)

        # Count the number of deletions w.r.t. query.
        deletion_vec = []
        deletion_count = 0
        for seq_res, query_res in zip(sequence, query):
            if seq_res != "-" or query_res != "-":
                if query_res == "-":
                    deletion_count += 1
                else:
                    deletion_vec.append(deletion_count)
                    deletion_count = 0
        deletion_matrix.append(deletion_vec)

    return msa, deletion_matrix, list(name_to_sequence.keys())


def _convert_sto_seq_to_a3m(
    query_non_gaps: Sequence[bool],
    sto_seq: str,
) -> Iterable[str]:
    for is_query_res_non_gap, sequence_res in zip(query_non_gaps, sto_seq):
        if is_query_res_non_gap:
            yield sequence_res
        elif sequence_res != "-":
            yield sequence_res.lower()


def convert_stockholm_to_a3m(
    stockholm_format: str,
    max_sequences: int | None = None,
) -> str:
    """Converts MSA in Stockholm format to the A3M format."""
    descriptions = {}
    sequences = {}
    reached_max_sequences = False

    for line in stockholm_format.splitlines():
        reached_max_sequences = max_sequences and len(sequences) >= max_sequences
        if line.strip() and not line.startswith(("#", "//")):
            # Ignore blank lines, markup and end symbols.
            # Remainder are alignment¸sequence parts.
            seqname, aligned_seq = line.split(maxsplit=1)
            if seqname not in sequences:
                if reached_max_sequences:
                    continue
                sequences[seqname] = ""
            sequences[seqname] += aligned_seq

    for line in stockholm_format.splitlines():
        if line[:4] == "#=GS":
            columns = line.split(maxsplit=3)
            seqname, feature = columns[1:3]
            value = columns[3] if len(columns) == 4 else ""
            if feature != "DE":
                continue
            if reached_max_sequences and seqname not in sequences:
                continue
            descriptions[seqname] = value
            if len(descriptions) == len(sequences):
                break

    # Convert sto format to a3m line by line
    a3m_sequences = {}
    # `query_sequence` is assumed to be the first sequence
    query_sequence = next(iter(sequences.values()))
    query_non_gaps = [res != "-" for res in query_sequence]
    for seqname, sto_sequence in sequences.items():
        a3m_sequences[seqname] = "".join(
            _convert_sto_seq_to_a3m(query_non_gaps, sto_sequence)
        )

    fasta_chunks = (
        f">{k} {descriptions.get(k, '')}\n{a3m_sequences[k]}" for k in a3m_sequences
    )
    return "\n".join(fasta_chunks) + "\n"  # Include terminating newline.
