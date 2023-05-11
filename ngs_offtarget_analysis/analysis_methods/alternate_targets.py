from utils.constants import TSD, GENOME_BUFFER


def get_likely_target(
    reads_site, genome, spacer, read_is_fwd_strand, target_distance=40
):
    # this function find the most likely target sequence for a given read
    # looking at a 3bp window within 40bp upstream of the tLR integration site
    # using similarity to spacer as criteria, particularly in seed region
    read_position = reads_site + GENOME_BUFFER
    SPACER_LEN = len(spacer)

    debug = False

    best_match = None
    # evaluate targets
    for i in range(-3, 4, 1):
        if read_is_fwd_strand:
            # The read coming from the left side was forward strand, so
            # tLR implies the target is downstream of the read, with the
            # reverse complement
            test_site = read_position + target_distance - i
            test_seq = genome[test_site : test_site + SPACER_LEN].reverse_complement()
        else:
            # The read coming from the left side was reverse strand, so
            # tLR implies the target is upstream of the read
            test_site = read_position - target_distance - SPACER_LEN + i
            test_seq = genome[test_site : test_site + SPACER_LEN]

        best_match = evaluate_best_match(
            best_match, test_seq, test_site, True, spacer, debug=debug
        )
    return best_match


def evaluate_best_match(
    existing_best,
    tgt_seq,
    position,
    is_fwd_strand,
    spacer,
    without_5ths=False,
    debug=False,
):
    effective_tgt = tgt_seq
    effective_spacer = spacer
    if without_5ths:
        effective_spacer = "".join(
            [s for (idx, s) in enumerate(spacer) if ((idx + 1) % 6 != 0)]
        )
        effective_tgt = "".join(
            [s for (idx, s) in enumerate(tgt_seq) if ((idx + 1) % 6 != 0)]
        )

    matches = list(
        tgt_bp == spacer_bp
        for tgt_bp, spacer_bp in zip(effective_tgt, effective_spacer)
    )
    group_size = 7 if without_5ths else 8

    seed_match = sum(matches[:group_size])
    total_match = sum(matches)
    distal_match = sum(matches[-group_size:])
    best_match = existing_best

    if (
        not best_match
        or best_match["seed_match"] < seed_match
        or (
            best_match["seed_match"] == seed_match
            and best_match["total_match"] < total_match
        )
    ):
        best_match = {
            "seed_match": seed_match,
            "total_match": total_match,
            "distal_match": distal_match,
            "seq": tgt_seq,
            "pos": position,
            "strand": "fwd" if is_fwd_strand else "rv",
        }
    return best_match
