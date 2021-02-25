import numpy as np

from alignment_item import AlignmentItem
from needleman_wunsch import needleman_wunsch, needleman_wunsch_score_only


def make_score_matrix(
    align_profiles, del_cost, ins_cost, match_cost, mismatch_cost
):
    n = len(align_profiles)
    scores = np.zeros((n, n))
    for i in range(n):
        scores[i, i] = np.nan
        for j in range(i + 1, n):
            scores[i, j] = scores[j, i] = needleman_wunsch_score_only(
                align_profiles[i],
                align_profiles[j],
                del_cost,
                ins_cost,
                match_cost,
                mismatch_cost,
            )
    return scores


def nanargmax2d(x):
    n = x.shape[0]
    idx = np.nanargmax(x)
    i = idx // n
    j = idx - i * n
    assert i < j
    return i, j


def multiple_sequence_alignment(
    seqs, del_cost, ins_cost, match_cost, mismatch_cost
):
    n = len(seqs)
    align_profiles = [
        AlignmentItem.from_seq(seq, i) for i, seq in enumerate(seqs)
    ]
    scores = make_score_matrix(
        align_profiles, del_cost, ins_cost, match_cost, mismatch_cost
    )
    for i in range(n - 2):
        i, j = nanargmax2d(scores)
        score, align_profiles[i] = needleman_wunsch(
            align_profiles[i],
            align_profiles[j],
            del_cost,
            ins_cost,
            match_cost,
            mismatch_cost,
        )
        scores[..., j] = scores[j, ...] = np.nan

        for j in range(n):
            if not np.isnan(scores[i, j]):
                scores[i, j] = scores[j, i] = needleman_wunsch_score_only(
                    align_profiles[i],
                    align_profiles[j],
                    del_cost,
                    ins_cost,
                    match_cost,
                    mismatch_cost,
                )

    i, j = nanargmax2d(scores)
    ans_score, ans_alignment = needleman_wunsch(
        align_profiles[i],
        align_profiles[j],
        del_cost,
        ins_cost,
        match_cost,
        mismatch_cost,
    )

    # reorder as in input
    ans_alignment.reorder_seqs()

    return ans_score, ans_alignment.seqs
