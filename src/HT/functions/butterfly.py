"""Functions for butterfly calculations."""

# --- external imports
from fractions import Fraction


def chern(pval, qval):
    r"""Compute Chern numbers using the Diophantine equation.

    For a rational flux, the energy gaps in the Hofstadter spectrum are characterized by the integers :math:`s` and :math:`t`, which are related by the following Diophantine equation,

    .. math::
        r = qs_r + pt_r, \;\;\; |t_r|\leq\frac{q}{2}, \;\;\; s_r,t_r\in\mathbb{Z},

    where :math:`r` denotes the :math:`r`-th energy gap and :math:`t_r=\sum_{i=0}^r C_i` is the cumulative Chern number or Hall conductivity. :cite:`DiColandrea22`

    Parameters
    ----------
    pval: int
        The numerator of the flux density.
    qval: int
        The denominator of the flux density.

    Returns
    -------
    Chern_list: list
        The list of Chern numbers for each band (length M).
    tr_list: list
        The list of cumulative Chern numbers for each band gap (length M-1).
    """

    nphi = Fraction(pval, qval)
    p = nphi.numerator
    q = nphi.denominator

    # determine r and s
    sr_list, tr_list = [], []

    for r in range(q+1):
        if q % 2 == 0 and r == q/2:
            continue
        for tr in range(-int(q/2), int(q/2)+1):
            for sr in range(-q, q+1):
                if r == q*sr + p*tr:
                    sr_list.append(sr)
                    tr_list.append(tr)
                    break
            else:
                continue  # only executed if the inner loop did NOT break
            break  # only executed if the inner loop DID break

    Chern_list = []
    if q % 2 != 0:
        numb_band_groups = q
    else:
        numb_band_groups = q-1

    for i in range(numb_band_groups):
        Chern_list.append(tr_list[i+1] - tr_list[i])

    if q % 2 == 0:
        Chern_list.insert(q//2-1, Chern_list[q//2-1])

    return Chern_list, tr_list
