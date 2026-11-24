"""
Checks and corrections applied to a Chimera parcellation before it is saved.

In a Chimera parcellation the cortical white matter region ``wm-<hemi>-X`` has
the code of the cortical region ``ctx-<hemi>-X`` plus 3000. The codes are
derived here from the names of the final parcellation, so they do not depend on
counts of regions made while the parcellation is being assembled.
"""

from collections import Counter

import numpy as np

WM_OFFSET = 3000
HEMIS = ("lh", "rh")


def wm_code_map(
    names: list[str],
    codes: list[int],
    offset: int = WM_OFFSET,
    hemis: tuple[str, ...] = HEMIS,
) -> dict[int, int]:
    """
    Codes that the cortical white matter regions should have.

    Parameters
    ----------
    names : list of str
        Region names of the parcellation.

    codes : list of int
        Region codes, in the same order as ``names``.

    offset : int, optional
        Difference between the code of ``wm-<hemi>-X`` and the code of
        ``ctx-<hemi>-X``. Default is 3000.

    hemis : tuple of str, optional
        Hemispheres to consider. Default is ("lh", "rh").

    Returns
    -------
    dict
        ``{current code: correct code}`` for the white matter regions whose code
        is not the code of their cortical region plus ``offset``. White matter
        regions without a cortical region of the same name are not included.
    """
    code_of = dict(zip(names, (int(c) for c in codes), strict=True))
    remap = {}
    for hemi in hemis:
        wm_tag, ctx_tag = f"wm-{hemi}-", f"ctx-{hemi}-"
        for name, code in code_of.items():
            if not name.startswith(wm_tag):
                continue
            ctx_code = code_of.get(ctx_tag + name[len(wm_tag) :])
            if ctx_code is not None and code != ctx_code + offset:
                remap[code] = ctx_code + offset
    return remap


def fix_wm_codes(parc, offset: int = WM_OFFSET) -> dict[int, int]:
    """
    Give the cortical white matter regions the code of their cortical region plus ``offset``.

    The voxels are relabeled in one pass from the original codes, so codes that
    are shifted onto each other (the right code of one region being the wrong
    code of another) do not cascade. Names, colors and opacities are kept.

    Parameters
    ----------
    parc : clabtoolkit.parcellationtools.Parcellation
        Parcellation to correct in place.

    offset : int, optional
        Difference between the white matter and cortical codes. Default is 3000.

    Returns
    -------
    dict
        The applied ``{old code: new code}`` mapping (empty if nothing changed).

    Raises
    ------
    ValueError
        If a corrected code is already used by a region that keeps its code.
        The parcellation is not changed in that case.
    """
    codes = [int(c) for c in parc.index]
    remap = wm_code_map(list(parc.name), codes, offset=offset)
    if not remap:
        return remap

    kept = set(codes) - set(remap)
    clash = sorted(set(remap.values()) & kept)
    if clash:
        raise ValueError(
            f"The corrected white matter codes {clash[:10]} are already used by "
            "other regions."
        )

    keys = np.array(sorted(remap))
    vals = np.array([remap[k] for k in keys], dtype=parc.data.dtype)
    mask = np.isin(parc.data, keys)
    parc.data[mask] = vals[np.searchsorted(keys, parc.data[mask])]
    parc.index = [remap.get(c, c) for c in codes]

    # Keep the table sorted by code, as the rest of the pipeline expects
    order = np.argsort(parc.index, kind="stable")
    for attr in ("index", "name", "color", "opacity"):
        values = getattr(parc, attr, None)
        if values is not None and len(values) == len(order):
            setattr(parc, attr, [values[i] for i in order])
    if hasattr(parc, "parc_range"):
        parc.parc_range()
    return remap


def check_parcellation(parc, offset: int = WM_OFFSET) -> list[str]:
    """
    Look for inconsistencies between a parcellation and its color table.

    Parameters
    ----------
    parc : clabtoolkit.parcellationtools.Parcellation
        Parcellation to check.

    offset : int, optional
        Expected difference between the white matter and cortical codes.
        Default is 3000.

    Returns
    -------
    list of str
        One message per problem found; empty when the parcellation is consistent.
        The checks are: the table lists have the same length, no name or code
        is repeated, every code of the volume is in the table, and every
        ``wm-<hemi>-X`` region has a ``ctx-<hemi>-X`` region and its code plus
        ``offset``.
    """
    problems = []
    names, codes = list(parc.name), [int(c) for c in parc.index]

    lengths = {
        attr: len(getattr(parc, attr))
        for attr in ("index", "name", "color", "opacity")
        if getattr(parc, attr, None) is not None
    }
    if len(set(lengths.values())) > 1:
        problems.append(f"The color table lists have different lengths: {lengths}")
        return problems

    repeated_names = sorted(n for n, k in Counter(names).items() if k > 1)
    if repeated_names:
        problems.append(f"Repeated region names: {repeated_names[:10]}")
    repeated_codes = sorted(c for c, k in Counter(codes).items() if k > 1)
    if repeated_codes:
        problems.append(f"Repeated region codes: {repeated_codes[:10]}")

    in_volume = np.unique(parc.data)
    missing = sorted(set(in_volume[in_volume != 0].tolist()) - set(codes))
    if missing:
        problems.append(
            f"Codes in the volume without a color table entry: {missing[:10]}"
        )

    code_of = dict(zip(names, codes, strict=True))
    for hemi in HEMIS:
        wm_tag, ctx_tag = f"wm-{hemi}-", f"ctx-{hemi}-"
        wm_names = [n for n in names if n.startswith(wm_tag)]
        no_ctx = [n for n in wm_names if ctx_tag + n[len(wm_tag) :] not in code_of]
        if no_ctx:
            problems.append(
                f"{len(no_ctx)} {wm_tag}X region(s) without a {ctx_tag}X region: "
                f"{no_ctx[:10]}"
            )
        wrong = [
            f"{n} ({code_of[n]} instead of {code_of[ctx_tag + n[len(wm_tag):]] + offset})"
            for n in wm_names
            if n not in no_ctx
            and code_of[n] != code_of[ctx_tag + n[len(wm_tag) :]] + offset
        ]
        if wrong:
            problems.append(
                f"{len(wrong)} {wm_tag}X region(s) without the code of {ctx_tag}X + "
                f"{offset}: {wrong[:10]}"
            )
    return problems


def validate_parcellation(parc, offset: int = WM_OFFSET, label: str = "") -> None:
    """
    Raise a ValueError listing every problem found by ``check_parcellation``.

    Parameters
    ----------
    parc : clabtoolkit.parcellationtools.Parcellation
        Parcellation to check.

    offset : int, optional
        Expected difference between the white matter and cortical codes.
        Default is 3000.

    label : str, optional
        Name of the parcellation (e.g. its file name), used in the message.

    Raises
    ------
    ValueError
        If the parcellation is not consistent with its color table.
    """
    problems = check_parcellation(parc, offset=offset)
    if problems:
        where = f" {label}" if label else ""
        raise ValueError(
            f"The parcellation{where} is not consistent and was not saved:\n  - "
            + "\n  - ".join(problems)
        )
