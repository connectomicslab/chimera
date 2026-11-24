"""Tests of the checks and corrections applied before saving a Chimera parcellation."""

import numpy as np
import pytest
from clabtoolkit.parcellationtools import Parcellation

from chimera.validation import (
    check_parcellation,
    fix_wm_codes,
    validate_parcellation,
    wm_code_map,
)


def _parc(regions: dict, row: int = 0, shape=(16, 8, 4)) -> Parcellation:
    """
    Parcellation with the given {name: code} regions. Region k fills voxels of
    data[k, row, :], so parcellations built with different rows do not overlap.
    """
    data = np.zeros(shape, dtype=np.int32)
    for k, code in enumerate(regions.values()):
        data[k, row, : 1 + k % 3] = code
    n = len(regions)
    return Parcellation(
        data,
        color_table={
            "index": list(regions.values()),
            "name": list(regions.keys()),
            "color": [f"#00{k:02x}00" for k in range(n)],
            "opacity": [1.0] * n,
        },
    )


def _codes(parc) -> dict:
    return dict(zip(parc.name, parc.index, strict=True))


CORRECT = {
    "ctx-rh-a": 1,
    "ctx-rh-b": 2,
    "subc-rh-x": 3,
    "ctx-lh-a": 4,
    "ctx-lh-b": 5,
    "subc-lh-x": 6,
    "wm-brain-whitematter": 3000,
    "wm-rh-a": 3001,
    "wm-rh-b": 3002,
    "wm-lh-a": 3004,
    "wm-lh-b": 3005,
}


def test_correct_parcellation_has_no_problems():
    parc = _parc(CORRECT)
    assert check_parcellation(parc) == []
    assert wm_code_map(parc.name, parc.index) == {}
    assert fix_wm_codes(parc) == {}


def test_fix_wm_codes_shifted_left_hemisphere():
    # The left white matter codes are shifted by one: the right code of
    # wm-lh-a (3004) is free, and the right code of wm-lh-b (3005) is the
    # wrong code of wm-lh-a, so the relabeling must not cascade.
    shifted = dict(CORRECT, **{"wm-lh-a": 3005, "wm-lh-b": 3006})
    parc = _parc(shifted)
    n_a, n_b = (parc.data == 3005).sum(), (parc.data == 3006).sum()
    color_a = parc.color[parc.name.index("wm-lh-a")]

    problems = check_parcellation(parc)
    assert len(problems) == 1 and "wm-lh-" in problems[0]

    assert fix_wm_codes(parc) == {3005: 3004, 3006: 3005}
    assert _codes(parc) == CORRECT
    assert (parc.data == 3004).sum() == n_a and (parc.data == 3005).sum() == n_b
    assert parc.color[parc.name.index("wm-lh-a")] == color_a
    assert parc.index == sorted(parc.index)
    assert check_parcellation(parc) == []


def test_assembly_with_missing_right_regions():
    """
    Replay the assembly of ParcellationBuilder._assemble_with_cortex when the
    right non-cortical table lists more regions than the volume contains: the
    provisional left white matter codes (3000 + nrh_ctx + nrh_subc) are then
    too high, as in the parcellations of a cohort that had to be corrected.
    """
    rh_ctx = _parc({"ctx-rh-a": 1, "ctx-rh-b": 2}, row=0)
    rh_parc = _parc({"subc-rh-x": 1, "subc-rh-y": 2}, row=1)
    lh_ctx = _parc({"ctx-lh-a": 1, "ctx-lh-b": 2}, row=2)
    rh_wm = _parc({"wm-rh-a": 1, "wm-rh-b": 2}, row=3)
    lh_wm = _parc({"wm-lh-a": 1, "wm-lh-b": 2}, row=4)

    # Two more entries in the right table than regions in its volume
    rh_parc.index = rh_parc.index + [3, 4]
    rh_parc.name = rh_parc.name + ["subc-rh-gone1", "subc-rh-gone2"]
    rh_parc.color = rh_parc.color + ["#ffffff", "#ffffff"]
    rh_parc.opacity = rh_parc.opacity + [1.0, 1.0]
    nrh_ctx, nrh_subc = len(rh_ctx.index), len(rh_parc.index)

    # Same steps and offsets as the builder
    rh_wm.rearrange(offset=3000)
    lh_wm.rearrange(offset=3000 + nrh_ctx + nrh_subc)
    chim = Parcellation(np.zeros((16, 8, 4), dtype=np.int32))
    chim.add_parcellation(rh_ctx, append=True)
    chim.add_parcellation(rh_parc, append=True)
    chim.add_parcellation(lh_ctx, append=True)
    chim.add_parcellation(rh_wm, append=False)
    chim.add_parcellation(lh_wm, append=False)

    codes = _codes(chim)
    # The bug: ctx-lh-a is 5 but wm-lh-a is 3007 instead of 3005
    assert codes["ctx-lh-a"] == 5 and codes["wm-lh-a"] == 3007
    with pytest.raises(ValueError, match="wm-lh-"):
        validate_parcellation(chim, label="test.nii.gz")

    fix_wm_codes(chim)
    codes = _codes(chim)
    assert codes["wm-lh-a"] == codes["ctx-lh-a"] + 3000
    assert codes["wm-lh-b"] == codes["ctx-lh-b"] + 3000
    assert codes["wm-rh-a"] == 3001 and codes["wm-rh-b"] == 3002
    validate_parcellation(chim)


def test_fix_wm_codes_refuses_a_taken_code():
    # The right code of wm-lh-a (3004) is used by a region that keeps its code
    parc = _parc(dict(CORRECT, **{"wm-lh-a": 3010, "other": 3004, "wm-lh-b": 3005}))
    before = parc.data.copy()
    with pytest.raises(ValueError, match="already used"):
        fix_wm_codes(parc)
    assert np.array_equal(parc.data, before)


def test_check_parcellation_finds_problems():
    parc = _parc(dict(CORRECT, **{"wm-lh-z": 3099}))
    parc.data[15, 7, 3] = 77  # Code without a color table entry
    parc.name[parc.name.index("subc-lh-x")] = "subc-rh-x"  # Repeated name
    problems = "\n".join(check_parcellation(parc))
    assert "without a color table entry: [77]" in problems
    assert "Repeated region names: ['subc-rh-x']" in problems
    assert "wm-lh-X region(s) without a ctx-lh-X region: ['wm-lh-z']" in problems

    parc.color = parc.color[:-1]
    assert "different lengths" in check_parcellation(parc)[0]


def test_validate_parcellation_message():
    parc = _parc(dict(CORRECT, **{"wm-rh-a": 3050}))
    with pytest.raises(ValueError, match=r"sub-01_dseg\.nii\.gz is not consistent"):
        validate_parcellation(parc, label="sub-01_dseg.nii.gz")
