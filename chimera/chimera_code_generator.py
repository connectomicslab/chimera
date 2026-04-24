#!/usr/bin/env python3
"""
chimera_code_generator.py
=========================
Interactive CLI tool to build a CHIMERA parcellation code from a JSON
configuration file.  For each supra-region the user is shown all available
atlases and prompted to select one by its single-letter key.

For atlases that expose multiple parcel variants (via ``seg`` and/or ``scale``
BIDS-like entities in the parcel name), a follow-up prompt allows the user to
select one or more values per dimension (comma-separated).  The final parcel
list is the cross-product of the selected seg(s) and scale(s).

The session ends with the composed code (e.g. "HFIIIISIFN") and a
ready-to-paste command line printed to the terminal.

The interactive workflow prompts for:
  1. One atlas per supra-region (letter key, Enter = default).
  2. Seg/scale variant selection for multi-resolution atlases.
  3. White-matter growing distance(s): non-negative integer(s) or ``wm``,
     comma-separated (e.g. ``0``, ``0,1``, ``0,wm``, ``wm``).  Default: ``0``.

Usage
-----
    python chimera_code_generator.py
    python chimera_code_generator.py parcellations.json
    python chimera_code_generator.py parcellations.json --non-interactive HFIIIISIFN
"""

import argparse
import json
import re
import os
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Terminal colours
# ---------------------------------------------------------------------------

from clabtoolkit.colorstools import bcolors


def _c(color: str, text: str) -> str:
    """Wrap *text* in a colour escape and reset."""
    return f"{color}{text}{bcolors.ENDC}"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_parcellations_info(parc_dict_file: str = None) -> dict:
    """
    Load the parcellation configuration from a JSON file.

    Parameters
    ----------
    parc_dict_file : str
        Path to the JSON configuration file.

    Returns
    -------
    dict
        Parsed JSON content.

    Raises
    ------
    SystemExit
        If the file cannot be found or parsed.
    """

    # Detecting the base directory
    chim_dir = os.path.dirname(os.path.abspath(__file__))

    if parc_dict_file is not None:
        if not os.path.isfile(parc_dict_file):
            raise ValueError("The parcellation dictionary file does not exist")
    else:
        parc_dict_file = os.path.join(
            chim_dir, "config", "supraregions_dictionary.json"
        )

    path = Path(parc_dict_file)
    if not path.is_file():
        print(_c(bcolors.OKRED, f"Error: file not found -> {parc_dict_file}"))
        sys.exit(1)
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except json.JSONDecodeError as exc:
        print(_c(bcolors.OKRED, f"Error: could not parse JSON -> {exc}"))
        sys.exit(1)


# ---------------------------------------------------------------------------
# Parcel variant helpers
# ---------------------------------------------------------------------------


def _parse_parcel_variants(parcels: list) -> dict:
    """
    Scan a list of parcel strings for BIDS-like ``seg`` and ``scale`` entities
    and collect the unique values for each dimension in their original order.

    Examples
    --------
    >>> _parse_parcel_variants(["atlas-HCPMM1_seg-5p3", "atlas-HCPMM1_seg-7p1"])
    {"seg": ["5p3", "7p1"], "scale": []}

    >>> _parse_parcel_variants(["atlas-Schaefer2018_seg-7n_scale-100",
    ...                          "atlas-Schaefer2018_seg-17n_scale-100"])
    {"seg": ["7n", "17n"], "scale": ["100"]}

    Parameters
    ----------
    parcels : list of str
        Raw parcel name strings from the JSON.

    Returns
    -------
    dict
        Keys ``"seg"`` and ``"scale"``, each a list of unique string values
        (empty list when the dimension is absent).
    """
    segs, scales = [], []
    for p in parcels:
        m_seg = re.search(r"seg-([^_]+)", p)
        m_scale = re.search(r"scale-([^_]+)", p)
        if m_seg and m_seg.group(1) not in segs:
            segs.append(m_seg.group(1))
        if m_scale and m_scale.group(1) not in scales:
            scales.append(m_scale.group(1))
    return {"seg": segs, "scale": scales}


def _has_variants(parcels) -> bool:
    """
    Return True when *parcels* is a list with more than one entry AND at least
    one entry contains a ``seg`` or ``scale`` entity.

    Parameters
    ----------
    parcels : str or list
        Value of the ``parcels`` field from the JSON.

    Returns
    -------
    bool
    """
    if not isinstance(parcels, list) or len(parcels) <= 1:
        return False
    v = _parse_parcel_variants(parcels)
    return bool(v["seg"] or v["scale"])


def _variant_hint(parcels: list) -> str:
    """
    Build a compact one-line hint summarising the parcel variants, e.g.:
    ``"3 variants  ->  seg: 5p3, 6p0, 7p1"``
    ``"30 variants  ->  seg: 7n, 17n, kong17n  x  scale: 100 ... 1000"``

    Parameters
    ----------
    parcels : list of str

    Returns
    -------
    str
    """
    variants = _parse_parcel_variants(parcels)
    parts = []
    if variants["seg"]:
        parts.append("seg: " + ", ".join(variants["seg"]))
    if variants["scale"]:
        sc = variants["scale"]
        scale_str = f"{sc[0]} ... {sc[-1]}" if len(sc) > 4 else ", ".join(sc)
        parts.append("scale: " + scale_str)
    dims = "  x  ".join(parts)
    return f"{len(parcels)} variants  ->  {dims}"


def _filter_parcels(
    parcels: list,
    selected_segs: list,
    selected_scales: list,
    has_seg: bool,
    has_scale: bool,
) -> list:
    """
    Return the subset of *parcels* whose ``seg`` and ``scale`` values are in
    the respective selection lists.

    Parameters
    ----------
    parcels : list of str
        Full parcel list for the atlas.
    selected_segs : list of str
        User-chosen seg values (empty when dimension is absent).
    selected_scales : list of str
        User-chosen scale values (empty when dimension is absent).
    has_seg : bool
        Whether the atlas uses the seg dimension at all.
    has_scale : bool
        Whether the atlas uses the scale dimension at all.

    Returns
    -------
    list of str
        Filtered parcel names.
    """
    result = []
    for p in parcels:
        m_seg = re.search(r"seg-([^_]+)", p)
        m_scale = re.search(r"scale-([^_]+)", p)
        seg_ok = (not has_seg) or (m_seg and m_seg.group(1) in selected_segs)
        scale_ok = (not has_scale) or (m_scale and m_scale.group(1) in selected_scales)
        if seg_ok and scale_ok:
            result.append(p)
    return result


# ---------------------------------------------------------------------------
# Default resolution
# ---------------------------------------------------------------------------


def _get_default_key(options: dict) -> str | None:
    """
    Return the single-letter key whose ``"default"`` field is ``True``, or
    ``None`` if no such entry exists.

    Parameters
    ----------
    options : dict
        Mapping of single-letter key -> atlas metadata dict for one region.

    Returns
    -------
    str or None
        The key marked as default, or None.
    """
    for key, meta in options.items():
        if meta.get("default") is True:
            return key
    return None


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

NONE_KEY = "N"
NONE_LABEL = "None  (skip this region)"


def _print_region_options(
    region_name: str, options: dict, default_key: str | None
) -> None:
    """
    Print all available atlas options for a single supra-region.
    The atlas marked as default is highlighted with a ``(default)`` tag.
    For atlases with multiple parcel variants a compact hint line is shown.

    Parameters
    ----------
    region_name : str
        Display name of the supra-region (e.g. "Thalamus").
    options : dict
        Mapping of single-letter key -> atlas metadata dict.
    default_key : str or None
        The key whose ``"default"`` field is True, or None.
    """
    print()
    print(
        _c(bcolors.BOLD, _c(bcolors.DARKCYAN, f"  {region_name}"))
        + _c(
            bcolors.DIM, f"  ({len(options)} option{'s' if len(options) != 1 else ''})"
        )
    )
    print(_c(bcolors.DIM, "  " + "-" * 62))

    for key, meta in options.items():
        name = meta.get("name", "")
        cita = meta.get("citation", "")
        src = meta.get("source", "")
        atype = meta.get("type", "")
        ref = meta.get("reference", "")
        parcels = meta.get("parcels", [])
        is_def = key == default_key

        default_tag = "  " + _c(bcolors.OKGREEN, "(default)") if is_def else ""
        print(
            "  "
            + _c(bcolors.OKGREEN + bcolors.BOLD, f"[{key}]")
            + "  "
            + _c(bcolors.BOLD, f"{name:<28}")
            + _c(bcolors.OKYELLOW, cita)
            + default_tag
        )
        print(
            "       "
            + _c(bcolors.DIM, f"source: {src:<14}  type: {atype:<10}  ref: {ref}")
        )
        if _has_variants(parcels):
            print("       " + _c(bcolors.OKBLUE, _variant_hint(parcels)))

    print(
        "  "
        + _c(bcolors.OKMAGENTA + bcolors.BOLD, f"[{NONE_KEY}]")
        + "  "
        + _c(bcolors.DIM, NONE_LABEL)
    )
    print()


def _print_variant_options(dim_name: str, values: list) -> None:
    """
    Print the available values for a single variant dimension (seg or scale).

    Parameters
    ----------
    dim_name : str
        Dimension label, e.g. ``"seg"`` or ``"scale"``.
    values : list of str
        Available values in their original order.
    """
    formatted = "  ".join(_c(bcolors.OKGREEN + bcolors.BOLD, f"[{v}]") for v in values)
    print(
        "    "
        + _c(bcolors.DARKCYAN + bcolors.BOLD, f"{dim_name:<8}")
        + "  "
        + formatted
    )


def _print_summary(data: dict, code: dict) -> None:
    """
    Print a formatted summary table of all chosen atlases and, where
    applicable, the selected parcel variants.

    Parameters
    ----------
    data : dict
        Full parcellation JSON.
    code : dict
        Mapping of region -> {"key": str, "selected_parcels": list | None}.
    """
    print()
    print(_c(bcolors.BOLD, "  -- Selection summary " + "-" * 40))

    for region, choice in code.items():
        key = choice["key"]
        sel_parcels = choice.get("selected_parcels")

        if key == NONE_KEY:
            label = _c(bcolors.DIM, "skipped")
            cita = ""
        else:
            meta = data[region][key]
            label = _c(bcolors.OKBLUE, meta.get("name", key))
            cita = _c(bcolors.OKYELLOW, meta.get("citation", ""))

        print(
            f"  {_c(bcolors.DARKCYAN + bcolors.BOLD, f'{region:<16}')}"
            f"  {_c(bcolors.OKGREEN + bcolors.BOLD, f'[{key}]')}  "
            f"{label}  {cita}"
        )

        if sel_parcels:
            for p in sel_parcels:
                print("       " + _c(bcolors.DIM, f"-> {p}"))

    print(_c(bcolors.DIM, "  " + "-" * 62))


# ---------------------------------------------------------------------------
# Core interaction
# ---------------------------------------------------------------------------


def _prompt_region(region_name: str, options: dict) -> str:
    """
    Display options for one region and return the validated single-letter key
    chosen by the user.  Pressing Enter selects the atlas marked as default
    in the JSON (``"default": true``).  Loops until a valid entry is given.

    Parameters
    ----------
    region_name : str
        Supra-region being configured.
    options : dict
        Available atlas keys -> metadata.

    Returns
    -------
    str
        Single uppercase letter key, or "N" for skip.
    """
    default_key = _get_default_key(options)
    valid_keys = set(options.keys()) | {NONE_KEY}
    _print_region_options(region_name, options, default_key)

    enter_hint = _c(bcolors.DIM, f"  (Enter = {default_key})") if default_key else ""

    while True:
        raw = (
            input(
                "  "
                + _c(bcolors.BOLD, f"Select parcellation for {region_name}")
                + _c(bcolors.DIM, f"  [{'/'.join(sorted(valid_keys))}]")
                + enter_hint
                + ": "
            )
            .strip()
            .upper()
        )

        if raw == "" and default_key:
            print("  " + _c(bcolors.DIM, f"-> {default_key}  (default)"))
            return default_key

        if raw in valid_keys:
            return raw

        print(
            _c(bcolors.OKRED, f"  x  '{raw}' is not available for {region_name}.")
            + _c(bcolors.DIM, f"  Valid options: {', '.join(sorted(valid_keys))}")
        )
        _print_region_options(region_name, options, default_key)


def _prompt_variant_dim(dim_name: str, values: list) -> list:
    """
    Prompt the user to select one or more values for a single parcel variant
    dimension.  Multiple values are entered comma-separated.  An empty entry
    selects all values.

    Parameters
    ----------
    dim_name : str
        ``"seg"`` or ``"scale"``.
    values : list of str
        All valid choices for this dimension.

    Returns
    -------
    list of str
        Validated subset of *values* chosen by the user, in original order.
    """
    valid = set(values)
    all_hint = _c(bcolors.DIM, "  (Enter = all)")

    while True:
        raw = input(
            "    "
            + _c(bcolors.BOLD, f"Select {dim_name}")
            + _c(bcolors.DIM, f"  [{', '.join(values)}]")
            + all_hint
            + ": "
        ).strip()

        if raw == "":
            print("    " + _c(bcolors.DIM, "-> all  (default)"))
            return list(values)

        chosen = [v.strip() for v in raw.split(",") if v.strip()]
        invalid = [v for v in chosen if v not in valid]

        if not invalid:
            return [v for v in values if v in set(chosen)]  # preserve order

        print(
            _c(bcolors.OKRED, f"    x  Unknown value(s): {', '.join(invalid)}.")
            + _c(bcolors.DIM, f"  Valid: {', '.join(values)}")
        )


def _prompt_parcel_selection(atlas_name: str, parcels: list) -> list:
    """
    Drive the sub-selection workflow for an atlas with multiple parcel variants.
    Prompts separately for ``seg`` and ``scale`` dimensions (only those present),
    then returns the matching subset of *parcels*.

    Parameters
    ----------
    atlas_name : str
        Human-readable atlas name shown in prompts.
    parcels : list of str
        Full parcel list for the atlas.

    Returns
    -------
    list of str
        Selected (filtered) parcel names.
    """
    variants = _parse_parcel_variants(parcels)
    has_seg = bool(variants["seg"])
    has_scale = bool(variants["scale"])

    print()
    print(
        "  "
        + _c(bcolors.BOLD, f"{atlas_name}")
        + _c(
            bcolors.DIM,
            f" has {len(parcels)} parcel variants -- select the ones you need:",
        )
    )
    print(_c(bcolors.DIM, "  " + "-" * 62))

    if has_seg:
        _print_variant_options("seg", variants["seg"])
    if has_scale:
        _print_variant_options("scale", variants["scale"])
    print()

    selected_segs = _prompt_variant_dim("seg", variants["seg"]) if has_seg else []
    selected_scales = (
        _prompt_variant_dim("scale", variants["scale"]) if has_scale else []
    )

    result = _filter_parcels(
        parcels, selected_segs, selected_scales, has_seg, has_scale
    )

    if not result:
        print(
            _c(
                bcolors.OKRED,
                "  x  No parcels match that combination -- selecting all.",
            )
        )
        result = list(parcels)

    print("\n  " + _c(bcolors.OKGREEN, f"+ {len(result)} parcel(s) selected:"))
    for p in result:
        print("     " + _c(bcolors.DIM, p))

    return result


def run_interactive(data: dict) -> dict:
    """
    Walk through every supra-region interactively and collect the user's
    choices.  For atlases with multi-variant parcel lists a follow-up prompt
    is shown to narrow down seg/scale variants.

    Parameters
    ----------
    data : dict
        Full parcellation JSON.

    Returns
    -------
    dict
        Ordered mapping of region ->
        {"key": str, "selected_parcels": list | None}.
    """
    print()
    print(_c(bcolors.BOLD + bcolors.OKBLUE, "  CHIMERA Parcellation Code Generator"))
    print(_c(bcolors.DIM, "  " + "=" * 62))
    print(
        _c(
            bcolors.DIM,
            "  For each supra-region select a parcellation by its letter key.\n"
            f"  Enter [{NONE_KEY}] to skip a region (it will appear as N in the code).",
        )
    )

    code = {}
    regions = list(data.keys())
    total = len(regions)

    for idx, region in enumerate(regions, start=1):
        print(_c(bcolors.DIM, f"\n  Step {idx}/{total}"))

        key = _prompt_region(region, data[region])

        selected_parcels = None
        if key != NONE_KEY:
            parcels = data[region][key].get("parcels", [])
            if _has_variants(parcels):
                atlas_name = data[region][key].get("name", key)
                selected_parcels = _prompt_parcel_selection(atlas_name, parcels)

        code[region] = {"key": key, "selected_parcels": selected_parcels}

    return code


def validate_code(data: dict, code_str: str) -> dict:
    """
    Parse and validate a pre-supplied code string (non-interactive mode).
    Parcel sub-selection is not performed in this mode.

    Parameters
    ----------
    data : dict
        Full parcellation JSON.
    code_str : str
        Code string whose length must match the number of supra-regions.

    Returns
    -------
    dict
        Mapping of region -> {"key": str, "selected_parcels": None}.

    Raises
    ------
    SystemExit
        On length mismatch or invalid key for any region.
    """
    regions = list(data.keys())
    code_str = code_str.strip().upper()

    if len(code_str) != len(regions):
        print(
            _c(
                bcolors.OKRED,
                f"Error: code length {len(code_str)} does not match "
                f"number of regions ({len(regions)}).",
            )
        )
        sys.exit(1)

    code = {}
    for region, key in zip(regions, code_str):
        valid = set(data[region].keys()) | {NONE_KEY}
        if key not in valid:
            print(
                _c(
                    bcolors.OKRED,
                    f"Error: key '{key}' is not valid for region '{region}'.\n"
                    f"       Valid keys: {', '.join(sorted(valid))}",
                )
            )
            sys.exit(1)
        code[region] = {"key": key, "selected_parcels": None}

    return code


# ---------------------------------------------------------------------------
# White-matter growing prompt
# ---------------------------------------------------------------------------


def _prompt_growwm() -> list:
    """
    Ask the user for white-matter growing distance(s).

    Accepts one or more comma-separated values; each token must be either a
    non-negative integer or the literal string ``"wm"``.  Pressing Enter
    selects the default value of ``["0"]`` (no growing).

    Returns
    -------
    list of str
        Validated growing values, e.g. ``["0"]``, ``["0", "1"]``,
        ``["wm"]``, or ``["0", "wm"]``.
    """
    print()
    print(_c(bcolors.BOLD + bcolors.DARKCYAN, "  White matter growing"))
    print(_c(bcolors.DIM, "  " + "-" * 62))
    print(
        _c(
            bcolors.DIM,
            "  Enter one or more comma-separated values.\n"
            "  Each value must be a non-negative integer or 'wm'.\n"
            "  Examples:  0   |   2   |   0,1   |   0,wm   |   wm",
        )
    )
    print()

    while True:
        raw = input(
            "  "
            + _c(bcolors.BOLD, "WM growing distance(s)")
            + _c(bcolors.DIM, "  [integer(s) or wm, comma-separated]")
            + _c(bcolors.DIM, "  (Enter = 0)")
            + ": "
        ).strip()

        if raw == "":
            print("  " + _c(bcolors.DIM, "-> 0  (default)"))
            return ["0"]

        tokens = [v.strip() for v in raw.split(",") if v.strip()]
        invalid = []
        result = []
        for t in tokens:
            if t.lower() == "wm":
                if "wm" not in result:
                    result.append("wm")
            else:
                try:
                    val = int(t)
                    if val < 0:
                        raise ValueError
                    s = str(val)
                    if s not in result:
                        result.append(s)
                except ValueError:
                    invalid.append(t)

        if invalid:
            print(
                _c(bcolors.OKRED, f"  x  Invalid value(s): {', '.join(invalid)}.")
                + _c(bcolors.DIM, "  Use non-negative integers or 'wm'.")
            )
            continue

        if not result:
            print(_c(bcolors.OKRED, "  x  Please enter at least one value."))
            continue

        return result


# ---------------------------------------------------------------------------
# Command-line summary
# ---------------------------------------------------------------------------


def _build_command_line(
    code_string: str, code: dict, growwm: list | None = None
) -> str:
    """
    Build a ready-to-paste command-line string summarising the full selection,
    e.g. ``--parcodes SFIIIISIFN --seg 7n,17n --scale 100,400 --growwm 0,1``.

    Seg and scale values are collected from all selected parcel names across
    every region, deduplicated, and kept in their original encounter order.

    Parameters
    ----------
    code_string : str
        The compact letter code (e.g. "SFIIIISIFN").
    code : dict
        Region -> {"key": str, "selected_parcels": list | None}.
    growwm : list of str, optional
        White-matter growing values (e.g. ``["0", "1"]``).  Omitted from
        the command line when ``None`` or equal to ``["0"]``.

    Returns
    -------
    str
        Space-joined argument string.
    """
    all_segs, all_scales = [], []
    for choice in code.values():
        for p in choice.get("selected_parcels") or []:
            m_seg = re.search(r"seg-([^_]+)", p)
            m_scale = re.search(r"scale-([^_]+)", p)
            if m_seg and m_seg.group(1) not in all_segs:
                all_segs.append(m_seg.group(1))
            if m_scale and m_scale.group(1) not in all_scales:
                all_scales.append(m_scale.group(1))

    parts = [f"--parcodes {code_string}"]
    if all_segs:
        parts.append(f"--seg {','.join(all_segs)}")
    if all_scales:
        parts.append(f"--scale {','.join(all_scales)}")
    parts.append(f"--growwm {','.join(growwm) if growwm else '0'}")
    return "  ".join(parts)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def build_code_string(code: dict) -> str:
    """Concatenate the single-letter keys into the final code string."""
    return "".join(choice["key"] for choice in code.values())


def save_code(code_string: str, output_path: str, code: dict, data: dict) -> None:
    """
    Write the compact code and a human-readable breakdown to a text file.

    Parameters
    ----------
    code_string : str
        The compact code (e.g. "HFIIIISIFN").
    output_path : str
        Destination file path.
    code : dict
        Region -> {"key": str, "selected_parcels": list | None}.
    data : dict
        Full parcellation JSON (used for labels).
    """
    lines = [
        "CHIMERA Parcellation Code",
        "=" * 50,
        f"Code: {code_string}",
        "",
        "Breakdown:",
        "-" * 50,
    ]
    for region, choice in code.items():
        key = choice["key"]
        sel_parcels = choice.get("selected_parcels")

        if key == NONE_KEY:
            lines.append(f"  {region:<16}  [{key}]  skipped")
        else:
            meta = data[region][key]
            lines.append(
                f"  {region:<16}  [{key}]  {meta.get('name', '')}  "
                f"{meta.get('citation', '')}"
            )
            if sel_parcels:
                for p in sel_parcels:
                    lines.append(f"                       ->  {p}")

    lines.append("")
    Path(output_path).write_text("\n".join(lines), encoding="utf-8")
    print(_c(bcolors.OKGREEN, f"  + Code saved to: {output_path}"))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="chimera_code_generator",
        description="Interactively build a CHIMERA parcellation code from a JSON config.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python chimera_code_generator.py\n"
            "  python chimera_code_generator.py parcellations.json\n"
            "  python chimera_code_generator.py parcellations.json --non-interactive HFIIIISIFN"
        ),
    )
    parser.add_argument(
        "json_file",
        metavar="JSON_FILE",
        nargs="?",
        default=None,
        help="Path to the parcellation configuration JSON file. "
        "Defaults to chimera's built-in supraregions_dictionary.json.",
    )
    parser.add_argument(
        "--non-interactive",
        "-n",
        metavar="CODE",
        default=None,
        dest="preset_code",
        help="Validate and display an existing code without prompting.",
    )
    args = parser.parse_args()

    data = load_parcellations_info(args.json_file)

    # -- Mode selection -------------------------------------------------------
    if args.preset_code:
        code = validate_code(data, args.preset_code)
        growwm = ["0"]
    else:
        try:
            code = run_interactive(data)
            growwm = _prompt_growwm()
        except (KeyboardInterrupt, EOFError):
            print(_c(bcolors.OKRED, "\n\n  Interrupted -- no code generated."))
            sys.exit(0)

    # -- Results --------------------------------------------------------------
    _print_summary(data, code)
    code_string = build_code_string(code)

    print(
        "  "
        + _c(bcolors.BOLD + bcolors.DARKCYAN, f"{'WM growing':<16}")
        + "  "
        + _c(bcolors.OKGREEN + bcolors.BOLD, f"[{', '.join(growwm)}]")
    )

    cmd_line = _build_command_line(code_string, code, growwm)
    print()
    print(
        "  "
        + _c(bcolors.BOLD, "Generated code:  ")
        + _c(bcolors.BOLD + bcolors.OKBLUE, f"  {code_string}  ")
    )
    print("  " + _c(bcolors.DIM, "Command line:    ") + _c(bcolors.DARKWHITE, cmd_line))
    print()


if __name__ == "__main__":
    main()
