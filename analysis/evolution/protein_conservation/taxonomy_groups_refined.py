#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Shared refined taxonomy configuration and normalization helpers."""

from __future__ import annotations

from copy import deepcopy
from typing import Optional

import pandas as pd


GROUP_ORDER = [
    "Fungi",
    "Invertebrates",
    "Basal_chordates",
    "Fish",
    "Amphibians",
    "Sauropsids",
    "Monotremes_Marsupials",
    "Mammals_non_primate",
    "Primates_Hominids",
]

GROUP_LABELS = {
    "Fungi": "Fungi",
    "Invertebrates": "Invertebrates",
    "Basal_chordates": "Basal chordates",
    "Fish": "Fish",
    "Amphibians": "Amphibians",
    "Sauropsids": "Sauropsids",
    "Monotremes_Marsupials": "Monotremes / Marsupials",
    "Mammals_non_primate": "Placental mammals",
    "Primates_Hominids": "Primates / Hominids",
    "Other": "Other",
}

GROUP_COLORS = {
    "Fungi": "#b07aa1",
    "Invertebrates": "#f28e2b",
    "Basal_chordates": "#76b7b2",
    "Fish": "#4e79a7",
    "Amphibians": "#7aa6c2",
    "Sauropsids": "#59a14f",
    "Monotremes_Marsupials": "#af7aa1",
    "Mammals_non_primate": "#9c755f",
    "Primates_Hominids": "#e15759",
    "Other": "#bab0ab",
}

FUNGI_GENERA = {"saccharomyces", "schizosaccharomyces", "neurospora", "aspergillus"}
INVERTEBRATE_GENERA = {
    "caenorhabditis",
    "pristionchus",
    "drosophila",
    "anopheles",
    "apis",
    "tribolium",
    "strongylocentrotus",
}
BASAL_CHORDATE_GENERA = {"ciona", "branchiostoma", "petromyzon", "eptatretus"}
EXPLICIT_FISH_GENERA = {"danio", "takifugu", "oryzias", "callorhinchus", "latimeria"}
AMPHIBIAN_GENERA = {"xenopus", "leptobrachium"}
SAUROPSID_GENERA = {
    "gallus",
    "meleagris",
    "taeniopygia",
    "anas",
    "anser",
    "aquila",
    "coturnix",
    "ficedula",
    "geospiza",
    "parus",
    "serinus",
    "strigops",
    "struthio",
    "anolis",
    "chelonia",
    "chelonoidis",
    "chrysemys",
    "crocodylus",
    "gopherus",
    "laticauda",
    "naja",
    "notechis",
    "pelodiscus",
    "podarcis",
    "pseudonaja",
    "salvator",
    "sphenodon",
    "terrapene",
}
PRIMATE_GENERA = {
    "aotus",
    "callithrix",
    "carlito",
    "cebus",
    "cercocebus",
    "chlorocebus",
    "gorilla",
    "homo",
    "macaca",
    "mandrillus",
    "microcebus",
    "nomascus",
    "otolemur",
    "pan",
    "papio",
    "pongo",
    "prolemur",
    "propithecus",
    "rhinopithecus",
    "saimiri",
}
MONOTREME_MARSUPIAL_GENERA = {
    "ornithorhynchus",
    "monodelphis",
    "sarcophilus",
    "notamacropus",
    "phascolarctos",
    "vombatus",
}
MAMMAL_GENERA = {
    "ailuropoda",
    "balaenoptera",
    "bison",
    "bos",
    "camelus",
    "canis",
    "capra",
    "catagonus",
    "cavia",
    "cervus",
    "chinchilla",
    "choloepus",
    "dasypus",
    "delphinapterus",
    "dipodomys",
    "echinops",
    "equus",
    "erinaceus",
    "felis",
    "heterocephalus",
    "ictidomys",
    "jaculus",
    "loxodonta",
    "marmota",
    "mesocricetus",
    "microtus",
    "monodon",
    "moschus",
    "mus",
    "mustela",
    "myotis",
    "nannospalax",
    "neovison",
    "ochotona",
    "octodon",
    "ovis",
    "oryctolagus",
    "panthera",
    "peromyscus",
    "phocoena",
    "physeter",
    "procavia",
    "pteropus",
    "rattus",
    "rhinolophus",
    "cricetulus",
    "sciurus",
    "sorex",
    "sus",
    "tupaia",
    "tursiops",
    "urocitellus",
    "ursus",
    "vicugna",
    "vulpes",
}

LEGACY_TO_REFINED = {
    "Fish_Amphibians": "Fish",
}

TREE_TEMPLATE = {
    "label": "root",
    "children": [
        {"group": "Fungi", "label": GROUP_LABELS["Fungi"]},
        {
            "label": "Metazoa",
            "children": [
                {"group": "Invertebrates", "label": GROUP_LABELS["Invertebrates"]},
                {
                    "label": "Chordates",
                    "children": [
                        {"group": "Basal_chordates", "label": GROUP_LABELS["Basal_chordates"]},
                        {
                            "label": "Vertebrates",
                            "children": [
                                {"group": "Fish", "label": GROUP_LABELS["Fish"]},
                                {"group": "Amphibians", "label": GROUP_LABELS["Amphibians"]},
                                {
                                    "label": "Amniotes",
                                    "children": [
                                        {"group": "Sauropsids", "label": GROUP_LABELS["Sauropsids"]},
                                        {
                                            "label": "Mammalia",
                                            "children": [
                                                {
                                                    "group": "Monotremes_Marsupials",
                                                    "label": GROUP_LABELS["Monotremes_Marsupials"],
                                                },
                                                {
                                                    "label": "Placentalia",
                                                    "children": [
                                                        {
                                                            "group": "Mammals_non_primate",
                                                            "label": GROUP_LABELS["Mammals_non_primate"],
                                                        },
                                                        {
                                                            "group": "Primates_Hominids",
                                                            "label": GROUP_LABELS["Primates_Hominids"],
                                                        },
                                                    ],
                                                },
                                            ],
                                        },
                                    ],
                                },
                            ],
                        },
                    ],
                },
            ],
        },
    ],
}


def clone_tree_template() -> dict:
    return deepcopy(TREE_TEMPLATE)


def infer_broad_group(species_name: str) -> str:
    species_name = str(species_name)
    genus = species_name.split("_", 1)[0].lower()

    if genus in FUNGI_GENERA:
        return "Fungi"
    if genus in INVERTEBRATE_GENERA:
        return "Invertebrates"
    if genus in BASAL_CHORDATE_GENERA:
        return "Basal_chordates"
    if genus in EXPLICIT_FISH_GENERA:
        return "Fish"
    if genus in AMPHIBIAN_GENERA:
        return "Amphibians"
    if genus in SAUROPSID_GENERA:
        return "Sauropsids"
    if genus in PRIMATE_GENERA:
        return "Primates_Hominids"
    if genus in MONOTREME_MARSUPIAL_GENERA:
        return "Monotremes_Marsupials"
    if genus in MAMMAL_GENERA:
        return "Mammals_non_primate"

    # In the current Compara set, remaining unmatched non-outgroup vertebrates are
    # overwhelmingly fishes. Keep this fallback to preserve coverage while still
    # separating explicit amphibians and avoiding known invertebrate misclassification.
    return "Fish"


def normalize_group_value(group: object, species_name: Optional[object] = None) -> str:
    group_str = "" if pd.isna(group) else str(group)
    if group_str == "Fish_Amphibians":
        if species_name is not None and not pd.isna(species_name):
            refined = infer_broad_group(str(species_name))
            if refined in {"Fish", "Amphibians"}:
                return refined
        return LEGACY_TO_REFINED[group_str]
    if group_str in GROUP_ORDER or group_str == "Other":
        return group_str
    if species_name is not None and not pd.isna(species_name):
        return infer_broad_group(str(species_name))
    return group_str


def normalize_species_meta(species_meta: pd.DataFrame, group_col: str = "taxonomy_class") -> pd.DataFrame:
    df = species_meta.copy()
    if "species" not in df.columns:
        return df
    if group_col in df.columns:
        df[group_col] = [
            normalize_group_value(group, species)
            for group, species in zip(df[group_col], df["species"])
        ]
    if "broad_group" in df.columns:
        df["broad_group"] = [
            normalize_group_value(group, species)
            for group, species in zip(df["broad_group"], df["species"])
        ]
    return df


def normalize_group_summary_index(df: pd.DataFrame) -> pd.DataFrame:
    renamed = []
    for group in df.index:
        refined = normalize_group_value(group)
        renamed.append(refined)
    out = df.copy()
    out.index = renamed
    if out.index.has_duplicates:
        out = out.groupby(level=0, sort=False).mean()
    return out
