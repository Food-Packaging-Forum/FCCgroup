"""
Chemical constants and reference data.

This module contains periodic table data, functional synonyms, list mappings,
and regex patterns used throughout the FCCgroup application.

Types:
    elements: Dict[str, Dict[str, int]] - Element groups and their atomic numbers
    synonym_lists: Dict[str, Dict[str, List[str]]] - Function categories and regex patterns
    individual_lists: Dict[str, List[str]] - List IDs and their names
    regex_combination_dictionary: Dict[str, List[str]] - Regex group combinations
"""

from typing import Final, Dict, List

# ============================================================================
# PERIODIC TABLE ELEMENTS
# ============================================================================
# Maps element groups to their symbols and atomic numbers
# Structure: {group_name: {symbol: atomic_number, ...}, ...}

elements: Final[Dict[str, Dict[str, int]]] = {
    "alkali_metals": {
        "Li": 3,
        "Na": 11,
        "K": 19,
        "Rb": 37,
        "Cs": 55,
        "Fr": 87,
        },

    "alkaline_metals": {
        "Be": 4,
        "Mg": 12,
        "Ca": 20,
        "Sr": 38,
        "Ba": 56,
        "Ra": 88,
        },

    "transition_metals": {
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
        "Y": 39,
        "Zr": 40,
        "Nb": 41,
        "Mo": 42,
        "Tc": 43,
        "Ru": 44,
        "Rh": 45,
        "Pd": 46,
        "Ag": 47,
        "Cd": 48,
        "Hf": 72,
        "Ta": 73,
        "W": 74,
        "Re": 75,
        "Os": 76,
        "Ir": 77,
        "Pt": 78,
        "Au": 79,
        "Hg": 80,
        "Rf": 104,
        "Db": 105,
        "Sg": 106,
        "Bh": 107,
        "Hs": 108,
        "Mt": 109,
        "Ds": 110,
        "Rg": 111,
        "Cn": 112,
        },

    "posttransition_metals": {
        "Al": 13,
        "Ga": 31,
        "In": 49,
        "Sn": 50,
        "Tl": 81,
        "Pb": 82,
        "Bi": 83,
        "Po": 84,
        "Nh": 113,
        "Fl": 114,
        "Mc": 115,
        "Lv": 116,
        "Ts": 117,
        },

    "lanthanoids": {
        "La": 57,
        "Ce": 58,
        "Pr": 59,
        "Nd": 60,
        "Pm": 61,
        "Sm": 62,
        "Eu": 63,
        "Gd": 64,
        "Tb": 65,
        "Dy": 66,
        "Ho": 67,
        "Er": 68,
        "Tm": 69,
        "Yb": 70,
        "Lu": 71,
        },

    "actinoids": {
        "Ac": 89,
        "Th": 90,
        "Pa": 91,
        "U": 92,
        "Np": 93,
        "Pu": 94,
        "Am": 95,
        "Cm": 96,
        "Bk": 97,
        "Cf": 98,
        "Es": 99,
        "Fm": 100,
        "Md": 101,
        "No": 102,
        "Lr": 103,
        },

    "metalloids": {
        "B": 5,
        "Si": 14,
        "Ge": 32,
        "As": 33,
        "Sb": 51,
        "Te": 52,
        "At": 85,
        "Ts": 117
        },

    "non_metals": {
        "H": 1,
        "C": 6,
        "N": 7,
        "O": 8,
        "P": 15,
        "S": 16,
        "Se": 34,
        # Halogens
        "F": 9,
        "Cl": 17,
        "Br": 35,
        "I": 53,
        # Noble Gases
        'He':2,
        'Ne': 10,
        'Ar': 18,
        'Kr': 36,
        'Xe': 54,
        'Rn': 86,
        'Og': 118,
        },

    "halogens": {
        "F": 9,
        "Cl": 17,
        "Br": 35,
        "I": 53
        },
    
    "noble_gases": {
        'He':2,
        'Ne': 10,
        'Ar': 18,
        'Kr': 36,
        'Xe': 54,
        'Rn': 86,
        'Og': 118,
        },
    
    "metals": {
        # Alkali Metals (Group 1)
        "Li": 3,
        "Na": 11,
        "K": 19,
        "Rb": 37,
        "Cs": 55,
        "Fr": 87,
        # Alkaline Earth Metals (Group 2)
        "Be": 4,
        "Mg": 12,
        "Ca": 20,
        "Sr": 38,
        "Ba": 56,
        "Ra": 88,
        # Transition Metals
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
        "Y": 39,
        "Zr": 40,
        "Nb": 41,
        "Mo": 42,
        "Tc": 43,
        "Ru": 44,
        "Rh": 45,
        "Pd": 46,
        "Ag": 47,
        "Cd": 48,
        "Hf": 72,
        "Ta": 73,
        "W": 74,
        "Re": 75,
        "Os": 76,
        "Ir": 77,
        "Pt": 78,
        "Au": 79,
        "Hg": 80,
        "Rf": 104,
        "Db": 105,
        "Sg": 106,
        "Bh": 107,
        "Hs": 108,
        "Mt": 109,
        "Ds": 110,
        "Rg": 111,
        "Cn": 112,
        # Post-transition Metals
        "Al": 13,
        "Ga": 31,
        "In": 49,
        "Sn": 50,
        "Tl": 81,
        "Pb": 82,
        "Bi": 83,
        "Po": 84,
        "Nh": 113,
        "Fl": 114,
        "Mc": 115,
        "Lv": 116,
        "Ts": 117,
        # Lanthanoids
        "La": 57,
        "Ce": 58,
        "Pr": 59,
        "Nd": 60,
        "Pm": 61,
        "Sm": 62,
        "Eu": 63,
        "Gd": 64,
        "Tb": 65,
        "Dy": 66,
        "Ho": 67,
        "Er": 68,
        "Tm": 69,
        "Yb": 70,
        "Lu": 71,
        # Actinoids
        "Ac": 89,
        "Th": 90,
        "Pa": 91,
        "U": 92,
        "Np": 93,
        "Pu": 94,
        "Am": 95,
        "Cm": 96,
        "Bk": 97,
        "Cf": 98,
        "Es": 99,
        "Fm": 100,
        "Md": 101,
        "No": 102,
        "Lr": 103
    },

}

# ============================================================================
# FUNCTIONAL GROUP SYNONYMS
# ============================================================================
# Maps function categories to regex patterns for matching
# Structure: {category_type: {function_name: [pattern1, pattern2, ...], ...}, ...}

synonym_lists: Final[Dict[str, Dict[str, List[str]]]] = {
    "function_synonyms_regex": {
        'Antioxidant': [
            'anti.?degrad',
            'anti.?ozon',
            'anti.?radical',
            'anti.?oxidant',
            'oxidation.?inhibit',
            'oxidation.?stabil',
            'peroxide.?decompos',
            'peroxide.?scaveng',
            'radical.?deactivat',
            'radical.?scaveng',
            'radical.?trap',
            '(?<!viscosity)(?<!friction)(?:^|[\n .,\\-])reduc.{0,10}agent'
            ],

        'Antistatic Agent': [
            'anti.?static',
            'charge.dissipat',
            'dissipat.{0,10}charge',
            'conduct.{0,10}additive',
            'conduct.{0,10}agent'
            ],

        'Biocide': [
            'anti.?bacteria',
            'anti.?biotic',
            'anti.?fung',
            'anti.?microb',
            'anti.?viral',
            'bactericid',
            'bakterizid',
            'biocid',
            'biozid',
            'fungicid',
            'fungizid'
            ],

        'Blowing Agent': [
            'blow.{0,10}agent',
            'expander',
            '(?<!anti)(?<!anti[-. ])foam.{0,10}agent',
            'pneumatogens'
            ],

        'Catalyst': ['cataly'],

        'Colorant': [
            'bleach.{0,10}agent',
            'brighten',
            'clarif',
            '(?<!(?:dis|non))(?:^|[ -,;.:])coloring',
            'colorant',
            '(?<!(?:dis|non))(?:^|[ -,;.:])colouring',
            'colourant',
            'dye',
            'paint',
            'pigment(?!ed)',
            'whiten'
            ],

        'Crosslinking Agent': [
            'accelerator',
            'cross.?link(?:.{0,10}agent|er)',
            'cur.{0,10}agent',
            'hardener',
            'vulcani(?:[sz]er|.{0,10}(?:agent|accel|activat))'
            ],

        'Filler': [
            'filler', 
            'filling'
            ],

        'Flame Retardant': [
            'anti.?fire',
            'anti.?flame',
            'anti.?ignit',
            'anti.?smoke',
            'fire.?protect',
            'fire.?retard',
            'fire.?stabil',
            'flame.?retard'
            ],

        'Heat Stabilizer': [
            'heat.?stabili[csz]er',
            '(?:therm|temp).{0,10}stabili[csz]er'
            ],

        'Impact Modifier': [
            'anti.?brittl',
            'durab.{0,10}agent',
            'impact.?modifier'
            ],

        'Initiator': ['initiat'],

        'Intermediates': [
            'intermed', 
            'intermediate'
            ],

        'Light Stabilizer': [
            'HALS',
            'light.?stabili[csz]er',
            'quencher',
            'ultra.?violet.?absor',
            'ultra.?violet.?stabil',
            'uv.?absorber',
            'uv.?stabili[csz]er'
            ],

        'Lubricant': [
            'anti.?block.{0,10}agent',
            'anti.?tack.?agent',
            'friction.?min',
            'min.{0,10}friction',
            'friction.*reduc',
            'reduc.{0,10}friction',
            'grease',
            '(?<!in )lubric',
            'mold.?release',
            'release.?agent',
            'slip.?agent',
            'slip.?promot'
            ],

        'Monomer': [
            'building.?block', 
            'monomer.(?!free)'
            ],

        'Nucleating Agent': [
            'crystalli[scz].{0,10}agent', 
            'nucleat'
            ],

        'Odor Agent': [
            'deodorant',
            'fragranc',
            'odor.{0,10}agent',
            'odorant'
            ],

        'Other Processing Aids': [
            'acid.?scavenger',
            'buffer',
            'ph.?regulat',
            'ph.?adjust',
            'neutral.{0,10}agent',
            'adhesi',
            '%emul%',
            'flocculant',
            'polish'
            ],

        'Plasticizer': [
            'distensibility',
            'flex.{0.10}mod',
            'plasticiser',
            'plasticizer'
            ],

        'Solvent': [
            '.{0,200}solvent (?! (?:dye|pigment))'
        ],

        'Viscosity Modifier': [
            'flexibilis',
            'flexibiliz',
            'rheolog.{0,10}modif',
            'viscos.{0,10}adj',
            'viscos.{0,10}control',
            'viscos.{0,10}red',
            'viscos.{0,10}reg',
            'viscos.{0,10}modif'
            ]
    },

    # FPF's Synonym library
    "function_synonyms": {
        "PFASs": [
            "Perfluoroalkyl and Polyfluoroalkyl Substances: Surfactants",
            "Per- and polyfluoroalkyl substances",
            "Per-/poly-fluorinated compounds (PFASs)",
            "Perfluoroalkyl and Polyfluoroalkyl Substances (PFASs)",
            "Perfluoroalkyl and Polyfluoroalkyl Substances (PFASs) "],
            
        "VOCs":[
            "Volatile Organic Compounds (VOC)",
            "Volatile organic compounds",
            "Aprotic solvents",
            "VOCs metabolites",
            "Volatile Organic Compounds (VOC) Metabolites",
            "Volatile organic compounds: Benzene",
            "Aldehydes",
            "Volatile organic compounds: Trihalomethanes",
            "adducts of hemoglobin",
            "Acrylamide"],
            
        "Metals":[
            "Metals and Metalloids",
            "Metals and trace elements",
            "Arsenic",
            "Heavy metals",
            "Metals",
            "Metals and trace elements: Mercury",
            "Cadmium",
            "Metals and trace elements: Arsenic",
            "Chromium",
            "Lead",
            "Mercury and its organic compounds"],
            
        "Flame retardants":[
            "Flame Retardant Metabolites",
            "Flame retardants",
            "Polybrominated Diphenyl Ethers (PBDEs)",
            "Polybrominated Diphenyl Ethers and PBB 153 (Pooled Samples after 2004)",
            "Flame retardants: Polybrominated diphenyl ethers (PBDEs)"],
            
        "Pesticides (herbicides, insecticides, fungicides, biocides, disinfection)": [
            "Carbamate Pesticide Metabolites",
            "Organochlorines",
            "Pesticides",
            "Pyrethroid pesticides metabolites",
            "Herbicides",
            "Pyrethroid Metabolites",
            "Organochlorines: Chlordane",
            "Pesticides (pyrethroids)",
            "Organochlorine Pesticides",
            "Sulfonyl Urea Herbicides",
            "Organochlorines: Dichlorodiphenyltrichloroethane (DDT)",
            "Organophosphate Pesticides",
            "Fungicides and Metabolites",
            "Organochlorines: Hexachlorocyclohexane (HCH)",
            "Other Pesticides",
            "Herbicides and Metabolites",
            "Organochlorines: Toxaphene",
            "Pyrethroid Pesticides",
            "Insect Repellent and Metabolites",
            "Pesticides: Atrazine",
            "Neonicotinoid Insecticides",
            "Pesticides: Carbamates",
            "Organochlorine Pesticide Metabolites ",
            "Organochlorine Pesticide Metabolites",
            "Pesticides: 2,4-Dichlorophenoxyacetic acid (2,4-D)",
            "Organochlorine Pesticides and Metabolites",
            "Organochlorine Pesticides and Metabolites ",
            "Pesticides: Ethylene bisdithiocarbamates",
            "Organophosphorus Insecticides: Dialkyl Phosphate Metabolites",
            "Pesticides: ortho-Phenylphenol (OPP)",
            "Organophosphorus Insecticides: Specific Metabolites",
            "Pesticides: Organophosphates",
            "Disinfection By-Products",
            "Pesticides: Pyrethroids",
            "Chlorophenols"],
            
        "PAHs":[
            "Polycyclic Aromatic Hydrocarbon Metabolites",
            "Polycyclic aromatic hydrocarbons: Benzo[a]pyrene",
            "Polycyclic Aromatic Hydrocarbons (PAHs)",
            "PAHs metabolites",
            "Polycyclic aromatic hydrocarbons: Chrysene",
            "Polycyclic aromatic hydrocarbons: Fluoranthene",
            "Polycyclic aromatic hydrocarbons: Fluorene",
            "Polycyclic aromatic hydrocarbons: Napthalene",
            "Polycyclic aromatic hydrocarbons: Phenanthrene",
            "Polycyclic aromatic hydrocarbons: Pyrene"],
            
        "Amines":[
            "Aromatic Amines",
            "Anilines and MOCA",
            "Aromatic Diamines",
            "Heterocyclic Amines",
            "Volatile N-Nitrosamines (VNAs)"],

        "Tobacco-related":[
            "Tobacco Alkaloids and Metabolites",
            "Nicotine",
            "Environmental tobacco smoke metabolites",
            "Tobacco-Specific Nitrosamines (TSNAs)"],

        "Natural origin & food related":[
            "Phytoestrogens and Metabolites",
            "Mycotoxins"],

        "Others":["Diesel exhaust"]
    }
}

# ============================================================================
# INDIVIDUAL LIST MAPPINGS  
# ============================================================================
# Maps function subcategories from specific lists
# Structure: {category_type: {function_id: [description, ...], ...}, ...}

individual_lists: Final[Dict[str, Dict[str, List[str]]]] = {
    "function_synonyms_g4": {
        # from G04
        'Hindered_phenols':                                 ["substances containing 4-TBP"],
        'simple_aliphatic_ketones_acyclic':                 ["Linear aliphatic ketones"],
        'simple_aliphatic_amides_acyclic':                  ["Aliphatic primary amides"],
        'Nitrates':                                         ["Alkyl nitrates"],
        'Aralkyladehydes':                                  ["Aralkylaldehydes"],
        'Ethers':                                           ["aromatic ethers"],
        'Cyclic_acetals':                                   ["Cyclic acetals from aldehydes"],
        'Cyclic_ethers':                                    ["Cyclic ethers"],
        'Aliphatic carboxylic acids':                       ["Branched/cyclic dialiphatic ethers (excluding alpha,beta-unsaturated ethers)"],
        'Dyes':                                             ["Diazo amino hydroxyl naphthalenedisulfonic acid dyes"],
        'Peroxides':                                        ["dibenzoyl peroxide derivatives"],
        'Dihydropurinediones':                              ["Dihydropurinedione derivatives"],
        'Ethanediols':                                      ["1,2-ethanediols and their carbonates"],
        'Iso-phthalates, tere-phthalates, trimellitates':   ["Isophthalates, Terephthalates and Trimellitates"],
        'contains_Manganese':                               ["Simple manganese compounds"],
        'Ortho-phthalates':                                 ["Ortho-phthalates"],
        'Parabens':                                         ["Paraben acid, salts and esters"],
        'Pyrazoles':                                        ["Pyrazoles"],
        'Salicylates':                                      ["Salicylic acid, its salts and alkylated derivatives, Salicylate esters"],
        'contains_Vanadium':                                ["simple vanadium compounds"],
        'PFAS':                                             ["Fluorinated aliphatic hydrocarbons"],
        'Bisphenols derivatives':                           ['Bisphenol A (BPA) derivatives','Bisphenol F (BPF) derivative','Miscellaneous bisphenols','Tetrabromobisphenol A derivatives','Bisphenol S (BPS) derivatives','Other aliphatic- or aryl-bridged bisphenol derivatives','Bisphenol AF (BPAF) derivatives'],
    },

    "PlastChem_columns": [
        'aromatic_amines', 
        'aralkyladehydes', 
        'alkylphenols', 
        'salicylate_esters', 
        'aromatic_ethers', 
        'bisphenols', 
        'orthophthalates', 
        'benzothiazole', 
        'benzotriazoles', 
        'organometallics', 
        'parabens', 
        'azodyes', 
        'acetophenones_benzophenones', 
        'chlorinated_paraffins', 
        'PFASs', 
        'UVCBs', 
        'polymers', 
        'mixtures', 
        'inorganic_compounds', 
        'DDT_DDE_DDD', 
        'dioxins', 
        'PBDEs', 
        'PCBs', 
        'PBDD_PBDF_PCDD_PCDF', 
        'aldehydes_simple', 
        'alkanes', 
        'alkenes', 
        'alkynes', 
        'alkane_ethers', 
        'aliphatic_ketones', 
        'aliphatic_primary_amides', 
        'alkyl_nitrates', 
        'aromatic_hydrocarbons', 
        'carboxylic_acids_salts', 
        'cyclic_acetals', 
        'cyclic_ethers', 
        'dialiphatic_ethers_excluding_unsatured', 
        'diazo_amino_hydroxyl_naphthalenedisulfonic_acid_dyes', 
        'dibenzoyl_peroxide_derivatives', 
        'dihydropurinediones', 
        'ethanediols', 
        'isophthalates_terephthalates_trimellitates', 
        'ketones_simple', 
        'organophosphates', 
        'phenolic_antioxidants', 
        'polychlorinated_naphthalenes', 
        'pyrazoles', 
        'salicyclic_acid', 
        'silanes_siloxanes_silicones', 
        'homo_CCO', 
        'homo_CF2', 
        'homo_CF2CF2O', 
        'homo_CF2O', 
        'homo_CH2', 
        'arsenic', 
        'cadmium', 
        'mercury', 
        'chromium', 
        'antimony', 
        'tin', 
        'bromo', 
        'chloro', 
        'fluoro', 
        'iodo', 
        'manganese', 
        'magnesium', 
        'barium', 
        'nickel', 
        'lead', 
        'vanadium', 
        'tellurium', 
        'thallium', 
        'beryllium', 
        'selenium'
    ]
}

# ============================================================================
# REGEX COMBINATION RULES
# ============================================================================
# Maps parent group names to required component patterns
# Structure: {parent_group: [component1, component2, ...], ...}

regex_combination_dictionary: Final[Dict[str, List[str]]] = {
    "Metal_Metalloid": [
        'Metal',
        'Metalloid'
    ],
    "OrganoMetallic": [
        'Metal_Metalloid',
        'Carbon'
    ],

    'OrganoMetallic_salt': [
        'OrganoMetallic', 
        'contains_salt'
    ],
        
    'Inorganic_salts': [
        'Inorganic_noC', 
        'contains_salt'
    ]
}

__all__: Final[List[str]] = ["elements", "synonym_lists", "individual_lists", "regex_combination_dictionary"]