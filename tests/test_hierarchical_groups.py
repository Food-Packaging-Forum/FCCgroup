import pytest

from fccgroup.constants import CAS_COLUMN
from tests.conftest import _repo_root

chemical_groups = {
'Not groupable': {
    'UVCBs_Mixtures_Polymers':{
        'UVCBs': {'UVCB - Biological origin', 'UVCB - Process based'},
        'Mixtures':{},
        'Polymers':{},
        },
    },
'Inorganic compound': {
    'Inorganic salt (no C)': {},
    'Elements (loose)': {},
    'Inorganic other': {},
    },
'Contains Organometallics': {
    'Organometallic salt': {},
    'Contains Organo B': {},
    'Contains Organo Si': {'Organosilanes'},
    'Contains Organo Sn': {},
    'Organometallic other': {},
    },
'Organic (only nonmetals)': {
    'Hydrocarbons',
    'Organo O (strict)',
    'Organo N (strict)',
    'Organo P (strict)',
    'Organo S (strict)',
    'Organo X (strict)',
    'Organic mixed'
    },
'Organic C': {
    'Contains Hydrocarbons':{
        'Hydrocarbons': {
            'Aliphatic hydrocarbons': {'Alkanes', 'Alkenes', 'Alkynes'},
            'Aromatic hydrocarbons': {
                'Fully aromatic hydrocarbons':{},
                'Biphenyls/Terphenyls':{}, 
                'PAH derivatives hydrocarbon':{'PAHs'},
                'Benzoids':{},
                },
            },
        'Biphenyl/Terphenyl derivatives': {'Biphenyls/Terphenyls'},
        'PAH derivatives':{
            'PAH derivatives hydrocarbon':{'PAHs'},
            },
        'Abietic acid derivatives (tricyclic diterpenes backbone)':{},
        'Steroid-like compounds':{},
        },		
    'Contains Organo N': {
        'Contains amine': {
            'Contains primary amine': {'Primary aromatic amines',},
            'Contains secondary amine': {'Secondary aromatic amines'},
            'Contains tertiary amine': {'Tertiary aromatic amines'},
            'Contains p-phenylenediamine':{},
            'Aromatic amines (loose)':{'Aromatic amines',},
            'Aliphatic amines (loose)':{'Aliphatic amines',},
            },
        'Contains quarternary ammonium ion': {},
        'Contains imine': {'Aromatic imines', 'Aliphatic imines'},
        'Contains nitrile': {'Aromatic nitriles', 'Aliphatic nitriles'},
        'Contains isocyanide':{},
        'Contains cyanate':{},
        'Contains isocyanate':{},
        'Contains amide':{
            'Aromatic amides (loose)':{'Aromatic amides',}, 
            'Aromatic amides (loose)':{'Aliphatic amides'},
             },
        'Contains urea':{},
        'Contains nitroso':{},
        'Contains nitrosamine':{},
        'Contains nitro':{},
        'Contains azo':{},
        'Contains amino acid derivatives':{'Contains amino acid'},
        'Contains N-heterocycle': {
            'Imidazole derivatives':{'Imidazole'},
            'Pyrazole derivatives':{'Pyrazoles'},
            'Triazines derivatives':{
                'Triazines':{'Triazines (1,2,3)', 'Triazines (1,2,4)', 'Triazines (1,3,5)'},
                'Triazines (1,2,3) derivatives':{},
                'Triazines (1,2,4) derivatives':{},
                'Triazines (1,3,5) derivatives':{},
                },
            'Triazoles derivatives': {'Triazoles (1,2,3) derivatives', 'Triazoles (1,2,4) derivatives'},
            'Benzotriazines derivatives':{
                'Benzotriazines':{'Benzotriazines (1,2,3)', 'Benzotriazines (1,2,4)'}, 
                'Benzotriazines (1,2,3) derivatives':{},
                'Benzotriazines (1,2,4) derivatives':{},
                },				
            'Benzotriazoles derivatives':{'Benzotriazoles'},
            'Benzooxazole derivatives': {'Benzooxazoles'},
            'Benzoisooxazole derivatives':{'Benzoisooxazoles'}, 
            'Acridanone derivatives': {'Acridanone'},
            },
        },
    'Contains Organo O': {
        'Organo O (strict)':{},
        'Contains C=O': {
            'Contains aldehyde':{'Aliphatic aldehydes','Aromatic aldehydes',},
            'Contains ketone': {
                'Benzoquinone backbone': {
                    'Benzoquinone 1,2 backbone': {},
                    'Benzoquinone 1,4 backbone': {},
                    'Hindered quinone derivatives': {},
                    'Benzoquinone derivatives': {
                        'Benzoquinone 1,2 derivatives',
                        'Benzoquinone 1,4 derivatives',
                        },
                    },
                'Aliphatic ketones':{},
                'Aromatic ketones':{},
                },
            'Aliphatic carbonyls (loose)':{
                'Aliphatic aldehydes':{},
                'Aliphatic ketones':{},
                },
            'Aromatic carbonyls (loose)':{
                'Aromatic aldehydes':{},
                'Aromatic ketones':{},
                },
            },
        'Contains alcohol': {
            "Aliphatic alcohols":{},
            "Aromatic alcohols":{},
            'Contains phenol':{'Alkyl phenols', 'Hindered phenols', 'Bisphenols'},
            },
        'Contains phenol derivatives':{
            'Contains phenol':{},
            'Alkyl phenol derivatives':{'Alkyl phenols'},
            'Hindered phenol derivatives':{'Hindered phenols'},
            'Bisphenol derivatives':{
                'Bisphenols': {'Bisphenol M', 'Bisphenol P', 'Bisphenol bridge', 'Bisphenol CS bridge'},
                'Bisphenol M derivatives':{'Bisphenol M', },
                'Bisphenol P derivatives':{'Bisphenol P',},
                'Bisphenol bridge derivatives':{'Bisphenol bridge',},
                },
            },
        'Contains ether':{
            'Contains epoxide': {},
            'Xanthenes derivatives': {},
            'Glycerol derivatives': {'Triglyceride'},
            'Aromatic ethers': {},
            'Aliphatic ethers': {},
            },
        'Aliphatic ethers or alcohols (loose)':{
            'Aliphatic ethers':{},
            'Aliphatic alcohols':{},
            },   
        'Aromatic ethers or alcohols (loose)':{
            'Aromatic ethers':{},
            'Aromatic alcohols':{},
            },    
        'Contains carboxylic acid derivatives': {
            'Contains carboxylic acid': {'Aliphatic carboxylic acids', 'Aromatic carboxylic acids'},
            'Contains ester': {'Triglyceride'},
            'Aliphatic carboxylic acid esters and salts': {},
            'Adipic acid esters': {},
            'Sebacic acid esters': {},
            'Citric acid esters': {},
            'Maleic and fumaric acid esters': {},
            'Stearic acid esters': {},
            'Aromatic carboxylic acid esters and salts': {},
            'Contains benzoic acid derivatives': {
                'Contains benzoic acid': {},
                'Contains benzoic ester': {},
                'Benzoates derivatives': {'Benzoates',},
                'Parabens derivatives': {'Parabens',},
                'Salicylates derivatives': {'Salicylates',},
                'Cyclohexane-1,2-dicarboxylates': {},
                'Ortho-phthalates derivatives (no rings)': {'Ortho-phthalates'},
                'Isophthalates derivatives (no rings)': {'Isophthalates'},
                'Terephthalates derivatives (no rings)': {'Terephthalates'},
                'Trimellitates': {},
                'Trimesitates': {},
                'Hemimellitates': {},
                'Pyromellitates': {},                      
                },
            },
        'Contains O-heterocycle':{
            'Contains epoxide': {},
            'Benzofuran derivatives': {'Benzofurans'},
            },
        },
    'Contains Organo P': {
        'Contains C~P':{},
        'Contains phosphates':{
            'Organophosphates':{'Organophosphates (alkyl)','Organophosphates (aryl)'},
            },
        'Contains phophonium':{},
        'Contains P-heterocycle':{},         
        'Organophosphites':{},
        'Organophosphites tautomers':{},
        'Organophosphonates':{},
        'Dialkyl phosphite':{},
        'Organothiophosphates':{},
        },
    'Contains Organo S': {
        'Organo S (strict)': {'Aromatic sulfides','Aliphatic sulfides','Aromatic thiols','Aliphatic thiols',
            'Aromatic thioaldehydes','Aliphatic thioaldehydes','Aromatic thioketones','Aliphatic thioketones',},
        'Contains C~S':{},
        'Contains thiol':{'Aromatic thiols','Aliphatic thiols'},
        'Contains sulfide':{'Thiodiproprionates','Aromatic sulfides','Aliphatic sulfides'},
        'Contains disulfides':{},
        'Contains thioaldehyde':{'Aromatic thioaldehydes','Aliphatic thioaldehydes'},
        'Contains thioketone':{'Aromatic thioketones','Aliphatic thioketones'},
        'Contains thioester':{},
        'Contains thiosulfinate':{},
        'Contains dithioester':{},
        'Contains dithiocarbamate':{},
        'Contains thioamide':{},
        'Contains thiourea':{},
        'Contains sulfoxide':{},
        'Contains sulfinate':{'Aromatic sulfinates','Aliphatic sulfinates'},
        'Contains sulfon':{'Aromatic sulfons','Aliphatic sulfons'},
        'Contains sulfonyl halide':{},
        'Contains sulfonamide':{},
        'Contains thiosulfonate':{},
        'Contains sulfate':{'Organosulfate'},
        'Contains sulfonate':{'Aromatic sulfonates','Aliphatic sulfonates'},
        'Contains sulfenes':{},
        'Contains thiocyanate':{},
        'Contains isothiocyanate':{},
        'Contains S-C-N':{},
        'Contains S-heterocycle':{
            'Thiophene derivatives':{'Thiophenes',},
            'Benzothiophene derivatives':{'Benzothiophenes',},
            'Benzothiazoles derivatives':{'Benzothiazoles',},
            'Benzisothiazoles derivatives':{'Benzisothiazoles',},
            'Thioxantenes derivatives':{},
            'Thioxantones derivatives':{},
            },
        'Triarylsulfonium derivatives':{},
        },
    'Contains Organo X':{
        'Contains Organo F':{
            'Organo F (strict)':{},
            'Contains C~F': {'PFAS'},
            },
        'Contains Organo Cl':{
            'Organo Cl (strict)':{'PCBs', 'DDT derivatives', 'Chlorinated alkynes', 'Chlorinated alkenes', 'Chlorinated alkanes'},
            'Contains C~Cl': {
                'Contains aliphatic C~Cl':{'Chlorinated alkynes', 'Chlorinated alkenes', 'Chlorinated alkanes'},
                'Contains aromatic c~Cl':{
                    'Chlorobenzenes':{},
                    'Chlorophenol derivatives': {'Chlorophenols'},
                    'PCBs derivatives': {'PCBs'},
                    'PCDDs derivatives': {'PCDDs'},
                    'PCDEs derivatives': {'PCDEs'},
                    'PCDFs derivatives':{'PCDFs'},
                    'DDT derivatives':{},},
                },
            },
        'Contains Organo Br':{
            'Organo Br (strict)':{'PBBs'},
            'Contains C~Br':{ 
                'Bromophenol derivatives': {'Bromophenols'},
                'PBBs derivatives': {'PBBs'},
                'PBDEs derivatives':{'PBDEs'},
                'PBDDs derivatives': {'PBDDs'},
                'PBDFs derivatives': {'PBDFs'},
                },
            },
        'Contains Organo I':{'Organo I (strict)', 'Contains C~I','Diaryliodonium derivatives'},
        },
    } 
}

def test_hierarchical_groups(universe):
    available_columns = {col[1] for col in universe.columns if isinstance(col, tuple) and len(col) > 1}
    required_groups = _collect_group_labels(chemical_groups)
    if not required_groups.issubset(available_columns):
        missing_count = len(required_groups - available_columns)
        pytest.skip(f"Universe fixture is missing {missing_count} hierarchy columns required by this test")

    not_matching_list = helper_hierarchical_groups(universe, chemical_groups)
    assert not not_matching_list, f"Some groups do not match their subgroups and are saved in: {save_not_matching_list_to_excel(not_matching_list)}"

def _collect_group_labels(chemical_hierarchy_dict):
    labels = set()
    for group, subgroups in chemical_hierarchy_dict.items():
        labels.add(group)
        if isinstance(subgroups, set):
            labels.update(subgroups)
        elif isinstance(subgroups, dict):
            labels.update(_collect_group_labels(subgroups))
    return labels


def helper_hierarchical_groups(df, chemical_hierarchy_dict, not_matching_list=None):
    if not_matching_list is None:
        not_matching_list = []

    # Check if the hierarchical groups are present in the universe DataFrame
    for group, subgroups in chemical_hierarchy_dict.items():
        # If children are a list, that means that we are on an end
        if isinstance(subgroups, set):
            # If list is empty
            if not subgroups:
                continue
            check_group_match_subgroups(df, group, subgroups, not_matching_list)
        # If children are a dict, we need to recursively check each subgroup
        elif isinstance(subgroups, dict):
            helper_hierarchical_groups(df, subgroups, not_matching_list)
        else:
            raise ValueError(f"Unexpected type for subgroups: {type(subgroups)} in group {group}")
    return not_matching_list

def check_group_match_subgroups(df, group, subgroups, not_matching_list):
    parent_column = find_column(df, group)
    filter1 = ~(df[parent_column] > 0)
    for subgroup in subgroups:
        child_columns = find_column(df, subgroup)
        filter2 = df[child_columns] > 0
        non_matching_df = df.loc[filter1 & filter2, [("Identifier", CAS_COLUMN)]]
        if not non_matching_df.empty:
            for _, row in non_matching_df.iterrows():
                not_matching_list.append({
                    "Parent": group,
                    "Child": subgroup,
                    "CASRN": row[("Identifier", CAS_COLUMN)]
                })

def find_column(df, column_name):
    """
    Find the column in the DataFrame that matches the given name.
    """
    for col in df.columns:
        if col[1] == column_name:
            return col
    raise ValueError(f"Column {column_name} not found in DataFrame.")

def save_not_matching_list_to_excel(not_matching_list):
    import pandas as pd
    df = pd.DataFrame(not_matching_list)
    df2 = df.groupby(["Parent", "Child"]).count()
    excel_filename = _repo_root() / "tests" / "Hierarchical_test_comparison.xlsx"
    with pd.ExcelWriter(excel_filename) as writer:
        df.to_excel(writer, sheet_name="Not Matching", index=False)
        df2.to_excel(writer, sheet_name="Summary")
    return excel_filename