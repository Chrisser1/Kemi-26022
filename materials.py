import re
import math
from collections import defaultdict

# --- UNIFIED DATA REPOSITORIES ---

LIGAND_DATA = {
    # Anionic Ligands (IUPAC 2005)
    'F': {'name': 'fluorido', 'charge': -1, 'type': 'simple'},
    'Cl': {'name': 'chlorido', 'charge': -1, 'type': 'simple'},
    'Br': {'name': 'bromido', 'charge': -1, 'type': 'simple'},
    'I': {'name': 'iodido', 'charge': -1, 'type': 'simple'},
    'CN': {'name': 'cyanido', 'charge': -1, 'type': 'simple'},
    'OH': {'name': 'hydroxido', 'charge': -1, 'type': 'simple'},
    'O': {'name': 'oxido', 'charge': -2, 'type': 'simple'},
    'SO4': {'name': 'sulfato', 'charge': -2, 'type': 'simple'},
    'CO3': {'name': 'carbonato', 'charge': -2, 'type': 'simple'},
    'C2O4': {'name': 'oxalato', 'charge': -2, 'type': 'complex_prefix'},
    'SCN': {'name': 'thiocyanato-S', 'charge': -1, 'type': 'simple'},
    'NCS': {'name': 'isothiocyanato', 'charge': -1, 'type': 'simple'},
    'NO2': {'name': 'nitrito-N', 'charge': -1, 'type': 'simple'},
    'ONO': {'name': 'nitrito-O', 'charge': -1, 'type': 'simple'},
    # Neutral Ligands
    'H2O': {'name': 'aqua', 'charge': 0, 'type': 'simple'},
    'NH3': {'name': 'ammine', 'charge': 0, 'type': 'simple'},
    'CO': {'name': 'carbonyl', 'charge': 0, 'type': 'simple'},
    'NO': {'name': 'nitrosyl', 'charge': 0, 'type': 'simple'},
    'en': {'name': 'ethylenediamine', 'charge': 0, 'type': 'complex_prefix'},
    'py': {'name': 'pyridine', 'charge': 0, 'type': 'simple'},
}

# Reverse mapping from name to formula for name-to-formula conversion
NAME_TO_LIGAND = {v['name']: {'formula': k, 'charge': v['charge'], 'type': v['type']} for k, v in LIGAND_DATA.items()}
# Add common/older names for broader compatibility
NAME_TO_LIGAND['nitro'] = NAME_TO_LIGAND['nitrito-N']
NAME_TO_LIGAND['nitrito'] = NAME_TO_LIGAND['nitrito-O']
NAME_TO_LIGAND['thiocyanato'] = NAME_TO_LIGAND['thiocyanato-S']
# Add common non-standard names to make the parser more flexible
NAME_TO_LIGAND['fluoro'] = NAME_TO_LIGAND['fluorido']
NAME_TO_LIGAND['chloro'] = NAME_TO_LIGAND['chlorido']
NAME_TO_LIGAND['bromo'] = NAME_TO_LIGAND['bromido']
NAME_TO_LIGAND['iodo'] = NAME_TO_LIGAND['iodido']
NAME_TO_LIGAND['cyano'] = NAME_TO_LIGAND['cyanido']
NAME_TO_LIGAND['hydroxo'] = NAME_TO_LIGAND['hydroxido']


METAL_DATA = {
    'Fe': {'name': 'iron', 'anion_name': 'ferrate'},
    'Cu': {'name': 'copper', 'anion_name': 'cuprate'},
    'Pb': {'name': 'lead', 'anion_name': 'plumbate'},
    'Ag': {'name': 'silver', 'anion_name': 'argentate'},
    'Au': {'name': 'gold', 'anion_name': 'aurate'},
    'Sn': {'name': 'tin', 'anion_name': 'stannate'},
    'Co': {'name': 'cobalt', 'anion_name': 'cobaltate'},
    'Pt': {'name': 'platinum', 'anion_name': 'platinate'},
    'Ni': {'name': 'nickel', 'anion_name': 'nickelate'},
    'Cr': {'name': 'chromium', 'anion_name': 'chromate'},
    'Zn': {'name': 'zinc', 'anion_name': 'zincate'},
    'Mn': {'name': 'manganese', 'anion_name': 'manganate'},
    'Al': {'name': 'aluminum', 'anion_name': 'aluminate'},
    'Ti': {'name': 'titanium', 'anion_name': 'titanate'},
}
# Reverse mapping for name-to-formula
NAME_TO_METAL = {v['name']: k for k, v in METAL_DATA.items()}
NAME_TO_METAL.update({v['anion_name']: k for k, v in METAL_DATA.items()})


COUNTER_ION_DATA = {
    'K': {'name': 'potassium', 'charge': 1},
    'Na': {'name': 'sodium', 'charge': 1},
    'NH4': {'name': 'ammonium', 'charge': 1},
    'H': {'name': 'hydrogen', 'charge': 1}, # For acids
    'Cl': {'name': 'chloride', 'charge': -1},
    'Br': {'name': 'bromide', 'charge': -1},
    'I': {'name': 'iodide', 'charge': -1},
    'SO4': {'name': 'sulfate', 'charge': -2},
    'NO3': {'name': 'nitrate', 'charge': -1},
}
# Reverse mapping for name-to-formula
NAME_TO_COUNTER_ION = {v['name']: {'formula': k, 'charge': v['charge']} for k, v in COUNTER_ION_DATA.items()}

PREFIXES = {1: '', 2: 'di', 3: 'tri', 4: 'tetra', 5: 'penta', 6: 'hexa'}
COMPLEX_PREFIXES = {1: '', 2: 'bis', 3: 'tris', 4: 'tetrakis', 5: 'pentakis', 6: 'hexakis'}
PREFIX_VALUES = {v: k for k, v in PREFIXES.items()}
PREFIX_VALUES.update({v: k for k, v in COMPLEX_PREFIXES.items()})
PREFIX_VALUES[''] = 1 # Handle mono case

# Spectrochemical series (weak-field to strong-field ligands)
SPECTROCHEMICAL_SERIES = [
    'I', 'Br', 'SCN', 'Cl', 'F', 'OH', 'H2O', 'NCS', 'py', 'NH3', 'en', 'NO2', 'CN', 'CO'
]

FULL_FORMULA_RE = re.compile(r"""
    ^                                       # Anchor to the start of the string.
    (?P<cation>                             # Start of the optional named group for the cation.
        (?:                                 # Non-capturing group for cation types.
            \([A-Z][A-Za-z0-9]*\)\d* | # Polyatomic cation, e.g., (NH4)2
            [A-Z][a-z]?\d* # Simple atomic cation, e.g., K2
        )
    )?                                      # End of the optional cation group.
    (?P<complex>                            # Start of the required named group for the complex.
        \[.+?\]                             # The coordination sphere, non-greedy.
    )                                       # End of the complex group.
    (?P<anion_charge>                       # Start of the optional named group for anion/charge.
        (?:                                 # Non-capturing group for anion/charge types.
            \([A-Z][A-Za-z0-9]*\)\d* | # Polyatomic anion, e.g., (SO4)3
            [A-Z][a-z]?\d* | # Simple atomic anion, e.g., Cl3
            \d*[+\-]                        # A charge indicator, e.g., 2+ or -
        )
    )?                                      # End of the optional anion/charge group.
    $                                       # Anchor to the end of the string.
""", re.VERBOSE)

# --- HELPER FUNCTIONS ---
def _to_roman(num: int) -> str:
    """Converts an integer to a Roman numeral string in parentheses."""
    if not isinstance(num, int) or num < 0: return "(?)"
    if num == 0: return "(0)"
    val_map = [(10, 'X'), (9, 'IX'), (5, 'V'), (4, 'IV'), (1, 'I')]
    res = ''
    for val, sym in val_map:
        while num >= val:
            res += sym
            num -= val
    return f"({res})"

def _from_roman(s: str) -> int:
    """Converts a Roman numeral string (e.g., 'III') to an integer."""
    s = s.strip().upper()
    roman_map = {'I': 1, 'V': 5, 'X': 10}
    num = 0
    for i in range(len(s)):
        if i > 0 and roman_map[s[i]] > roman_map[s[i-1]]:
            num += roman_map[s[i]] - 2 * roman_map[s[i-1]]
        else:
            num += roman_map[s[i]]
    return num

def _parse_ion_formula(formula: str) -> tuple[str, int]:
    """Parses an ion formula like 'Cl3' or '(NH4)2' into its symbol and count."""
    if formula.startswith('('):
        match = re.match(r"\(([A-Z][A-Za-z0-9]*)\)(\d*)", formula)
        if not match: raise ValueError(f"Could not parse polyatomic ion: {formula}")
        ion, count_str = match.groups()
        return ion, int(count_str) if count_str else 1
    else:
        match = re.match(r"([A-Z][a-z]?)(\d*)", formula)
        if not match: raise ValueError(f"Could not parse simple ion: {formula}")
        ion, count_str = match.groups()
        return ion, int(count_str) if count_str else 1

def _parse_complex_name(name: str) -> tuple[str, int]:
    """Parses the name of a single complex ion and returns its formula and charge."""
    name = name.strip()
    
    # 1. Parse metal and oxidation state
    ox_state_match = re.search(r"\((\w+)\)$", name)
    if not ox_state_match: raise ValueError(f"Oxidation state not found in name: {name}")
    oxidation_state = _from_roman(ox_state_match.group(1))
    
    metal_part_end = ox_state_match.start()
    ligands_str = name[:metal_part_end]
    
    # Find the metal name by checking against known metal names (longest first)
    metal_name, metal_symbol = "", ""
    for m_name in sorted(NAME_TO_METAL.keys(), key=len, reverse=True):
        if ligands_str.endswith(m_name):
            metal_symbol = NAME_TO_METAL[m_name]
            ligands_str = ligands_str[:-len(m_name)]
            break
    if not metal_symbol: raise ValueError(f"Unknown metal in: {name}")

    # 2. Parse ligands from the string
    ligand_counts = defaultdict(int)
    sorted_ligand_names = sorted(NAME_TO_LIGAND.keys(), key=len, reverse=True)
    i = 0
    while i < len(ligands_str):
        found = False
        prefix_match = re.match(r"^(bis|tris|tetrakis)\((.+?)\)", ligands_str[i:])
        if prefix_match:
            prefix, ligand_name = prefix_match.groups()
            if ligand_name in NAME_TO_LIGAND:
                ligand_counts[ligand_name] += PREFIX_VALUES[prefix]
                i += prefix_match.end()
                found = True
        if not found:
            for lig_name in sorted_ligand_names:
                for prefix_str in sorted(PREFIX_VALUES.keys(), key=len, reverse=True):
                    if ligands_str[i:].startswith(prefix_str + lig_name):
                        ligand_counts[lig_name] += PREFIX_VALUES[prefix_str]
                        i += len(prefix_str) + len(lig_name)
                        found = True
                        break
                if found: break
        if not found and i < len(ligands_str):
            raise ValueError(f"Could not parse ligand string at: '{ligands_str[i:]}'")
        elif not found: break

    # 3. Calculate complex charge
    total_ligand_charge = sum(NAME_TO_LIGAND[name]['charge'] * count for name, count in ligand_counts.items())
    complex_charge = oxidation_state + total_ligand_charge
    
    # 4. Assemble the formula string (ligands sorted alphabetically by formula)
    ligand_formulas = []
    for name, count in sorted(ligand_counts.items(), key=lambda item: NAME_TO_LIGAND[item[0]]['formula']):
        info = NAME_TO_LIGAND[name]
        formula, count_str = info['formula'], str(count) if count > 1 else ""
        if info.get('type') == 'complex_prefix' or len(formula) > 2 or (len(formula) > 1 and formula.isalpha()):
             ligand_formulas.append(f"({formula}){count_str}")
        else:
             ligand_formulas.append(f"{formula}{count_str}")
    
    complex_formula = f"[{metal_symbol}{''.join(ligand_formulas)}]"
    return complex_formula, complex_charge

# --- CORE FUNCTION: FORMULA TO NAME ---
def name_from_formula(formula: str) -> str:
    """Generates the IUPAC name for a coordination compound from its formula."""
    match = FULL_FORMULA_RE.match(formula)
    if not match: raise ValueError(f"Could not parse formula: {formula}")
    parts = match.groupdict()

    complex_inner = parts['complex'][1:-1]
    metal_match = re.match(r"([A-Z][a-z]?)", complex_inner)
    metal_symbol = metal_match.group(1)
    ligands_str = complex_inner[len(metal_symbol):]

    ligand_counts = defaultdict(int)
    sorted_ligand_keys = sorted(LIGAND_DATA.keys(), key=len, reverse=True)
    i = 0
    while i < len(ligands_str):
        found_ligand = False
        paren_match = re.match(r"\((\w+)\)(\d*)", ligands_str[i:])
        if paren_match and paren_match.group(1) in sorted_ligand_keys:
            key, count_str = paren_match.groups()
            ligand_counts[key] += int(count_str) if count_str else 1
            i += paren_match.end()
            found_ligand = True
        if not found_ligand:
            for key in sorted_ligand_keys:
                if ligands_str[i:].startswith(key):
                    i += len(key)
                    count_match = re.match(r"(\d+)", ligands_str[i:])
                    count = 1
                    if count_match:
                        count = int(count_match.group(1))
                        i += len(count_match.group(1))
                    ligand_counts[key] += count
                    found_ligand = True
                    break
        if not found_ligand: raise ValueError(f"Unknown ligand in: {ligands_str[i:]}")

    complex_charge, cation_name, anion_name = 0, None, None
    if parts['cation']:
        cation_formula, cation_count = _parse_ion_formula(parts['cation'])
        cation_name = COUNTER_ION_DATA[cation_formula]['name']
        complex_charge = -(COUNTER_ION_DATA[cation_formula]['charge'] * cation_count)
    if parts['anion_charge']:
        charge_str = parts['anion_charge']
        if charge_str.endswith(('+', '-')):
            val_str = charge_str[:-1]
            if charge_str.endswith('+'): complex_charge = int(val_str) if val_str else 1
            else: complex_charge = -int(val_str) if val_str else -1
        else:
            anion_formula, anion_count = _parse_ion_formula(parts['anion_charge'])
            anion_name = COUNTER_ION_DATA[anion_formula]['name']
            complex_charge = -(COUNTER_ION_DATA[anion_formula]['charge'] * anion_count)

    total_ligand_charge = sum(LIGAND_DATA[lig]['charge'] * count for lig, count in ligand_counts.items())
    oxidation_state = complex_charge - total_ligand_charge

    ligand_name_parts = []
    for formula, count in sorted(ligand_counts.items(), key=lambda item: LIGAND_DATA[item[0]]['name']):
        info = LIGAND_DATA[formula]
        prefix = COMPLEX_PREFIXES.get(count, '') if info['type'] == 'complex_prefix' else PREFIXES.get(count, '')
        name_part = f"({info['name']})" if info['type'] == 'complex_prefix' else info['name']
        ligand_name_parts.append(f"{prefix}{name_part}")
    ligand_part = "".join(ligand_name_parts)

    is_anionic = complex_charge < 0 or bool(parts['cation'])
    metal_info = METAL_DATA.get(metal_symbol, {})
    metal_name = metal_info.get('anion_name', metal_symbol.lower() + 'ate') if is_anionic else metal_info.get('name', metal_symbol.lower())
    metal_part = f"{metal_name}{_to_roman(oxidation_state)}"
    complex_name = f"{ligand_part}{metal_part}"

    if cation_name: return f"{cation_name} {complex_name}"
    if anion_name: return f"{complex_name} {anion_name}"
    return f"{complex_name} ion"

# --- CORE FUNCTION: NAME TO FORMULA ---
def formula_from_name(name: str) -> str:
    """Generates the chemical formula for a coordination compound from its IUPAC name."""
    name = name.strip()
    ox_state_matches = re.findall(r"\((?:[IVX]+|0)\)", name)
    
    # Case 1: Two complex ions (e.g., complex cation + complex anion)
    if len(ox_state_matches) == 2:
        split_point = name.rfind(ox_state_matches[0]) + len(ox_state_matches[0])
        cation_name, anion_name = name[:split_point], name[split_point:]
        
        cation_formula, cation_charge = _parse_complex_name(cation_name)
        anion_formula, anion_charge = _parse_complex_name(anion_name)
        
        if cation_charge <= 0 or anion_charge >= 0:
            raise ValueError("Name must represent a complex cation followed by a complex anion.")
            
        common_divisor = math.gcd(cation_charge, abs(anion_charge))
        cation_count = abs(anion_charge) // common_divisor
        anion_count = cation_charge // common_divisor
        
        cation_count_str = str(cation_count) if cation_count > 1 else ""
        anion_count_str = str(anion_count) if anion_count > 1 else ""
        return f"{cation_formula}{cation_count_str}{anion_formula}{anion_count_str}"

    # Case 2: One complex ion (with optional counter ions or as a standalone ion)
    elif len(ox_state_matches) == 1:
        words = name.replace(" ion", "").split()
        cation_name, anion_name, complex_part_str = None, None, name
        
        if words[0] in NAME_TO_COUNTER_ION and NAME_TO_COUNTER_ION[words[0]]['charge'] > 0:
            cation_name, complex_part_str = words[0], " ".join(words[1:])
        elif words[-1] in NAME_TO_COUNTER_ION and NAME_TO_COUNTER_ION[words[-1]]['charge'] < 0:
            anion_name, complex_part_str = words[-1], " ".join(words[:-1])

        complex_formula, complex_charge = _parse_complex_name(complex_part_str.replace(" ion", ""))

        if cation_name:
            info = NAME_TO_COUNTER_ION[cation_name]
            if complex_charge >= 0: raise ValueError("Charge mismatch for counter-cation.")
            num = abs(complex_charge) // info['charge']
            count = str(num) if num > 1 else ""
            formula_str = f"({info['formula']})" if len(info['formula']) > 2 else info['formula']
            return f"{formula_str}{count}{complex_formula}"
        elif anion_name:
            info = NAME_TO_COUNTER_ION[anion_name]
            if complex_charge <= 0: raise ValueError("Charge mismatch for counter-anion.")
            num = complex_charge // abs(info['charge'])
            count = str(num) if num > 1 else ""
            formula_str = f"({info['formula']})" if len(info['formula']) > 2 else info['formula']
            return f"{complex_formula}{formula_str}{count}"
        else:
            if complex_charge == 0: return complex_formula
            sign = '+' if complex_charge > 0 else '-'
            val = abs(complex_charge)
            charge_str = f"{val}{sign}" if val > 1 else sign
            return f"{complex_formula}{charge_str}"
    else:
        raise ValueError("Could not determine compound structure from name (invalid oxidation states).")

# --- Example Usage ---
if __name__ == '__main__':
    # --- Test cases for name_from_formula ---
    print("--- Testing Formula to Name ---")
    formulas_to_test = [
        "K2[TiF6]",
        "(NH4)2[Ni(CN)4]",
        "[Fe(NH3)4Cl2]+",
        "[Co(NH3)6]Cl3",
        "Na3[Co(NO2)6]",
        "[Cr(en)3]Cl3",
    ]
    for formula in formulas_to_test:
        try:
            name = name_from_formula(formula)
            print(f"{formula:<20} -> {name}")
        except (ValueError, KeyError) as e:
            print(f"{formula:<20} -> ERROR: {e}")

    # --- Test cases for formula_from_name ---
    print("\n--- Testing Name to Formula ---")
    names_to_test = [
        "potassium hexafluoridotitanate(IV)",
        "tetraamminedichlorochromium(III) ion",
        "tetraamminedichlorochromium(III)", # Test case that was failing
        "tris(ethylenediamine)cobalt(III) chloride",
        # New test case for double complex salt
        "tetraamminedichlorochromium(III)tetrabromomanganate(II)"
    ]
    for name in names_to_test:
        try:
            formula = formula_from_name(name)
            print(f"{name:<60} -> {formula}")
        except (ValueError, KeyError) as e:
            print(f"{name:<60} -> ERROR: {e}")

def get_oxidation_state_from_formula(formula: str) -> int:
    """Calculates the oxidation state of the central metal in a coordination compound."""
    match = FULL_FORMULA_RE.match(formula)
    if not match: raise ValueError(f"Could not parse formula: {formula}")
    parts = match.groupdict()

    complex_inner = parts['complex'][1:-1]
    ligands_str = re.split(r"([A-Z][a-z]?)", complex_inner, maxsplit=1)[2]

    ligand_counts = defaultdict(int)
    sorted_ligand_keys = sorted(LIGAND_DATA.keys(), key=len, reverse=True)
    i = 0
    while i < len(ligands_str):
        found_ligand = False
        paren_match = re.match(r"\((\w+)\)(\d*)", ligands_str[i:])
        if paren_match and paren_match.group(1) in sorted_ligand_keys:
            key, count_str = paren_match.groups()
            ligand_counts[key] += int(count_str) if count_str else 1
            i += paren_match.end()
            found_ligand = True
        if not found_ligand:
            for key in sorted_ligand_keys:
                if ligands_str[i:].startswith(key):
                    i += len(key)
                    count_match = re.match(r"(\d+)", ligands_str[i:])
                    count = int(count_match.group(1)) if count_match else 1
                    i += len(count_match.group(1) if count_match else "")
                    ligand_counts[key] += count
                    found_ligand = True
                    break
        if not found_ligand: raise ValueError(f"Unknown ligand in: {ligands_str[i:]}")

    complex_charge = 0
    if parts['cation']:
        cation_formula, cation_count = _parse_ion_formula(parts['cation'])
        complex_charge = -(COUNTER_ION_DATA[cation_formula]['charge'] * cation_count)
    if parts['anion_charge']:
        charge_str = parts['anion_charge']
        if charge_str.endswith(('+', '-')):
            val_str = charge_str[:-1]
            if charge_str.endswith('+'): complex_charge = int(val_str) if val_str else 1
            else: complex_charge = -int(val_str) if val_str else -1
        else:
            anion_formula, anion_count = _parse_ion_formula(parts['anion_charge'])
            complex_charge = -(COUNTER_ION_DATA[anion_formula]['charge'] * anion_count)

    total_ligand_charge = sum(LIGAND_DATA[lig]['charge'] * count for lig, count in ligand_counts.items())
    return complex_charge - total_ligand_charge

def compare_ligand_field_splitting(complexes: list[str]) -> str:
    """Ranks complexes based on ligand field splitting energy (Δ), considering geometry."""
    ligand_strength = {ligand: i for i, ligand in enumerate(SPECTROCHEMICAL_SERIES)}
    complex_scores = {}

    for formula in complexes:
        try:
            match = FULL_FORMULA_RE.match(formula)
            if not match:
                complex_scores[formula] = -1
                continue

            complex_inner = match.group('complex')[1:-1]
            ligands_str = re.split(r"([A-Z][a-z]?)", complex_inner, maxsplit=1)[2]
            
            ligand_counts = defaultdict(int)
            sorted_ligand_keys = sorted(LIGAND_DATA.keys(), key=len, reverse=True)
            i = 0
            while i < len(ligands_str):
                found_ligand = False
                paren_match = re.match(r"\((\w+)\)(\d*)", ligands_str[i:])
                if paren_match and paren_match.group(1) in sorted_ligand_keys:
                    key, count_str = paren_match.groups()
                    ligand_counts[key] += int(count_str) if count_str else 1
                    i += paren_match.end()
                    found_ligand = True
                if not found_ligand:
                    for key in sorted_ligand_keys:
                        if ligands_str[i:].startswith(key):
                            i += len(key)
                            count_match = re.match(r"(\d+)", ligands_str[i:])
                            count = int(count_match.group(1)) if count_match else 1
                            i += len(count_match.group(1) if count_match else "")
                            ligand_counts[key] += count
                            found_ligand = True
                            break
                if not found_ligand:
                    i += 1
            
            if not ligand_counts:
                complex_scores[formula] = -1
                continue

            # Calculate the average ligand strength score
            total_score = sum(ligand_strength.get(lig, -1) * count for lig, count in ligand_counts.items())
            coordination_number = sum(ligand_counts.values())
            avg_score = total_score / coordination_number if coordination_number > 0 else -1

            # Apply a penalty for tetrahedral geometry (CN=4)
            if coordination_number == 4:
                avg_score *= (4/9) # Approximate reduction for tetrahedral vs octahedral
            
            complex_scores[formula] = avg_score

        except (ValueError, IndexError):
            complex_scores[formula] = -1 # Assign a low score if parsing fails

    sorted_complexes = sorted(complex_scores, key=complex_scores.get)
    return " < ".join(sorted_complexes)