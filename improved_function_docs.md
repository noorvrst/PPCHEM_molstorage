# MolStorage Package Documentation

## Overview
Our package `molstorage` contains 8 main functions for chemical compound management, safety analysis, and storage organization.

---

### 1. `get_compound_safety_data`

**Purpose**: Retrieves safety data for a chemical compound from PubChem.

**Input**:
- `name: str` - The molecule's common name

**Output**: 
- `Tuple[str, List[str], List[str]]` - Contains:
  - PubChem CID (ID)
  - GHS pictograms 
  - Hazard statements

**Description**: Connects to PubChem's REST API to fetch safety information for the specified compound.

---

### 2. `get_name_and_smiles`

**Purpose**: Retrieves identification information for a chemical compound.

**Input**:
- `cid: str` - PubChem compound ID (CID)

**Output**:
- `Tuple[str, str, str]` - Contains:
  - Record Title (generic name)
  - IUPAC name
  - SMILES string

**Description**: Uses the PubChemPy library to retrieve standardized chemical identifiers.

---

### 3. `classify_acid_base`

**Purpose**: Determines the acid/base classification of a compound.

**Input**:
- `name: str` - The generic name of the compound
- `iupac_name: str` - The IUPAC name of the compound
- `smiles: str` - SMILES string representing the compound's chemical structure
- `ghs_statements: List[str]` - List of GHS hazard statements related to the compound

**Output**:
- `Union[str, Tuple[str, ...]]` - Classification such as "acid", "base", "neutral", or "amphoteric"

**Description**: Analyzes compound properties to determine its acid/base behavior for compatibility assessment.

---

### 4. `get_mp_bp`

**Purpose**: Retrieves melting and boiling points for a compound.

**Input**:
- `name: str` - The compound's name

**Output**:
- `Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]` - Contains:
  - Melting point (°C)
  - Boiling point (°C)
  - Melting point (°F)
  - Boiling point (°F)

**Description**: Extracts temperature data from PubChem, returning `None` for any values that cannot be found.

---

### 5. `compound_state`

**Purpose**: Predicts the physical state of a compound at room temperature.

**Input**:
- `mp_c: Optional[float]` - Melting point in Celsius
- `bp_c: Optional[float]` - Boiling point in Celsius
- `mp_f: Optional[float]` - Melting point in Fahrenheit
- `bp_f: Optional[float]` - Boiling point in Fahrenheit

**Output**:
- `str` - Physical state ("solid", "liquid", "gas", or "unknown")

**Description**: Determines compound state at room temperature (20°C / 68°F) based on melting and boiling points.

---

### 6. `prioritize_pictograms`

**Purpose**: Sorts GHS pictograms by hazard severity.

**Input**:
- `pictograms: List[str]` - List of GHS pictogram names

**Output**:
- `List[str]` - Sorted list of pictograms

**Description**: Orders pictograms according to predefined hazard severity priority (lower numbers = more severe hazards).

---

### 7. `is_chemically_compatible`

**Purpose**: Determines whether two chemicals can be safely stored together.

**Input**:
- `existing_pictograms: List[str]` - GHS pictograms for the existing chemical
- `new_pictograms: List[str]` - GHS pictograms for the new chemical
- `existing_acid_base_class: str` - Acid/base classification of the existing chemical
- `new_acid_base_class: str` - Acid/base classification of the new chemical
- `existing_state: str` - Physical state of the existing chemical
- `new_state: str` - Physical state of the new chemical
- `group_name: str` - Storage group used for applying group-specific compatibility rules

**Output**:
- `bool` - `True` if compatible, `False` otherwise

**Description**: Evaluates compatibility based on pictograms, acid/base classifications, physical states, and storage group.

---

### 8. `chemsort_multiple_order_3`

**Purpose**: Sorts compounds into compatible storage groups.

**Input**:
- `compounds: List[Dict[str, Any]]` - List of compounds to sort, each represented as a dictionary with keys:
  - `name`
  - `sorted_pictograms`
  - `hazard_statements`
  - `acid_base_class`
  - `state_room_temp`
- `storage_groups: Dict[str, Dict[str, List[Dict[str, Any]]]]` - Existing storage groups dictionary

**Output**:
- `Dict[str, Dict[str, List[Dict[str, Any]]]]` - Updated dictionary of storage groups with sorted compounds

**Description**: Processes compounds by assigning them to predefined groups if compatible, or creates new custom groups as needed. Each storage group contains sub-categories for solids, liquids, and gases.
