"""Storage of hazardous chemicals based on their pictograms, hazards statements, acid/base classification, state"""

from __future__ import annotations

__version__ = "0.0.1"

from .molstorage import get_compound_safety_data, get_name_and_smiles, classify_acid_base, get_mp_bp, compound_state, prioritize_pictograms, is_chemically_compatible, default_group, default_group_gas, initialize_storage_groups, chemsort_multiple_order_3 

from .app_run import App
