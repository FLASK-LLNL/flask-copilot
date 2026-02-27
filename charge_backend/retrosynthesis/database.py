from dataclasses import dataclass
from fastapi import WebSocket
import os
from typing import Any, Callable, Literal

from charge_backend.retrosynthesis.context import RetrosynthesisContext
from lc_conductor.callback_logger import CallbackLogger
from backend_helper_funcs import (
    Node,
    Reaction,
    ReactionAlternative,
    PathwayStep,
    FlaskRunSettings,
)
from charge_backend.moleculedb.molecule_naming import (
    smiles_to_html,
)
from charge_backend.moleculedb.purchasable import is_purchasable
from charge_backend.moleculedb.dynamic_import import import_from_path
from charge_backend.moleculedb.reactiondb_query import ReactionDatabaseReader

try:
    from rdkit import Chem
except ImportError:
    Chem = None


@dataclass
class ReactionDBEntry:
    """
    A parsed entry in the reaction database
    """

    name: str
    text: str
    components: list[dict[Literal["role", "name", "smiles", "inchi"], str]]
    reaction_yield: float = -1.0


def generate_hover_info(entry: ReactionDBEntry) -> str:
    return f"""# {entry.name}
* Yield: {"Not reported" if entry.reaction_yield < 0 else entry.reaction_yield}
## Description
{entry.text or "Not reported"}
"""


ReactionParserFunction = Callable[[str, dict[str, Any]], ReactionDBEntry]


REACTIONDB_PATH = os.getenv("FLASK_REACTION_DB", "/data/db/reactions.db")
REACTIONDB_PARSER_PATH = os.getenv(
    "FLASK_REACTION_DB_PARSER",
    os.path.join(os.path.dirname(REACTIONDB_PATH), "parse_entry.py"),
)

# Globals
REACTIONDB_HANDLE = None
db_entry_to_reaction: ReactionParserFunction | None = None


def _smiles_to_inchi(smiles: str) -> str | None:
    if Chem is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    inchi = Chem.MolToInchi(mol)
    return str(inchi)


async def find_exact_reactions(
    product: Node,
    context: RetrosynthesisContext,
    clogger: CallbackLogger,
    websocket: WebSocket,
    run_settings: FlaskRunSettings,
) -> Reaction | None:
    # Load state
    global REACTIONDB_HANDLE
    global db_entry_to_reaction
    if REACTIONDB_HANDLE is None:
        if os.path.exists(REACTIONDB_PATH):
            REACTIONDB_HANDLE = ReactionDatabaseReader(REACTIONDB_PATH)
        else:
            await clogger.warning(f"Cannot load database at {REACTIONDB_PATH}")
            return None
    if db_entry_to_reaction is None:
        if not os.path.exists(REACTIONDB_PARSER_PATH):
            await clogger.warning(
                f"Database entry parser not found at {REACTIONDB_PARSER_PATH}"
            )
            return None
        mod = import_from_path("parse_entry", REACTIONDB_PARSER_PATH)
        if mod is None:
            await clogger.error(
                f"Cannot import database entry parser at {REACTIONDB_PARSER_PATH}"
            )
            return None
        db_entry_to_reaction = mod.db_entry_to_reaction
        if db_entry_to_reaction is None:
            await clogger.error(
                f"Cannot load database entry parser function at {REACTIONDB_PARSER_PATH}"
            )
            return None

    product_inchi = _smiles_to_inchi(product.smiles)
    if product_inchi is None:
        return None

    entries = REACTIONDB_HANDLE.get(product_inchi)
    if not entries:
        return None

    # Create an output reaction
    alternatives: list[ReactionAlternative] = []
    result = Reaction(
        "dbentry", "", label="Exact", highlight="normal", alternatives=alternatives
    )

    # Process database entries
    processed = [db_entry_to_reaction(product_inchi, entry) for entry in entries]
    # Sort entries by ones having a description coming first
    processed_entries = [e for e in processed if e.text] + [
        e for e in processed if not e.text
    ]

    # Optionally filter reactions where the product appears as a reactant
    processed_entries = [
        e
        for e in processed_entries
        if not any(
            c["inchi"] == product_inchi and c["role"] != "Product"
            for c in e.components
            if "role" in c and "inchi" in c
        )
    ]

    # Sort so entries where product appears first come first
    processed_entries.sort(
        key=lambda e: not next(
            (
                c["inchi"] == product_inchi
                for c in e.components
                if "role" in c and "inchi" in c
            ),
            False,
        )
    )

    for i, processed in enumerate(processed_entries):

        # Create pathway with reagents
        reactant_smiles: list[str] = []
        reactant_labels: list[str] = []
        reactant_roles: list[str] = []
        product_name = product.label
        for component in processed.components:
            if "role" not in component:
                continue
            if component["role"] == "Product" and "name" in component:
                product_name = component["name"]
                continue
            if component["role"] != "Product" and "smiles" in component:
                if "name" in component:
                    name = component["name"] + f" ({component['role']})"
                else:
                    name = ""
                reactant_smiles.append(component["smiles"])
                reactant_labels.append(name)
                reactant_roles.append(component["role"])
        pathway = [
            PathwayStep([product.smiles], [product_name], [-1]),
            PathwayStep(
                reactant_smiles,
                [
                    label or smiles_to_html(smiles, run_settings.molecule_name_format)
                    for label, smiles in zip(reactant_labels, reactant_smiles)
                ],
                [0] * len(reactant_smiles),
            ),
        ]

        alt = ReactionAlternative(
            f"dbentry_{i}",
            processed.name,
            type="exact",
            status="available",
            pathway=pathway,
            hoverInfo=generate_hover_info(processed),
        )
        alternatives.append(alt)

        if i == 0:
            alt.status = "active"
            result.id = alt.id
            result.hoverInfo = alt.hoverInfo
            # result.label = alt.name
            if processed.reaction_yield >= 0:
                product.yield_ = processed.reaction_yield
            if product_name != product.label:
                product.hoverInfo += f"\n\n**Given name in database**: {product_name}"

            # Transmit nodes from first database entry back to UI
            for smiles, label, role in zip(
                reactant_smiles, reactant_labels, reactant_roles
            ):
                mol_sources = is_purchasable(smiles)
                if mol_sources:
                    purchasable_str = f"Yes (via {', '.join(mol_sources)})"
                else:
                    purchasable_str = "No"
                node = Node(
                    id=context.new_node_id(),
                    smiles=smiles,
                    label=label
                    or smiles_to_html(smiles, run_settings.molecule_name_format),
                    hoverInfo=f"""# {role}
{("  * Given name in database: " + label) if label else ""}
  * SMILES: {smiles}
  * Purchasable? {purchasable_str}
""",
                    level=product.level + 1,
                    parentId=product.id,
                    purchasable=(len(mol_sources) > 0),
                )
                await context.add_node(node, product, websocket)

    return result
