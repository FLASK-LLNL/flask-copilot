from loguru import logger
from fastapi import WebSocket
import asyncio
import json
from typing import Any, Dict, Optional, Literal, Tuple
from dataclasses import dataclass, asdict
from collections import defaultdict

from charge.clients.autogen import AutoGenAgent
import charge.servers.AiZynthTools as aizynth_funcs
from charge.servers import SMILES_utils
from charge.servers.molecular_property_utils import get_density
from callback_logger import CallbackLogger


@dataclass
class PathwayStep:
    smiles: list[str]
    label: list[str]
    parents: list[int]

    def json(self):
        return asdict(self)


@dataclass
class ReactionAlternative:
    id: str
    name: str
    type: Literal["exact", "template", "ai"]
    status: Literal["active", "available", "computing"]
    pathway: list[PathwayStep]
    hoverInfo: str
    disabled: Optional[bool] = None
    disabledReason: Optional[str] = None

    def json(self):
        return asdict(self)


@dataclass
class Reaction:
    id: str
    hoverInfo: str
    highlight: str = "normal"
    label: Optional[str] = None
    alternatives: Optional[list[ReactionAlternative]] = None
    templatesSearched: bool = False
    mappedReaction: Optional[Dict[str, Any]] = None

    def json(self):
        return asdict(self)


@dataclass
class Node:
    id: str
    smiles: str
    label: str
    hoverInfo: str
    level: int
    parentId: Optional[str] = None
    x: Optional[int] = None
    y: Optional[int] = None
    highlight: Optional[str] = "normal"
    # Properties
    cost: Optional[float] = None
    bandgap: Optional[float] = None
    yield_: Optional[float] = None
    density: Optional[float] = None
    sascore: Optional[float] = None
    purchasable: Optional[bool] = None
    # Reaction properties
    reaction: Optional[Reaction] = None

    def json(self):
        ret = asdict(self)
        ret["yield"] = ret["yield_"]
        del ret["yield_"]
        return ret


@dataclass
class Edge:
    id: str
    fromNode: str
    toNode: str
    status: Literal["computing", "complete"]
    label: Optional[str] = None

    def json(self):
        return asdict(self)


@dataclass
class ModelMessage:
    message: str
    smiles: Optional[str]

    def json(self):
        ret = asdict(self)
        if self.smiles is None:
            if "smiles" in ret:
                del ret["smiles"]
        return ret


@dataclass(frozen=True)
class RdkitjsMolPayload:
    """Best-practice payload for rdkit.js rendering.

    - Use a MolBlock to preserve atom ordering across Python -> JS.
    - Provide highlight atom indices for direct rdkit.js highlighting.
    - Also include 1-based indices (highlight_atom_mapnums) as a stable debugging aid.
      These are NOT atom-map numbers.
    """

    molblock: str
    smiles: str
    highlight_atom_idxs: list[int]
    highlight_atom_mapnums: list[int]


def get_price(smiles: str) -> float:
    """Mock function to get price of a molecule given its SMILES string."""
    # In a real implementation, this would query a database or an API.
    # Here, we return a random price for demonstration purposes.
    import random

    return round(random.uniform(10.0, 100.0), 2)


def get_yield(parent, child_smiles_list):
    """Mock function to get yield of a reaction given parent node and child SMILES strings."""
    # In a real implementation, this would query a database or an API.
    # Here, we return a random yield for demonstration purposes.
    import random

    return round(random.uniform(50.0, 100.0), 2)


def get_bandgap(smiles: str) -> float:
    """Mock function to get bandgap of a molecule given its SMILES string."""
    # In a real implementation, this would query a database or an API.
    # Here, we return a random bandgap for demonstration purposes.
    import random

    return round(random.uniform(1.0, 5.0), 2)


class CallbackHandler:
    def __init__(self, websocket: WebSocket, name: Optional[str] = None):
        self.websocket = websocket
        self.name = name
        self.clogger = CallbackLogger(websocket)

    async def send(self, assistant_message):
        send = self.websocket.send_json
        source = self.name if self.name else "FLASK-Generic"
        if assistant_message.type == "UserMessage":
            message = f"[{source}] User: {assistant_message.content}"
            logger.info(message)
            await send({"type": "response", "message": {"message": message}})
        elif assistant_message.type == "AssistantMessage":

            if assistant_message.thought is not None:
                _str = f"[{source}] Model thought: {assistant_message.thought}"
                output = ModelMessage(message=_str, smiles=None)
                logger.info(_str)
                await send({"type": "response", "message": output.json()})
            if isinstance(assistant_message.content, list):
                for item in assistant_message.content:
                    if hasattr(item, "name") and hasattr(item, "arguments"):
                        name = item.name
                        _str = f"[{source}] Calling {item.name}: with args {item.arguments}"
                        logger.info(_str)
                        if "log_msg" in item.arguments:
                            str_to_dict = json.loads(item.arguments)
                            if "log_msg" in str_to_dict:
                                _str = str_to_dict["log_msg"]

                        msg = {
                            "type": "response",
                            "message": {"source": name, "message": _str},
                        }
                        if "smiles" in item.arguments:
                            str_to_dict = json.loads(item.arguments)
                            if "smiles" in str_to_dict:
                                msg["smiles"] = str_to_dict["smiles"]
                        await send(msg)
                    else:
                        logger.info(f"[{source}] Model: {item}")
        elif assistant_message.type == "FunctionExecutionResultMessage":

            for result in assistant_message.content:
                if result.is_error:
                    await self.clogger.error(
                        f"[{source}] Function {result.name} errored with output: {result.content}",
                        source=result.name,
                    )
                else:
                    await self.clogger.info(
                        f"[{source}] Returning {result.name}: {result.content}",
                        source=result.name,
                    )
        else:
            message = f"[{source}] Model: {assistant_message.message.content}"
            await self.clogger.info(message)

    def __call__(self, assistant_message):
        asyncio.create_task(self.send(assistant_message))


def calculate_positions(nodes: list[Node], y_offset: int = 0):
    """
    Calculate positions for all nodes (matching frontend logic).
    Operates in-place.
    """
    BOX_WIDTH = 270  # Must match with javascript!
    BOX_GAP = 160  # Must match with javascript!
    level_gap = BOX_WIDTH + BOX_GAP
    node_spacing = 150

    # Group by level
    levels = {}
    for node in nodes:
        level = node.level
        if level not in levels:
            levels[level] = []
        levels[level].append(node)

    # Position nodes
    for node in nodes:
        level_nodes = levels[node.level]
        index_in_level = level_nodes.index(node) + y_offset

        node.x = 100 + node.level * level_gap
        node.y = 100 + index_in_level * node_spacing


async def highlight_node(node: Node, websocket: WebSocket, highlight: bool):
    await websocket.send_json(
        {
            "type": "node_update",
            "node": {
                "id": node.id,
                "highlight": "yellow" if highlight else "normal",
            },
        }
    )


async def loop_executor(executor, func, *args, **kwargs):
    loop = asyncio.get_event_loop()
    return await loop.run_in_executor(executor, func, *args, **kwargs)


def post_process_lmo_smiles(
    smiles: str, parent_id: int, node_id: int, tool_properties: Optional[dict] = None
) -> Dict:
    """
    Post-process LMO SMILES, preferring tool-calculated properties.

    Args:
        tool_properties: If provided, properties are taken from here
                        (avoids recalculation)
    """
    tool_properties = tool_properties or {}
    canonical_smiles = SMILES_utils.canonicalize_smiles(smiles)
    density = float(tool_properties.get("density", get_density(canonical_smiles)))
    sascore = float(
        tool_properties.get(
            "synthesizability", SMILES_utils.get_synthesizability(canonical_smiles)
        )
    )
    bandgap = float(tool_properties.get("bandgap", get_bandgap(canonical_smiles)))
    return {
        "smiles": canonical_smiles,
        "parent_id": parent_id,
        "node_id": node_id,
        "density": density,
        "sascore": sascore,
        "bandgap": bandgap,
    }
