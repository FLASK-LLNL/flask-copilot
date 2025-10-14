from loguru import logger
from fastapi import WebSocket
import asyncio
import json
from typing import Dict, Optional, Literal, Tuple
from dataclasses import dataclass, asdict
from collections import defaultdict

from charge.clients.autogen import AutoGenClient
import charge.servers.AiZynthTools as aizynth_funcs


# TODO: Put this on the top level package and make it reusable
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
    # Properties
    cost: Optional[float] = None
    bandgap: Optional[float] = None
    yield_: Optional[float] = None
    highlight: Optional[str] = "normal"

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
    def __init__(self, websocket: WebSocket):
        self.websocket = websocket

    async def send(self, assistant_message):
        send = self.websocket.send_json
        if assistant_message.type == "UserMessage":
            message = f"User: {assistant_message.content}"
            logger.info(message)
            await send({"type": "response", "message": message})
        elif assistant_message.type == "AssistantMessage":

            if assistant_message.thought is not None:
                _str = f"Model thought: {assistant_message.thought}"
                output = ModelMessage(message=_str, smiles=None)
                logger.info(_str)
                await send({"type": "response", **output.json()})
            if isinstance(assistant_message.content, list):
                for item in assistant_message.content:
                    if hasattr(item, "name") and hasattr(item, "arguments"):
                        _str = f"Function call: {item.name} with args {item.arguments}"
                        logger.info(_str)
                        if "log_msg" in item.arguments:
                            str_to_dict = json.loads(item.arguments)
                            if "log_msg" in str_to_dict:
                                _str = str_to_dict["log_msg"]

                        msg = {"type": "response", "message": _str}
                        if "smiles" in item.arguments:
                            str_to_dict = json.loads(item.arguments)
                            if "smiles" in str_to_dict:
                                msg["smiles"] = str_to_dict["smiles"]
                        await send(msg)

                    else:

                        logger.info(f"Model: {item}")
        elif assistant_message.type == "FunctionExecutionResultMessage":

            for result in assistant_message.content:
                if result.is_error:
                    message = (
                        f"Function {result.name} errored with output: {result.content}"
                    )
                    logger.error(message)
                else:
                    message = f"Function {result.name} returned: {result.content}"
                    logger.info(message)
        else:
            message = f"Model: {assistant_message.message.content}"
            logger.info(message)

    def __call__(self, assistant_message):
        asyncio.create_task(self.send(assistant_message))


class RetrosynthesisContext:
    """
    Manages a retrosynthesis experiment
    """

    def __init__(self):
        self.node_ids: dict[str, Node] = {}
        self.node_id_to_planner: dict[str, aizynth_funcs.RetroPlanner] = {}
        self.node_id_to_charge_client: dict[str, AutoGenClient] = {}
        self.azf_nodes: dict[str, aizynth_funcs.Node] = {}
        self.nodes_per_level: dict[int, int] = defaultdict(int)
        self.parents: dict[str, str] = {}

    def reset(self):
        self.node_ids.clear()
        self.node_id_to_planner.clear()
        self.node_id_to_charge_client.clear()
        self.azf_nodes.clear()
        self.nodes_per_level.clear()
        self.parents.clear()


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
            "id": node.id,
            "highlight": "yellow" if highlight else "normal",
        }
    )
