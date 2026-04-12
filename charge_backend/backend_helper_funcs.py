from loguru import logger
from fastapi import WebSocket
import asyncio
import json
from pydantic.dataclasses import dataclass
from pydantic import Field
from typing import Any, Dict, Optional, Literal, Tuple
from dataclasses import asdict

from flask_tools.chemistry import smiles_utils
from flask_tools.lmo.molecular_property_utils import get_density
from lc_conductor import CallbackLogger, RunSettings
from charge_backend.moleculedb.molecule_naming import MolNameFormat


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
    x: Optional[float] = None
    y: Optional[float] = None
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


@dataclass
class FlaskRunSettings(RunSettings):
    molecule_name_format: MolNameFormat = Field(alias="moleculeName", default="brand")
    # AI-based retrosynthesis flag
    use_ai_based: bool = Field(alias="useAiBased", default=True)
    # RSA settings
    use_rsa: bool = Field(alias="useRsa", default=False)
    rsa_n: int = Field(alias="rsaN", default=8)
    rsa_k: int = Field(alias="rsaK", default=4)
    rsa_t: int = Field(alias="rsaT", default=3)
    rsa_mode: str = Field(alias="rsaMode", default="standalone")


@dataclass(frozen=True)
class RdkitjsMolPayload:
    """Best-practice payload for rdkit.js rendering.

    - Provide highlight atom indices for direct rdkit.js highlighting.
    - Also include 1-based indices (highlight_atom_mapnums) as a stable debugging aid.
      These are NOT atom-map numbers.
    """

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
        self.pending_tool_calls: dict[str, dict[str, Any]] = {}
        self.pending_tool_calls_by_name: dict[str, list[dict[str, Any]]] = {}
        self._pending_send_tasks: set[asyncio.Task[Any]] = set()

    @staticmethod
    def _parse_tool_arguments(
        arguments: Any,
    ) -> tuple[str, Optional[str], Optional[str]]:
        if isinstance(arguments, str):
            arguments_str = arguments
            try:
                parsed = json.loads(arguments)
            except Exception:
                parsed = None
        elif arguments is None:
            arguments_str = ""
            parsed = None
        else:
            parsed = arguments
            arguments_str = json.dumps(arguments)

        log_msg = None
        smiles = None
        if isinstance(parsed, dict):
            raw_log_msg = parsed.get("log_msg")
            raw_smiles = parsed.get("smiles") or parsed.get("product_smiles")
            log_msg = raw_log_msg if isinstance(raw_log_msg, str) else None
            smiles = raw_smiles if isinstance(raw_smiles, str) else None
        return arguments_str, log_msg, smiles

    @staticmethod
    def _merge_arguments_str(existing: str, new: str) -> str:
        if not existing:
            return new
        if not new or new == existing:
            return existing
        return existing + new

    @staticmethod
    def _try_parse_json(value: Any) -> Any:
        if isinstance(value, str):
            try:
                return json.loads(value)
            except Exception:
                return None
        return value

    @classmethod
    def _format_arguments_str(cls, arguments_str: str) -> str:
        parsed = cls._try_parse_json(arguments_str)
        if parsed is None:
            return arguments_str
        if isinstance(parsed, dict):
            lines = []
            for key, value in parsed.items():
                if isinstance(value, (dict, list)):
                    formatted_value = json.dumps(
                        value, indent=2, sort_keys=True, default=str
                    )
                else:
                    formatted_value = str(value)
                lines.append(f"{key}: {formatted_value}")
            return "\n".join(lines)
        if isinstance(parsed, list):
            return "\n".join(str(item) for item in parsed)
        return str(parsed)

    @classmethod
    def _format_result_text(cls, result: Any) -> tuple[str, str]:
        parsed = cls._try_parse_json(result)
        if parsed is not None:
            return "json", json.dumps(parsed, indent=2, sort_keys=True, default=str)

        if isinstance(result, str):
            return "", result
        return "json", json.dumps(result, indent=2, sort_keys=True, default=str)

    async def on_tool_call(
        self,
        tool_name: str,
        arguments: Any = None,
        *,
        source: Optional[str] = None,
        call_id: Optional[str] = None,
    ) -> None:
        source_name = source or self.name or "FLASK-Generic"
        arguments_str, log_msg, smiles = self._parse_tool_arguments(arguments)
        call_info: dict[str, Any] = {
            "tool_name": tool_name,
            "source_name": source_name,
            "arguments_str": arguments_str,
            "log_msg": log_msg,
            "smiles": smiles,
        }
        if call_id is not None:
            existing_call_info = self.pending_tool_calls.get(call_id)
            if existing_call_info is not None:
                existing_call_info["tool_name"] = tool_name
                existing_call_info["source_name"] = source_name
                existing_call_info["arguments_str"] = self._merge_arguments_str(
                    existing_call_info.get("arguments_str", ""),
                    arguments_str,
                )
                if log_msg is not None:
                    existing_call_info["log_msg"] = log_msg
                if smiles is not None:
                    existing_call_info["smiles"] = smiles
            else:
                self.pending_tool_calls[call_id] = call_info
        else:
            self.pending_tool_calls_by_name.setdefault(tool_name, []).append(call_info)

    def _pop_pending_tool_call(
        self,
        tool_name: str,
        call_id: Optional[str] = None,
    ) -> Optional[dict[str, Any]]:
        if call_id is not None:
            return self.pending_tool_calls.pop(call_id, None)

        pending_calls = self.pending_tool_calls_by_name.get(tool_name)
        if not pending_calls:
            return None

        call_info = pending_calls.pop(0)
        if not pending_calls:
            self.pending_tool_calls_by_name.pop(tool_name, None)
        return call_info

    async def on_tool_result(
        self,
        tool_name: str,
        result: Any,
        *,
        is_error: bool = False,
        source: Optional[str] = None,
        call_id: Optional[str] = None,
    ) -> None:
        send = self.websocket.send_json
        call_info = self._pop_pending_tool_call(tool_name, call_id)
        source_name = (
            call_info["source_name"]
            if call_info is not None
            else source or self.name or "FLASK-Generic"
        )
        result_lang, formatted_result_text = self._format_result_text(result)
        fence_suffix = result_lang if result_lang else ""
        if call_info is not None and call_info["log_msg"]:
            message = (
                f"{call_info['log_msg']}\nArguments:\n```\n"
                f"{self._format_arguments_str(call_info['arguments_str'])}\n```\n"
                f"Result:\n```{fence_suffix}\n{formatted_result_text}\n```"
            )
        elif call_info is not None:
            message = (
                f"[{source_name}] `{tool_name}` called with:\n```\n"
                f"{self._format_arguments_str(call_info['arguments_str'])}\n```\n"
                f"Returned:\n```{fence_suffix}\n{formatted_result_text}\n```"
            )
        elif is_error:
            message = (
                f"[{source_name}] `{tool_name}` errored:\n"
                f"```{fence_suffix}\n{formatted_result_text}\n```"
            )
        else:
            message = (
                f"[{source_name}] `{tool_name}` returned:\n"
                f"```{fence_suffix}\n{formatted_result_text}\n```"
            )

        payload: dict[str, Any] = {
            "type": "response",
            "message": {"source": tool_name, "message": message},
        }
        if call_info is not None and call_info["smiles"] is not None:
            payload["smiles"] = call_info["smiles"]
        if call_id is not None:
            payload["toolCallId"] = call_id

        logger.info(message)
        await send(payload)

    async def send(self, assistant_message):
        send = self.websocket.send_json
        source = self.name if self.name else "FLASK-Generic"
        try:
            message_type = getattr(assistant_message, "type", None)
            if message_type == "UserMessage":
                message = (
                    f"[{source}] User: {getattr(assistant_message, 'content', '')}"
                )
                logger.info(message)
                await send({"type": "response", "message": {"message": message}})
            elif message_type == "AssistantMessage":
                thought = getattr(assistant_message, "thought", None)
                if thought is not None:
                    _str = f"[{source}] Model thought: {thought}"
                    output = ModelMessage(message=_str, smiles=None)
                    logger.info(_str)
                    await send({"type": "response", "message": output.json()})
                content = getattr(assistant_message, "content", None)
                if isinstance(content, list):
                    for item in content:
                        if hasattr(item, "name") and hasattr(item, "arguments"):
                            await self.on_tool_call(
                                item.name,
                                item.arguments,
                                source=source,
                                call_id=getattr(item, "id", None),
                            )
                        else:
                            logger.info(f"[{source}] Model: {item}")
            elif message_type == "FunctionExecutionResultMessage":
                content = getattr(assistant_message, "content", None) or []
                for result in content:
                    if result.is_error:
                        await self.on_tool_result(
                            result.name,
                            result.content,
                            is_error=True,
                            source=source,
                            call_id=getattr(result, "call_id", None),
                        )
                    else:
                        await self.on_tool_result(
                            result.name,
                            result.content,
                            source=source,
                            call_id=getattr(result, "call_id", None),
                        )
            else:
                nested_message = getattr(assistant_message, "message", None)
                nested_content = getattr(nested_message, "content", assistant_message)
                message = f"[{source}] Model: {nested_content}"
                await self.clogger.info(message)
        except Exception:
            logger.exception("CallbackHandler failed while sending assistant message")
            raise

    def __call__(self, assistant_message):
        task = asyncio.create_task(self.send(assistant_message))
        self._pending_send_tasks.add(task)
        task.add_done_callback(self._pending_send_tasks.discard)
        return task

    async def drain(self) -> None:
        if not self._pending_send_tasks:
            return
        results = await asyncio.gather(
            *list(self._pending_send_tasks), return_exceptions=True
        )
        for result in results:
            if isinstance(result, Exception):
                logger.error(
                    f"CallbackHandler background send task failed: {type(result).__name__}: {result}"
                )


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
    # For properties where there are "fast" approximate values (from RDKit), estimate them
    canonical_smiles = smiles_utils.canonicalize_smiles(smiles)
    density = float(tool_properties.get("density", get_density(canonical_smiles)[1]))
    sascore = float(
        tool_properties.get(
            "synthesizability", smiles_utils.get_synthesizability(canonical_smiles)
        )
    )
    bandgap = float(tool_properties.get("bandgap", get_bandgap(canonical_smiles)))
    lmo_props = {
        "smiles": canonical_smiles,
        "parent_id": parent_id,
        "node_id": node_id,
        "density": density,
        "sascore": sascore,
        "bandgap": bandgap,
    }
    for k, v in tool_properties.items():
        if k not in lmo_props:
            lmo_props[k] = v
    return lmo_props


def post_process_smiles(smiles: str, parent_id: int, node_id: int) -> dict:
    """
    Post-process a solution SMILES string, add additional properties and return
    a dictionary that can be appended to the known molecules JSON file.

    Args:
        smiles (str): The input SMILES string.
    Returns:
        dict: The post-processed dictionary.
    """
    canonical_smiles = smiles_utils.canonicalize_smiles(smiles)
    sascore = smiles_utils.get_synthesizability(canonical_smiles)
    density = get_density(canonical_smiles)

    return {"smiles": canonical_smiles, "sascore": sascore, "density": density}
