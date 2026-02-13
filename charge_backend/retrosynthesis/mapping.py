from __future__ import annotations

from typing import Any, Dict, Optional, Sequence

from loguru import logger

from rdkitjs_payload import build_rdkitjs_mapped_reaction, mapped_reaction_to_json_dict


def build_mapped_reaction_dict_or_none(
    *,
    reactants: Sequence[str],
    products: Sequence[str],
    log_msg: str,
    **log_kwargs: Any,
) -> Optional[Dict[str, Any]]:
    """Best-effort rdkit.js mapping builder.

    Returns a JSON-serializable mapping dict on success, otherwise None.
    """

    if not reactants or not products:
        return None

    try:
        mapped_reaction = build_rdkitjs_mapped_reaction(
            reactants=list(reactants),
            products=list(products),
        )
        return mapped_reaction_to_json_dict(mapped_reaction)
    except Exception as e:
        logger.opt(exception=e).warning(log_msg, **log_kwargs)
        return None
