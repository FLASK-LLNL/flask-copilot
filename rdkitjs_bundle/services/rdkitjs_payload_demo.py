from __future__ import annotations

import json
import os
import sys


# Allow `python3 services/rdkitjs_payload_demo.py` from repo root.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from services.rdkitjs_payload import build_rdkitjs_reaction_payload, reaction_payload_to_json_dict


def main() -> None:
    # Example from earlier discussion
    reactants = ["CN(N)c1nnnn1C", "BrBr"]
    products = ["CN(/C=N/N(C)c1nnnn1C)c1nnnn1C"]

    payload = build_rdkitjs_reaction_payload(
        reactants=reactants,
        products=products,
        highlight_rgb=(220, 0, 0),
        highlight_alpha=0.45,
    )
    print(json.dumps(reaction_payload_to_json_dict(payload), indent=2))


if __name__ == "__main__":
    main()
