################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from dataclasses import dataclass, asdict
from typing import Optional, Literal


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
    density: Optional[float] = None
    sascore: Optional[float] = None

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
