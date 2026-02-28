from dataclasses import dataclass
from typing import List, Optional
from flask_tools.chemistry.smiles_utils import verify_smiles
from charge.tasks.task import Task

try:
    from rdkit import Chem
except ModuleNotFoundError:
    Chem = None


def _smiles_to_inchi(smiles: str) -> str:
    if Chem is None:
        return smiles
    mol = Chem.MolFromSmiles(smiles)
    return str(Chem.MolToInchi(mol))


@dataclass
class ReactionTask(Task):
    constrained_mols: List[str]
    product: str
    constrained_mols_are_inchi: bool = False

    def __post_init__(self):
        if Chem is not None:
            for i, smiles in enumerate(self.constrained_mols):
                inchi = _smiles_to_inchi(smiles)
                self.constrained_mols[i] = inchi
            self.constrained_mols_are_inchi = True

    def validate_reaction(
        self,
        reactants: List[str],
        reagents: Optional[List[str]],
    ) -> List[str]:
        """
        Checks a given reaction for reagent validity, feasibility of reaction, and other constraints given by the user.
        All participating molecules are given as SMILES strings.
        Returns a list of problems with the reaction if found. A valid reaction returns an empty list.

        :param reactants: A list of reactants involved in this reaction.
        :param reagents: An optional list of other reagents (e.g., solvents, catalysts) used in this reaction.
        :return: Returns an empty list on success, or a list of problems with the reaction in case of failure.
        """
        reagents = reagents or []

        issues: List[str] = []

        # Check for validity of each SMILES string and constraint
        def _check_list(molecules: List[str], mtype: str):
            for m in molecules:
                if not verify_smiles(m):
                    issues.append(f"Invalid SMILES string for {mtype}: {m}")
                if self.constrained_mols_are_inchi:
                    if _smiles_to_inchi(m) in self.constrained_mols:
                        issues.append(f"Disallowed {mtype} {m} due to constraints")
                else:
                    if m in self.constrained_mols:
                        issues.append(f"Disallowed {mtype} {m} due to constraints")

        _check_list(reactants, "reactant")
        _check_list(reagents, "reagent")

        if issues:
            return issues

        # TODO: Check for reaction validity with neural network
        # if invalid_reaction(reactants, reagents, self.product):
        #     issues.append("Reaction is unlikely to yield target product due to atomic structure")

        return issues
