from fastapi import APIRouter
from typing import List

from charge_backend.database.deps import (
    GetSession,
    CurrentUser,
    ValidatedProject,
    ValidatedExperiment,
)
from charge_backend.database.models import (
    Experiment,
    ExperimentCreate,
    ExperimentResponse,
    ExperimentUpdate,
)

# The copilot interaction always carries the project_id in the
# {java,type}script interfaces. So I may as well use the nested route.
router = APIRouter(prefix="/projects/{project_id}/experiments", tags=["experiments"])


@router.post("", response_model=ExperimentResponse)
async def create_experiment(
    *,
    session: GetSession,
    experiment_in: ExperimentCreate,
    validated_project: ValidatedProject,
):

    db_experiment = Experiment(
        project_id=validated_project.id, data=experiment_in.model_dump(by_alias=False)
    )
    session.add(db_experiment)
    await session.commit()
    await session.refresh(db_experiment)
    return db_experiment


# FIXME (trb): I don't know that this is necessary (but it was trivial
# to implement)
@router.get("", response_model=List[ExperimentResponse])
async def get_experiments(
    *,
    session: GetSession,
    validated_project: ValidatedProject,
):
    """
    Retrieve all Experiments associated with 'project_id'
    """
    await session.refresh(validated_project)
    return validated_project.experiments


@router.get(
    "/{experiment_id}",
    response_model=ExperimentResponse,
)
async def get_experiment(*, db_experiment: ValidatedExperiment):
    return db_experiment


@router.put(
    "/{experiment_id}",
    response_model=ExperimentResponse,
)
async def update_experiment(
    *,
    session: GetSession,
    db_experiment: ValidatedExperiment,
    exp_update: ExperimentUpdate,
):
    update_dict = exp_update.model_dump(by_alias=False, exclude_unset=True)

    if "data" in update_dict:
        new_data = update_dict.pop("data")

        # Merge dictionaries with update semantics
        db_experiment.data = db_experiment.data | new_data
        # NOTE (trb): CANNOT update in-place because sqlalchemy will
        # not recognize that subfields in the JSON have changed! To
        # remedy this, we could use a MutableDict type or just
        # reassign the field.

    # Other fields use assignment semantics
    for key, value in update_dict.items():
        setattr(db_experiment, key, value)

    session.add(db_experiment)
    await session.commit()
    await session.refresh(db_experiment)

    return db_experiment


@router.delete("/{experiment_id}")
async def delete_experiment(
    *,
    session: GetSession,
    db_experiment: ValidatedExperiment,
):
    await session.delete(db_experiment)
    await session.commit()
    return {"ok": True}
