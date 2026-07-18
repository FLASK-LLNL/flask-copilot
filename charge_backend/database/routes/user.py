from fastapi import APIRouter, HTTPException
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from charge_backend.database.deps import GetSession, CurrentUser
from charge_backend.database.models import User, UserResponse, UserCreate

router = APIRouter(prefix="/users", tags=["users"])


async def create_user(session: AsyncSession, user: UserCreate) -> User:
    db_user = User(**user.model_dump())
    session.add(db_user)
    await session.commit()
    await session.refresh(db_user)
    return db_user


async def get_user_by_username(session: AsyncSession, username: str) -> User:
    statement = select(User).where(User.name == username)
    db_user = await session.scalar(statement)

    return db_user


@router.post("/", response_model=UserResponse)
async def create_user_endpoint(*, session: GetSession, user: UserCreate):
    existing_user = await get_user_by_username(session, user.name)
    if existing_user is not None:
        raise HTTPException(status_code=409, detail="username already exists")
    return await create_user(session, user)


# NOTE (trb): There's nothing to update for Users (they're just a
# username tied to a database identifier), so we omit that API
# entirely. This can obviously change if there's a reason going
# forward.


@router.delete("/me/")
async def delete_user_me(session: GetSession, user: CurrentUser):
    await session.delete(user)
    await session.commit()
    return {"ok": True}
