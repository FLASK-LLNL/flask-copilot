from sqlalchemy.ext.asyncio import create_async_engine, async_sessionmaker
from .models import Base

sql_db_filename = "database.db"
sqlite_url = f"sqlite+aiosqlite:///{sql_db_filename}"

connect_args = {"check_same_thread": False}
engine = create_async_engine(
    sqlite_url, echo=True, connect_args=connect_args, pool_pre_ping=True
)
AsyncSessionLocal = async_sessionmaker(autocommit=False, autoflush=False, bind=engine)


async def get_session():
    async with AsyncSessionLocal() as session:
        try:
            yield session
        except:
            session.close()
