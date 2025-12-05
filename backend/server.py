################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################
"""
General server and ``/`` routing infrastructure.
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager

from backend.database.engine import engine, Base
from backend.routers import webui

@asynccontextmanager
async def lifespan(app: FastAPI):
    if engine is not None:
        # Startup: Create tables
        async with engine.begin() as conn:
            await conn.run_sync(Base.metadata.create_all)

    yield

    if engine is not None:
        # Shutdown: Close connections
        await engine.dispose()


app = FastAPI(title="Flask Copilot Backend", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/health")
async def health_check():
    return {"status": "healthy"}

app.include_router(webui.router)

if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8001)
