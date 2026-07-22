from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from .routes import user, project, experiment


app = FastAPI(title="Tom Plays with CRUD")

# Enable Cross-Origin Resource Sharing (CORS)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.post("/make-coffee")
def make_tom_a_coffee():
    raise HTTPException(status_code=418, detail="I am a teapot")


app.include_router(user.router)
app.include_router(project.router)
app.include_router(experiment.router)
