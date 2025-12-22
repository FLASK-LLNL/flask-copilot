################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################
"""
Web UI (``/``) routing
"""
from fastapi import Request, APIRouter
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
import os
from loguru import logger

router = APIRouter(tags=["webui"])

if "FLASK_APPDIR" in os.environ:
    DIST_PATH = os.environ["FLASK_APPDIR"]
else:
    DIST_PATH = os.path.join(os.path.dirname(__file__), "flask-app", "dist")
ASSETS_PATH = os.path.join(DIST_PATH, "assets")

if os.path.exists(ASSETS_PATH):
    # Serve the frontend
    router.mount("/assets", StaticFiles(directory=ASSETS_PATH), name="assets")
    router.mount("/rdkit", StaticFiles(directory=os.path.join(DIST_PATH, "rdkit")), name="rdkit")

    @router.get("/")
    async def root(request: Request):
        logger.info(f"Request for Web UI received. Headers: {str(request.headers)}")
        with open(os.path.join(DIST_PATH, "index.html"), "r") as fp:
            html = fp.read()

        html = html.replace(
            "<!-- APP CONFIG -->",
            f"""
           <script>
           window.APP_CONFIG = {{
               WS_SERVER: '{os.getenv("WS_SERVER", "ws://localhost:8001/ws")}',
               VERSION: '{os.getenv("SERVER_VERSION", "")}'
           }};
           </script>""",
        )
        return HTMLResponse(html)
