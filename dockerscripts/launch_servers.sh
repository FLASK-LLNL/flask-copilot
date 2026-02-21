#!/bin/sh

. /venv/bin/activate

# Launch main server
export PYTHONPATH=$PYTHONPATH:/app/charge_backend:/app
uvicorn --host 0.0.0.0 --port 8001 --workers 8 charge_server:app
