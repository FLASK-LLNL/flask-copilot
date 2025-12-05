#!/bin/sh

. /venv/bin/activate

# Launch main server
export PYTHONPATH=$PYTHONPATH:/app/charge_backend:/app
uvicorn --host 0.0.0.0 --port 8001 --workers 8 charge_server:app &

# Wait for server to start up
sleep 10

# Launch MCP servers
python /app/charge_backend/retro_tool_server.py --config /data/config.yml --copilot-host localhost --copilot-port 8001 &
python /app/charge_backend/lmo_tool_server.py --copilot-host localhost --copilot-port 8001
