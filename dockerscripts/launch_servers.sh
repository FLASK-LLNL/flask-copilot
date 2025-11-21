#!/bin/sh

# Launch main server
export PYTHONPATH=$PYTHONPATH:/app/charge-backend
uvicorn --host 0.0.0.0 --port 8001 --workers 8 charge_server:app &

# Wait for server to start up
sleep 10

# Launch MCP servers
python3 /app/charge-backend/retro_tool_server.py --config /aizynth/config.yml --copilot-host localhost --copilot-port 8001 &
python3 /app/charge-backend/lmo_tool_server.py --copilot-host localhost --copilot-port 8001
