#!/bin/sh

# Launch main server
uvicorn $@ -- --json_file known_molecules.json --config /aizynth/config.yml --backend openai --model gpt-5-nano &

# Wait for server to start up
sleep 10

# Launch MCP servers
python3 /app/charge-backend/retro_tool_server.py --config /aizynth/config.yml --copilot-host localhost --copilot-port 8001 &
python3 /app/charge-backend/lmo_tool_server.py --copilot-host localhost --copilot-port 8001
