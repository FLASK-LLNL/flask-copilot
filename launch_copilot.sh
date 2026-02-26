#!/bin/bash
set -e

REPO="/Users/marathe1/work/flask/flask-copilot"

# --- DB target: "local" (default) or "remote" --------------------------
DB_TARGET="${1:-local}"
case "$DB_TARGET" in
  local)
    echo "Using LOCAL MariaDB config (.env-local)"
    cp "$REPO/.env-local" "$REPO/.env"
    ;;
  remote)
    echo "Using REMOTE LaunchIT MariaDB config (.env-launchit)"
    cp "$REPO/.env-launchit" "$REPO/.env"
    ;;
  *)
    echo "Usage: $0 [local|remote]  (default: local)"
    exit 1
    ;;
esac

# Kill any stale processes on the ports we need before starting.
# This ensures Vite always binds to 5173 (not a fallback port) and
# that the backend gets port 8001 cleanly.
echo "Clearing stale processes on ports 8001 and 5173..."
lsof -ti:8001 | xargs kill -9 2>/dev/null || true
lsof -ti:5173 | xargs kill -9 2>/dev/null || true
pkill -9 -f "vite" 2>/dev/null || true
pkill -9 -f "npm run dev" 2>/dev/null || true
sleep 1

# Confirm ports are free
if lsof -i:8001 | grep -q LISTEN 2>/dev/null; then
  echo "ERROR: Port 8001 still occupied after kill attempt. Aborting."
  exit 1
fi
if lsof -i:5173 | grep -q LISTEN 2>/dev/null; then
  echo "ERROR: Port 5173 still occupied after kill attempt. Aborting."
  exit 1
fi

# Launch frontend in background with output to fixed log file
echo "Starting frontend..."
cd "$REPO/flask-app" && npm run dev > "$REPO/npm_dev.log" 2>&1 &
FRONTEND_PID=$!
echo "Frontend started (PID: $FRONTEND_PID). Logs: $REPO/npm_dev.log"

# Wait for Vite to be ready and confirm it bound to 5173 (not a fallback port)
echo "Waiting for frontend on port 5173..."
for i in $(seq 1 15); do
  if lsof -i:5173 | grep -q LISTEN 2>/dev/null; then
    echo "Frontend is up on port 5173."
    break
  fi
  if [ "$i" -eq 15 ]; then
    echo "WARNING: Frontend did not bind to 5173 after 15s. Check $REPO/npm_dev.log"
  fi
  sleep 1
done

# Launch backend in foreground (output directly to terminal)
echo "Starting backend server (ctrl-c or closing terminal will stop both)..."
cd "$REPO" && source .venv/bin/activate && python3 mock_server.py

# If backend exits, clean up frontend
echo "Backend stopped. Cleaning up frontend..."
kill $FRONTEND_PID 2>/dev/null || true
wait $FRONTEND_PID 2>/dev/null || true
