"""
Realistic end-to-end test:
1. Create an experiment in MariaDB (simulating what the frontend does)
2. Run a WS computation that sets serverSessionId matching the experiment
3. Disconnect mid-computation (simulating browser close)
4. Verify the backend saved nodes/edges to the DB
5. Query the REST API to verify the experiment state is loadable
"""
import asyncio
import json
import time
import requests
import websockets


BASE_URL = "http://127.0.0.1:8001"
WS_URL = "ws://127.0.0.1:8001/ws"
SMILES = "CCO"


def create_experiment():
    """Create an experiment via the save_session API endpoint.
    Returns the actual experiment ID assigned by the backend."""
    payload = {
        "state": {
            "smiles": SMILES,
            "problemType": "retrosynthesis",
            "isComputing": True,
            "treeNodes": [],
            "edges": [],
            "metricsHistory": [],
            "visibleMetrics": {},
            "autoZoom": True,
        },
    }

    resp = requests.post(f"{BASE_URL}/api/sessions/save", json=payload)
    data = resp.json()
    session_id = data.get("sessionId")
    print(f"Created experiment {session_id}: HTTP {resp.status_code}")
    if resp.status_code != 200:
        print(f"  Body: {resp.text[:200]}")
    return session_id


async def run_computation_then_disconnect(experiment_id, disconnect_after=4):
    """
    Start a WS computation. Once we get the session_started message,
    update the experiment's serverSessionId in the DB (as the frontend does).
    Then disconnect after receiving disconnect_after node messages.
    """
    async with websockets.connect(WS_URL) as ws:
        await ws.send(json.dumps({
            "action": "compute",
            "smiles": SMILES,
            "problemType": "retrosynthesis",
            "depth": 2,
        }))

        node_count = 0
        server_session_id = None

        while True:
            msg = await asyncio.wait_for(ws.recv(), timeout=30)
            data = json.loads(msg)

            if data.get("type") == "session_started":
                server_session_id = data.get("sessionId")
                print(f"WS session started: {server_session_id}")

                # Update experiment's serverSessionId (as frontend does via save)
                payload = {
                    "sessionId": experiment_id,
                    "state": {
                        "smiles": SMILES,
                        "problemType": "retrosynthesis",
                        "isComputing": True,
                        "treeNodes": [],
                        "edges": [],
                        "metricsHistory": [],
                        "visibleMetrics": {},
                        "serverSessionId": server_session_id,
                        "autoZoom": True,
                    },
                }
                resp = requests.post(f"{BASE_URL}/api/sessions/save", json=payload)
                print(f"  Updated serverSessionId in DB: HTTP {resp.status_code}")

            elif data.get("type") == "node":
                nid = data.get("node", {}).get("id", "?")
                node_count += 1
                print(f"  Received node #{node_count}: {nid}")

                if node_count >= disconnect_after:
                    print(f"  Disconnecting after {node_count} nodes (simulating browser close)...")
                    break

            elif data.get("type") in ("edge", "edge_update", "response"):
                pass
            elif data.get("type") == "complete":
                print("  Computation completed before disconnect target")
                break

    # Give the backend a moment to run save_session_to_db
    await asyncio.sleep(2)
    return server_session_id


def verify_db_state(experiment_id):
    """Check what's in the DB for this experiment."""
    import pymysql
    conn = pymysql.connect(
        host="127.0.0.1", port=3306,
        user="marathe1", password="marathe1",
        database="flaskcopilot",
    )
    cur = conn.cursor()
    cur.execute(
        "SELECT id, is_running, JSON_LENGTH(tree_nodes), JSON_LENGTH(edges), "
        "JSON_EXTRACT(graph_state, '$.serverSessionId') "
        "FROM experiments WHERE id = %s",
        (experiment_id,),
    )
    row = cur.fetchone()
    conn.close()
    if row:
        print(f"\nDB state for {experiment_id}:")
        print(f"  is_running: {row[1]}")
        print(f"  tree_nodes count: {row[2]}")
        print(f"  edges count: {row[3]}")
        print(f"  serverSessionId: {row[4]}")
        return row[2], row[3]
    else:
        print(f"\nNo DB row found for {experiment_id}")
        return 0, 0


def verify_api_response():
    """Check that the REST API returns the saved state."""
    resp = requests.get(f"{BASE_URL}/api/sessions/latest")
    if resp.status_code == 200:
        data = resp.json()
        state = data.get("state", {})
        tn = state.get("treeNodes") or []
        ed = state.get("edges") or []
        print(f"\nAPI /sessions/latest response:")
        print(f"  sessionId: {data.get('sessionId')}")
        print(f"  treeNodes: {len(tn)}")
        print(f"  edges: {len(ed)}")
        print(f"  isComputing: {state.get('isComputing')}")
    else:
        print(f"API error: {resp.status_code}")


async def main():
    print("=" * 60)
    print("END-TO-END TEST: Backend DB Persistence on WS Disconnect")
    print("=" * 60)

    # Step 1: Create experiment
    print("\n--- Step 1: Create experiment in DB ---")
    exp_id = create_experiment()

    # Step 2: Run computation, disconnect after 4 nodes
    print("\n--- Step 2: Run WS computation, disconnect after 4 nodes ---")
    server_sid = await run_computation_then_disconnect(exp_id, disconnect_after=4)

    # Step 3: Verify DB state
    print("\n--- Step 3: Verify DB state ---")
    node_count, edge_count = verify_db_state(exp_id)

    # Step 4: Verify API
    print("\n--- Step 4: Verify API response ---")
    verify_api_response()

    # Summary
    print("\n" + "=" * 60)
    if node_count and node_count > 0:
        print(f"SUCCESS: Backend saved {node_count} nodes and {edge_count} edges to DB on disconnect")
    else:
        print("FAILURE: No nodes were saved to DB on disconnect")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())
