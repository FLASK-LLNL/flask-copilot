#!/usr/bin/env python3
"""
Quick test: verify experiments are saved to DB and can be loaded correctly.
Simulates the browser flow using direct API + WebSocket calls.
"""
import asyncio
import json
import sys
import httpx
import websockets

HTTP = "http://localhost:8001"
WS = "ws://localhost:8001/ws"
MOLECULE = "CC.O=[N+]([O-])C(COc1nc(OCC([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])nc(C([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])n1)([N+](=O)[O-])[N+](=O)[O-].[Na+].[N-]=[N+]=[N-].O"
USER = "marathe1"
HEADERS = {"X-Forwarded-User": USER}


async def create_project(client: httpx.AsyncClient) -> str:
    resp = await client.post(f"{HTTP}/api/projects", json={"name": "Test Project"}, headers=HEADERS)
    resp.raise_for_status()
    proj = resp.json()
    print(f"  Created project: {proj['id']}")
    return proj["id"]


async def create_experiment(client: httpx.AsyncClient, project_id: str, name: str) -> str:
    resp = await client.post(f"{HTTP}/api/projects/{project_id}/experiments", json={"name": name}, headers=HEADERS)
    resp.raise_for_status()
    exp = resp.json()
    print(f"  Created experiment: {exp['id']} ({name})")
    return exp["id"]


async def run_experiment(ws, experiment_id: str, label: str):
    """Send compute command and wait for completion."""
    msg = {
        "action": "compute",
        "smiles": MOLECULE,
        "problemType": "retrosynthesis",
        "experimentId": experiment_id,
    }
    await ws.send(json.dumps(msg))
    print(f"  [{label}] Compute started for experiment {experiment_id}")

    nodes = []
    edges = []
    sidebar_msgs = []
    while True:
        raw = await asyncio.wait_for(ws.recv(), timeout=60)
        data = json.loads(raw)
        t = data.get("type")
        if t == "node":
            nodes.append(data["node"]["id"])
        elif t == "edge":
            edges.append(data["edge"]["id"])
        elif t == "response":
            sidebar_msgs.append(data.get("message", {}).get("message", "")[:60])
        elif t == "complete":
            print(f"  [{label}] Complete: {len(nodes)} nodes, {len(edges)} edges, {len(sidebar_msgs)} sidebar msgs")
            return nodes, edges, sidebar_msgs
        elif t == "session_started":
            print(f"  [{label}] Session started: {data.get('sessionId')}")
        # Ignore other types


async def check_db_state(client: httpx.AsyncClient, project_id: str) -> list:
    resp = await client.get(f"{HTTP}/api/projects", headers=HEADERS)
    resp.raise_for_status()
    projects = resp.json()
    proj = next((p for p in projects if p["id"] == project_id), None)
    if not proj:
        print("  ERROR: Project not found in DB!")
        return []
    
    exps = proj["experiments"]
    for e in exps:
        tn = e.get("treeNodes") or []
        ed = e.get("edges") or []
        sb = (e.get("sidebarState") or {}).get("messages") or []
        print(f"  DB: {e['name']} ({e['id']}): {len(tn)} nodes, {len(ed)} edges, {len(sb)} sidebar msgs, isRunning={e.get('isRunning')}")
    return exps


async def main():
    async with httpx.AsyncClient(follow_redirects=True) as client:
        # Step 1: Create project + 3 experiments
        print("=== Step 1: Create project and experiments ===")
        project_id = await create_project(client)
        exp1_id = await create_experiment(client, project_id, "Experiment 1")
        exp2_id = await create_experiment(client, project_id, "Experiment 2")
        exp3_id = await create_experiment(client, project_id, "Experiment 3")

        # Step 2: Connect WS and run experiments sequentially
        print("\n=== Step 2: Run experiments via WebSocket ===")
        async with websockets.connect(WS, additional_headers=HEADERS) as ws:
            # Drain initial messages (tool list, username, etc.)
            # We need to send list-tools and get-username first
            await ws.send(json.dumps({"action": "list-tools"}))
            await ws.send(json.dumps({"action": "get-username"}))

            # Drain responses
            for _ in range(2):
                raw = await asyncio.wait_for(ws.recv(), timeout=5)
                data = json.loads(raw)
                print(f"  Init: {data.get('type')}")

            # Run experiment 1
            nodes1, edges1, msgs1 = await run_experiment(ws, exp1_id, "Exp1")
            await asyncio.sleep(1)

            # Run experiment 2 (this detaches exp1 on the server)
            nodes2, edges2, msgs2 = await run_experiment(ws, exp2_id, "Exp2")
            await asyncio.sleep(1)

            # Run experiment 3
            nodes3, edges3, msgs3 = await run_experiment(ws, exp3_id, "Exp3")
            await asyncio.sleep(1)

        # WS is now disconnected — triggers save_session_to_db for all sessions
        print("\n=== Step 3: WebSocket disconnected, waiting for DB saves ===")
        await asyncio.sleep(3)

        # Step 4: Check DB state
        print("\n=== Step 4: Check DB state after disconnect ===")
        exps = await check_db_state(client, project_id)

        # Step 5: Load via /api/sessions/latest
        print("\n=== Step 5: Check /api/sessions/latest ===")
        resp = await client.get(f"{HTTP}/api/sessions/latest", headers=HEADERS)
        if resp.status_code == 200:
            data = resp.json()
            if data:
                state = data.get("state", {})
                tn = state.get("treeNodes") or []
                sm = state.get("sidebarMessages") or []
                print(f"  Latest session: {data['name']} ({data['sessionId']})")
                print(f"    {len(tn)} nodes, {len(sm)} sidebar messages")
            else:
                print("  Latest session: null")
        else:
            print(f"  Latest session request failed: {resp.status_code}")

        # Step 6: Verify all experiments have data
        print("\n=== Step 6: Verification ===")
        ok = True
        for e in exps:
            tn = e.get("treeNodes") or []
            if len(tn) == 0:
                print(f"  FAIL: {e['name']} has 0 nodes!")
                ok = False
            else:
                print(f"  OK: {e['name']} has {len(tn)} nodes")
        
        if ok:
            print("\n  ALL EXPERIMENTS HAVE DATA")
        else:
            print("\n  SOME EXPERIMENTS ARE MISSING DATA")
            sys.exit(1)


if __name__ == "__main__":
    asyncio.run(main())
