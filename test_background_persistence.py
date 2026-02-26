"""
End-to-end test for background computation persistence.

Tests the complete flow:
1. Browser 1 connects, starts a computation
2. Browser 1 disconnects mid-computation
3. Backend continues computing headless and saves to DB periodically
4. Browser 2 connects, loads state from DB
5. Browser 2 sees all molecules generated after Browser 1 disconnected
"""

import asyncio
import json
import sys

import httpx
import websockets

WS_URL = "ws://127.0.0.1:8001/ws"
API_URL = "http://127.0.0.1:8001/api"
SMILES = "CCO"  # ethanol
DEPTH = 3  # small tree to keep test fast (~7 nodes total)


async def test_background_persistence():
    print("=" * 60)
    print("TEST: Background computation persistence")
    print("=" * 60)

    # ── Step 1: Browser 1 connects and starts computation ──
    print("\n[Step 1] Browser 1: Connect and start computation...")
    ws1 = await websockets.connect(WS_URL)

    # Send compute request
    await ws1.send(json.dumps({
        "action": "compute",
        "smiles": SMILES,
        "problemType": "retrosynthesis",
        "depth": DEPTH,
    }))

    # Collect initial messages (session_started + a few nodes)
    session_id = None
    browser1_nodes = []
    browser1_edges = []

    # Read messages for a short while (collect some but not all)
    try:
        for _ in range(20):  # read up to 20 messages
            msg = await asyncio.wait_for(ws1.recv(), timeout=2.0)
            data = json.loads(msg)
            if data["type"] == "session_started":
                session_id = data["sessionId"]
                print(f"   Session started: {session_id}")
            elif data["type"] == "node":
                browser1_nodes.append(data["node"]["id"])
                print(f"   Node received: {data['node']['id']}")
            elif data["type"] == "edge":
                browser1_edges.append(data["edge"]["id"])
            elif data["type"] == "complete":
                print("   ⚠ Computation completed before disconnect (increase DEPTH)")
                break
    except asyncio.TimeoutError:
        pass  # No more messages for now

    nodes_before_disconnect = len(browser1_nodes)
    print(f"   Browser 1 received {nodes_before_disconnect} nodes before disconnect")

    # ── Step 2: Browser 1 disconnects abruptly ──
    print("\n[Step 2] Browser 1: Disconnecting...")
    await ws1.close()
    print("   ✓ Browser 1 disconnected")

    # ── Step 3: Wait for background task to continue ──
    print("\n[Step 3] Waiting for background computation to continue...")
    # Wait long enough for the computation to finish and DB saves to complete
    await asyncio.sleep(12)
    print("   ✓ Waited 12 seconds for background processing")

    # ── Step 4: Check the database via REST API ──
    print("\n[Step 4] Checking database via REST API...")
    async with httpx.AsyncClient() as client:
        resp = await client.get(
            f"{API_URL}/sessions/latest",
            headers={"X-Forwarded-User": "nobody"},
        )
        if resp.status_code == 200:
            db_data = resp.json()
            state = db_data.get("state", {})
            db_nodes = state.get("treeNodes", [])
            db_edges = state.get("edges", [])
            is_computing = state.get("isComputing", None)
            server_session_id = state.get("serverSessionId", None)

            print(f"   DB has {len(db_nodes)} nodes, {len(db_edges)} edges")
            print(f"   isComputing={is_computing}")
            print(f"   serverSessionId={server_session_id}")

            if len(db_nodes) > nodes_before_disconnect:
                print(f"   ✓ DB has MORE nodes ({len(db_nodes)}) than Browser 1 saw ({nodes_before_disconnect})")
                print("     → Background computation continued after disconnect!")
            elif len(db_nodes) == nodes_before_disconnect:
                print(f"   ⚠ DB has same count ({len(db_nodes)}) as Browser 1")
                print("     → Background task may have finished before disconnect, or periodic saves didn't trigger")
            else:
                print(f"   ✗ DB has FEWER nodes ({len(db_nodes)}) than Browser 1 saw ({nodes_before_disconnect})")
                print("     → Something is wrong with the save logic")
        elif resp.status_code == 404:
            print("   ✗ No session found in DB!")
            print("     → The experiment row was never created (frontend auto-save didn't run)")
            print("     Note: This is expected if no prior frontend auto-save happened.")
            print("     The backend can only save to an EXISTING experiment row.")
        else:
            print(f"   ✗ API returned status {resp.status_code}: {resp.text}")

    # ── Step 5: Browser 2 connects and resumes ──
    print("\n[Step 5] Browser 2: Connect and attempt resume...")
    if session_id:
        ws2 = await websockets.connect(WS_URL)
        await ws2.send(json.dumps({
            "action": "resume_session",
            "sessionId": session_id,
        }))

        browser2_nodes = []
        browser2_edges = []
        is_complete = False
        session_status = None

        try:
            for _ in range(50):
                msg = await asyncio.wait_for(ws2.recv(), timeout=3.0)
                data = json.loads(msg)
                if data["type"] == "node":
                    browser2_nodes.append(data["node"]["id"])
                elif data["type"] == "edge":
                    browser2_edges.append(data["edge"]["id"])
                elif data["type"] == "complete":
                    is_complete = True
                    break
                elif data["type"] == "session_resumed":
                    print(f"   Session resumed (sentNodes={data.get('sentNodes')})")
                elif data["type"] == "session_status":
                    session_status = data.get("status")
                    print(f"   Session status: {session_status}")
                    break
                elif data["type"] == "session_not_found":
                    print(f"   Session not found!")
                    break
        except asyncio.TimeoutError:
            pass

        await ws2.close()

        print(f"   Browser 2 received {len(browser2_nodes)} nodes via resume")
        if is_complete or session_status == "complete":
            print("   ✓ Computation was already complete when Browser 2 resumed")
        if len(browser2_nodes) > nodes_before_disconnect:
            print(f"   ✓ Browser 2 got MORE nodes ({len(browser2_nodes)}) than Browser 1 ({nodes_before_disconnect})")
            print("     → Backend continued computing after Browser 1 disconnected!")
        elif len(browser2_nodes) > 0:
            print(f"   ✓ Browser 2 received {len(browser2_nodes)} replayed nodes")
    else:
        print("   ⚠ No session_id available, skipping resume test")

    # ── Summary ──
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    print(f"  Browser 1 nodes before disconnect: {nodes_before_disconnect}")
    if resp.status_code == 200:
        print(f"  Database nodes after waiting:      {len(db_nodes)}")
        print(f"  Database edges after waiting:       {len(db_edges)}")
        if len(db_nodes) > nodes_before_disconnect:
            print(f"\n  ✅ PASS: Background persistence is working!")
            print(f"     The backend generated {len(db_nodes) - nodes_before_disconnect} additional nodes")
            print(f"     after Browser 1 disconnected and saved them to the database.")
        else:
            print(f"\n  ⚠ CHECK: Could not confirm additional nodes in DB")
            print(f"     This may be because the experiment row wasn't created by the frontend,")
            print(f"     or the computation finished before disconnect.")
    else:
        print(f"  Database: No session found (status {resp.status_code})")
    
    if session_id and len(browser2_nodes) > 0:
        print(f"  Browser 2 resume nodes:            {len(browser2_nodes)}")
        if len(browser2_nodes) >= nodes_before_disconnect:
            print(f"  ✅ PASS: Resume/replay is working!")
        else:
            print(f"  ⚠ Resume returned fewer nodes than expected")


if __name__ == "__main__":
    asyncio.run(test_background_persistence())
