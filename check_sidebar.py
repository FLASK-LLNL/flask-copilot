#!/usr/bin/env python3
"""Quick script to check sidebar_state in DB for experiment 4."""
import pymysql, ssl, json

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

pw = "Eked1c2OWATXtD0YhHKP5CUKh5FGlbIkTaIDGtl1vKMHSB5lrW1FmA8RJB5k0V4x0lgxNSkMAhYbxo4f"

conn = pymysql.connect(
    host='127.0.0.1', port=32636, user='marathe1',
    password=pw, database='flaskcopilot', ssl=ctx,
)
cur = conn.cursor()

# Check all 4 recent experiments
exp_ids = [
    ('exp_1771964788350_1o1zry0d6', 'Experiment 1'),
    ('exp_1771964789266_8k8v8liw7', 'Experiment 2'),
    ('exp_1771964790243_6xymtx8uh', 'Experiment 3'),
    ('exp_1771964791095_l3240j2z3', 'Experiment 4'),
]

for exp_id, label in exp_ids:
    cur.execute("SELECT sidebar_state FROM experiments WHERE id = %s", (exp_id,))
    row = cur.fetchone()
    if row and row[0]:
        sidebar = json.loads(row[0]) if isinstance(row[0], str) else row[0]
        messages = sidebar.get('messages', [])
        retro = [m for m in messages if 'Retrosynthesis Complete' in m.get('message', '')]
        optim = [m for m in messages if 'Optimization Complete' in m.get('message', '')]
        has_complete = bool(retro or optim)
        print(f"\n=== {label} ({exp_id}) ===")
        print(f"  Total messages: {len(messages)}")
        print(f"  Has completion message: {has_complete}")
        if has_complete:
            for m in (retro + optim):
                print(f"  -> {m.get('message', '')[:100]}")
        else:
            print("  MISSING completion message!")
            # Show last 3 messages
            print("  Last 3 messages:")
            for m in messages[-3:]:
                src = m.get('source', '?')
                txt = m.get('message', '')[:100]
                print(f"    [{src}] {txt}")
    else:
        print(f"\n=== {label} ({exp_id}) ===")
        print("  No sidebar_state found!")

conn.close()
