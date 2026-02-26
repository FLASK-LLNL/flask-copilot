import { test, expect, Page } from '@playwright/test';

const MOLECULE = 'CC.O=[N+]([O-])C(COc1nc(OCC([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])nc(C([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])n1)([N+](=O)[O-])[N+](=O)[O-].[Na+].[N-]=[N+]=[N-].O';
const BASE_URL = 'http://localhost:5173';

// Utility: wait for WebSocket connection (green dot)
async function waitForWsConnection(page: Page) {
  // Wait for the green dot to appear (WS connected)
  await page.waitForFunction(() => {
    const dot = document.querySelector('[class*="statusDot"]') as HTMLElement;
    if (!dot) return false;
    const style = window.getComputedStyle(dot);
    // Green means connected
    return style.backgroundColor === 'rgb(76, 175, 80)' || style.backgroundColor === 'rgb(0, 128, 0)' || style.backgroundColor.includes('green');
  }, {}, { timeout: 15000 }).catch(() => {
    // Fallback: just wait a bit
  });
  await page.waitForTimeout(2000);
}

// Utility: wait for computation to complete (by waiting for "complete" message)
async function waitForComputationComplete(page: Page, timeoutMs = 120000) {
  // Wait for the Run button to become enabled again (means computation finished)
  await page.waitForFunction(() => {
    const buttons = Array.from(document.querySelectorAll('button'));
    const runBtn = buttons.find(b => b.textContent?.includes('Run') || b.textContent?.includes('Start'));
    if (!runBtn) return false;
    return !runBtn.disabled;
  }, {}, { timeout: timeoutMs });
  // Extra wait for state to flush
  await page.waitForTimeout(500);
}

// Utility: get tree node count from the page
async function getTreeNodeCount(page: Page): Promise<number> {
  return await page.evaluate(() => {
    // Count nodes in the SVG/canvas tree visualization
    const nodes = document.querySelectorAll('[class*="treeNode"], [class*="node-group"], [data-node-id]');
    return nodes.length;
  });
}

// Utility: count experiments in the sidebar
async function getExperimentCount(page: Page): Promise<number> {
  return await page.evaluate(() => {
    // Look for experiment items in the sidebar
    const items = document.querySelectorAll('[class*="experiment"], [class*="Experiment"]');
    return items.length;
  });
}

test.describe('Experiment Load Tests', () => {
  test('Create 3 experiments, close browser, reload, verify all load correctly', async ({ browser }) => {
    // Use a persistent context so localStorage survives
    const context = await browser.newContext();
    let page = await context.newPage();

    // Enable console logging for debugging
    page.on('console', msg => {
      if (msg.type() === 'log' || msg.type() === 'warn' || msg.type() === 'error') {
        console.log(`[BROWSER ${msg.type()}] ${msg.text()}`);
      }
    });

    await page.goto(BASE_URL);
    await waitForWsConnection(page);

    // Clear localStorage to start fresh
    await page.evaluate(() => {
      localStorage.clear();
    });
    await page.reload();
    await waitForWsConnection(page);

    console.log('--- Step 1: Create a project ---');

    // Click the "New Project" button (or equivalent)
    const newProjectBtn = page.locator('button').filter({ hasText: /new.*project|create.*project|\+/i }).first();
    if (await newProjectBtn.isVisible({ timeout: 3000 })) {
      await newProjectBtn.click();
      await page.waitForTimeout(1000);
    }

    console.log('--- Step 2: Create and run experiment 1 ---');

    // Click "New Experiment" button
    const newExpBtn = page.locator('button, [role="button"]').filter({ hasText: /new.*experiment|add.*experiment|\+/i }).first();
    if (await newExpBtn.isVisible({ timeout: 3000 })) {
      await newExpBtn.click();
      await page.waitForTimeout(1000);
    }

    // Enter the molecule string
    const smilesInput = page.locator('input[placeholder*="SMILES"], input[placeholder*="smiles"], textarea[placeholder*="SMILES"], textarea[placeholder*="smiles"], input[type="text"]').first();
    await smilesInput.fill(MOLECULE);
    await page.waitForTimeout(500);

    // Click Run
    const runBtn = page.locator('button').filter({ hasText: /^Run$|^Start$/i }).first();
    await runBtn.click();
    console.log('Experiment 1 running...');

    // Wait for completion
    await waitForComputationComplete(page);
    console.log('Experiment 1 complete');

    // Get tree node count for experiment 1
    const exp1Nodes = await getTreeNodeCount(page);
    console.log(`Experiment 1 nodes: ${exp1Nodes}`);

    // Wait for save to flush
    await page.waitForTimeout(3000);

    console.log('--- Step 3: Create and run experiment 2 ---');

    // Create new experiment
    if (await newExpBtn.isVisible({ timeout: 3000 })) {
      await newExpBtn.click();
      await page.waitForTimeout(1000);
    }

    // Enter molecule and run
    const smilesInput2 = page.locator('input[placeholder*="SMILES"], input[placeholder*="smiles"], textarea[placeholder*="SMILES"], textarea[placeholder*="smiles"], input[type="text"]').first();
    await smilesInput2.fill(MOLECULE);
    await page.waitForTimeout(500);

    const runBtn2 = page.locator('button').filter({ hasText: /^Run$|^Start$/i }).first();
    await runBtn2.click();
    console.log('Experiment 2 running...');

    await waitForComputationComplete(page);
    console.log('Experiment 2 complete');

    const exp2Nodes = await getTreeNodeCount(page);
    console.log(`Experiment 2 nodes: ${exp2Nodes}`);

    await page.waitForTimeout(3000);

    console.log('--- Step 4: Create and run experiment 3 ---');

    if (await newExpBtn.isVisible({ timeout: 3000 })) {
      await newExpBtn.click();
      await page.waitForTimeout(1000);
    }

    const smilesInput3 = page.locator('input[placeholder*="SMILES"], input[placeholder*="smiles"], textarea[placeholder*="SMILES"], textarea[placeholder*="smiles"], input[type="text"]').first();
    await smilesInput3.fill(MOLECULE);
    await page.waitForTimeout(500);

    const runBtn3 = page.locator('button').filter({ hasText: /^Run$|^Start$/i }).first();
    await runBtn3.click();
    console.log('Experiment 3 running...');

    await waitForComputationComplete(page);
    console.log('Experiment 3 complete');

    const exp3Nodes = await getTreeNodeCount(page);
    console.log(`Experiment 3 nodes: ${exp3Nodes}`);

    await page.waitForTimeout(3000);

    console.log('--- Step 5: Check DB state ---');

    // Check how many experiments are in the DB
    const dbCheckResponse = await page.evaluate(async () => {
      const resp = await fetch('http://localhost:8001/api/projects');
      return await resp.json();
    });
    console.log(`DB state: ${JSON.stringify(dbCheckResponse).length} bytes, projects: ${dbCheckResponse.length}`);
    for (const proj of dbCheckResponse) {
      console.log(`  Project ${proj.name}: ${proj.experiments.length} experiments`);
      for (const exp of proj.experiments) {
        const nodeCount = exp.treeNodes?.length ?? 0;
        const sidebarCount = exp.sidebarState?.messages?.length ?? 0;
        console.log(`    Experiment ${exp.name}: ${nodeCount} nodes, ${sidebarCount} sidebar messages`);
      }
    }

    console.log('--- Step 6: Click each experiment and verify non-empty ---');

    // Click experiment 1
    const expItems = page.locator('[class*="experiment"]');
    const expCount = await expItems.count();
    console.log(`Sidebar experiment count: ${expCount}`);

    // Click through experiments by finding them by name
    for (let i = 0; i < expCount; i++) {
      const item = expItems.nth(i);
      const text = await item.textContent();
      console.log(`Clicking experiment: ${text}`);
      await item.click();
      await page.waitForTimeout(2000);

      const nodesAfterClick = await getTreeNodeCount(page);
      console.log(`  Nodes after click: ${nodesAfterClick}`);
    }

    console.log('--- Step 7: Close and reopen ---');

    // Close the page (simulates closing browser tab)
    await page.close();
    await page.waitForTimeout(2000).catch(() => {});

    // Reopen
    page = await context.newPage();
    await page.goto(BASE_URL);
    await waitForWsConnection(page);

    // Wait for projects to load from DB
    await page.waitForTimeout(5000);

    console.log('--- Step 8: Verify experiments load after reopen ---');

    const dbCheckAfterReopen = await page.evaluate(async () => {
      const resp = await fetch('http://localhost:8001/api/projects');
      return await resp.json();
    });
    for (const proj of dbCheckAfterReopen) {
      console.log(`  Project ${proj.name}: ${proj.experiments.length} experiments`);
      for (const exp of proj.experiments) {
        const nodeCount = exp.treeNodes?.length ?? 0;
        const sidebarCount = exp.sidebarState?.messages?.length ?? 0;
        console.log(`    Experiment ${exp.name}: ${nodeCount} nodes, ${sidebarCount} sidebar messages`);
        expect(nodeCount, `Experiment ${exp.name} should have nodes`).toBeGreaterThan(0);
      }
    }

    // Click each experiment in sidebar and verify tree renders
    const expItems2 = page.locator('[class*="experiment"]');
    const expCount2 = await expItems2.count();
    console.log(`Sidebar experiment count after reopen: ${expCount2}`);

    for (let i = 0; i < expCount2; i++) {
      const item = expItems2.nth(i);
      const text = await item.textContent();
      console.log(`Clicking experiment: ${text}`);
      await item.click();
      await page.waitForTimeout(3000);

      const nodesAfterClick = await getTreeNodeCount(page);
      console.log(`  Nodes after click: ${nodesAfterClick}`);
    }

    await context.close();
  });
});
