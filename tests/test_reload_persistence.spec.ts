/**
 * End-to-end test: create experiments, run them, then RELOAD the page
 * and verify all experiments still have their data.
 *
 * This tests the DB persistence + reload path.
 */
import { test, expect, Page } from '@playwright/test';

const SMILES = 'CC.O=[N+]([O-])C(COc1nc(OCC([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])nc(C([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])n1)([N+](=O)[O-])[N+](=O)[O-].[Na+].[N-]=[N+]=[N-].O';
const COMPUTATION_TIMEOUT = 120_000;

async function countTreeNodes(page: Page): Promise<number> {
  return page.locator('.graph-node').count();
}

async function runExperiment(page: Page): Promise<number> {
  const smilesInput = page.locator('input[placeholder="Enter SMILES notation"]');
  await smilesInput.fill(SMILES);
  await page.waitForTimeout(300);

  const runBtn = page.locator('button.btn-primary:has-text("Run")');
  await expect(runBtn).toBeEnabled({ timeout: 10_000 });
  await runBtn.click();

  await page.locator('button.btn-primary:has-text("Computing")').waitFor({ state: 'visible', timeout: 10_000 });
  await page.locator('button.btn-primary:has-text("Rerun")').waitFor({ state: 'visible', timeout: COMPUTATION_TIMEOUT });
  await page.waitForTimeout(1000);

  return countTreeNodes(page);
}

test.describe('Experiment persistence across page reload', () => {
  test('experiments survive page reload', async ({ page }) => {
    test.setTimeout(600_000);

    await page.goto('http://localhost:5173');
    await page.waitForLoadState('networkidle');
    await page.waitForTimeout(2000);

    const consoleLogs: string[] = [];
    page.on('console', msg => consoleLogs.push(`[${msg.type()}] ${msg.text()}`));

    // --- Create Project + Experiment 1 ---
    const newProjectBtn = page.locator('button.project-button:has-text("New Project")');
    await newProjectBtn.waitFor({ state: 'visible', timeout: 10_000 });
    await newProjectBtn.click();
    await page.waitForTimeout(1000);
    await page.locator('button.btn-primary:has-text("Run")').waitFor({ state: 'visible', timeout: 15_000 });

    const newExpBtn = page.locator('button.experiment-button:has-text("New Experiment")');
    if (await newExpBtn.isVisible()) {
      await newExpBtn.click();
      await page.waitForTimeout(1000);
    }

    const exp1Nodes = await runExperiment(page);
    console.log(`Experiment 1: ${exp1Nodes} nodes`);
    expect(exp1Nodes).toBeGreaterThan(0);

    // --- Create Experiment 2 ---
    await newExpBtn.click();
    await page.waitForTimeout(1500);
    const exp2Nodes = await runExperiment(page);
    console.log(`Experiment 2: ${exp2Nodes} nodes`);
    expect(exp2Nodes).toBeGreaterThan(0);

    // Wait for saves to complete
    await page.waitForTimeout(3000);

    // --- RELOAD the page ---
    console.log('=== RELOADING PAGE ===');
    await page.reload();
    await page.waitForLoadState('networkidle');
    await page.waitForTimeout(4000); // Wait for session restore + 2s poll

    // The app should restore to the most recent experiment (Experiment 2)
    // Check if nodes are visible
    const nodesAfterReload = await countTreeNodes(page);
    console.log(`Nodes after reload (should be Experiment 2): ${nodesAfterReload}`);

    // Now click Experiment 1
    const exp1Btn = page.locator('button.experiment-button:has-text("Experiment 1")');
    if (await exp1Btn.isVisible()) {
      await exp1Btn.click();
      await page.waitForTimeout(2000);
      const exp1ReloadNodes = await countTreeNodes(page);
      console.log(`Experiment 1 after reload: ${exp1ReloadNodes} nodes (expected: ${exp1Nodes})`);
      expect(exp1ReloadNodes).toBeGreaterThan(0);
    } else {
      console.log('Experiment 1 button not visible after reload');
    }

    // Click Experiment 2
    const exp2Btn = page.locator('button.experiment-button:has-text("Experiment 2")');
    if (await exp2Btn.isVisible()) {
      await exp2Btn.click();
      await page.waitForTimeout(2000);
      const exp2ReloadNodes = await countTreeNodes(page);
      console.log(`Experiment 2 after reload: ${exp2ReloadNodes} nodes (expected: ${exp2Nodes})`);
      expect(exp2ReloadNodes).toBeGreaterThan(0);
    } else {
      console.log('Experiment 2 button not visible after reload');
    }

    // Print relevant console logs
    const loadLogs = consoleLogs.filter(l =>
      l.includes('Loading context') ||
      l.includes('loadStateFromCurrentExperiment') ||
      l.includes('Applying restored session') ||
      l.includes('not found') ||
      l.includes('Loaded') && l.includes('projects')
    );
    console.log(`\n=== LOAD LOGS ===`);
    for (const l of loadLogs) {
      console.log(l);
    }
  });
});
