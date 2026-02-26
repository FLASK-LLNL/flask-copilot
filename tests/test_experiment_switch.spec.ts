/**
 * End-to-end test: create 3 experiments, run each, then switch between them
 * and verify that each experiment's tree nodes are correctly loaded.
 *
 * This test reproduces the "empty experiment on sidebar click" bug.
 */
import { test, expect, Page } from '@playwright/test';

const SMILES = 'CC.O=[N+]([O-])C(COc1nc(OCC([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])nc(C([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])n1)([N+](=O)[O-])[N+](=O)[O-].[Na+].[N-]=[N+]=[N-].O';

// How long to wait for a computation to finish (mock server is fast).
const COMPUTATION_TIMEOUT = 120_000;

async function waitForWsConnected(page: Page): Promise<void> {
  // The Run button is enabled only when WS is connected and smiles is set.
  // We wait for it to NOT be disabled.
  await page.locator('button.btn-primary:has-text("Run")').waitFor({ state: 'visible', timeout: 15_000 });
}

async function countTreeNodes(page: Page): Promise<number> {
  return page.locator('.graph-node').count();
}

async function countMessages(page: Page): Promise<number> {
  return page.locator('.message-card').count();
}

/**
 * Run a single experiment: fill SMILES, click Run, wait for completion.
 * Returns the final node count.
 */
async function runExperiment(page: Page): Promise<{ nodeCount: number; messageCount: number }> {
  // Fill SMILES
  const smilesInput = page.locator('input[placeholder="Enter SMILES notation"]');
  await smilesInput.fill(SMILES);
  await page.waitForTimeout(300);

  // Click Run
  const runBtn = page.locator('button.btn-primary:has-text("Run")');
  await expect(runBtn).toBeEnabled({ timeout: 10_000 });
  await runBtn.click();

  // Wait for computation to complete (button changes from "Computing" back to "Rerun")
  await page.locator('button.btn-primary:has-text("Computing")').waitFor({ state: 'visible', timeout: 10_000 });
  await page.locator('button.btn-primary:has-text("Rerun")').waitFor({ state: 'visible', timeout: COMPUTATION_TIMEOUT });

  // Give a moment for final saves
  await page.waitForTimeout(1000);

  const nodeCount = await countTreeNodes(page);
  const messageCount = await countMessages(page);
  return { nodeCount, messageCount };
}

test.describe('Experiment switching', () => {
  test('create 3 experiments, run each, then verify all load correctly on switch', async ({ page }) => {
    // Increase default timeout for this long test
    test.setTimeout(600_000);

    // Navigate to app
    await page.goto('http://localhost:5173');
    await page.waitForLoadState('networkidle');
    await page.waitForTimeout(2000);

    // Check console for debugging
    const consoleLogs: string[] = [];
    page.on('console', msg => {
      const text = msg.text();
      consoleLogs.push(`[${msg.type()}] ${text}`);
    });

    // --- Create Project ---
    console.log('=== Creating project ===');
    const newProjectBtn = page.locator('button.project-button:has-text("New Project")');
    await newProjectBtn.waitFor({ state: 'visible', timeout: 10_000 });
    await newProjectBtn.click();
    await page.waitForTimeout(1000);

    // Verify the project was created (should be a project button with active class)
    const projectButtons = page.locator('button.project-button-active');
    await expect(projectButtons).toHaveCount(1, { timeout: 5000 });

    // Wait for WS connection
    await waitForWsConnected(page);

    // --- Create and run Experiment 1 ---
    console.log('=== Creating and running Experiment 1 ===');
    const newExpBtn = page.locator('button.experiment-button:has-text("New Experiment")');
    if (await newExpBtn.isVisible()) {
      await newExpBtn.click();
      await page.waitForTimeout(1000);
    }
    // There should be an active experiment
    const activeExp = page.locator('button.experiment-button-active');
    await expect(activeExp).toHaveCount(1, { timeout: 5000 });
    const exp1Name = await activeExp.locator('span.truncate').innerText();
    console.log(`Experiment 1: "${exp1Name}"`);

    const exp1Result = await runExperiment(page);
    console.log(`Experiment 1 completed: ${exp1Result.nodeCount} nodes, ${exp1Result.messageCount} messages`);
    expect(exp1Result.nodeCount).toBeGreaterThan(0);

    // --- Create and run Experiment 2 ---
    console.log('=== Creating Experiment 2 ===');
    await newExpBtn.click();
    await page.waitForTimeout(1500);

    // Verify the new experiment is selected
    const activeExp2 = page.locator('button.experiment-button-active');
    const exp2Name = await activeExp2.locator('span.truncate').innerText();
    console.log(`Experiment 2: "${exp2Name}"`);

    // The view should be reset (0 nodes)
    const nodesAfterCreate = await countTreeNodes(page);
    console.log(`Nodes after creating Experiment 2: ${nodesAfterCreate}`);
    expect(nodesAfterCreate).toBe(0);

    const exp2Result = await runExperiment(page);
    console.log(`Experiment 2 completed: ${exp2Result.nodeCount} nodes, ${exp2Result.messageCount} messages`);
    expect(exp2Result.nodeCount).toBeGreaterThan(0);

    // --- Create and run Experiment 3 ---
    console.log('=== Creating Experiment 3 ===');
    await newExpBtn.click();
    await page.waitForTimeout(1500);

    const activeExp3 = page.locator('button.experiment-button-active');
    const exp3Name = await activeExp3.locator('span.truncate').innerText();
    console.log(`Experiment 3: "${exp3Name}"`);

    const nodesAfterCreate3 = await countTreeNodes(page);
    console.log(`Nodes after creating Experiment 3: ${nodesAfterCreate3}`);
    expect(nodesAfterCreate3).toBe(0);

    const exp3Result = await runExperiment(page);
    console.log(`Experiment 3 completed: ${exp3Result.nodeCount} nodes, ${exp3Result.messageCount} messages`);
    expect(exp3Result.nodeCount).toBeGreaterThan(0);

    // --- Now switch between experiments and verify data loads ---
    console.log('=== Switching to Experiment 1 ===');
    const exp1Button = page.locator(`button.experiment-button:has-text("${exp1Name}")`);
    await exp1Button.click();
    await page.waitForTimeout(2000);  // Wait for load + render

    const exp1ReloadedNodes = await countTreeNodes(page);
    const exp1ReloadedMessages = await countMessages(page);
    console.log(`Experiment 1 reloaded: ${exp1ReloadedNodes} nodes, ${exp1ReloadedMessages} messages`);
    expect(exp1ReloadedNodes).toBeGreaterThan(0);

    // --- Switch to Experiment 2 ---
    console.log('=== Switching to Experiment 2 ===');
    const exp2Button = page.locator(`button.experiment-button:has-text("${exp2Name}")`);
    await exp2Button.click();
    await page.waitForTimeout(2000);

    const exp2ReloadedNodes = await countTreeNodes(page);
    const exp2ReloadedMessages = await countMessages(page);
    console.log(`Experiment 2 reloaded: ${exp2ReloadedNodes} nodes, ${exp2ReloadedMessages} messages`);
    expect(exp2ReloadedNodes).toBeGreaterThan(0);

    // --- Switch to Experiment 3 ---
    console.log('=== Switching to Experiment 3 ===');
    const exp3Button = page.locator(`button.experiment-button:has-text("${exp3Name}")`);
    await exp3Button.click();
    await page.waitForTimeout(2000);

    const exp3ReloadedNodes = await countTreeNodes(page);
    const exp3ReloadedMessages = await countMessages(page);
    console.log(`Experiment 3 reloaded: ${exp3ReloadedNodes} nodes, ${exp3ReloadedMessages} messages`);
    expect(exp3ReloadedNodes).toBeGreaterThan(0);

    // --- Switch back to Experiment 1 one more time ---
    console.log('=== Switching to Experiment 1 again ===');
    await exp1Button.click();
    await page.waitForTimeout(2000);

    const exp1FinalNodes = await countTreeNodes(page);
    console.log(`Experiment 1 final check: ${exp1FinalNodes} nodes`);
    expect(exp1FinalNodes).toBeGreaterThan(0);

    // --- Print summary ---
    console.log('\n=== SUMMARY ===');
    console.log(`Experiment 1: original=${exp1Result.nodeCount}, reloaded=${exp1ReloadedNodes}, final=${exp1FinalNodes}`);
    console.log(`Experiment 2: original=${exp2Result.nodeCount}, reloaded=${exp2ReloadedNodes}`);
    console.log(`Experiment 3: original=${exp3Result.nodeCount}, reloaded=${exp3ReloadedNodes}`);

    // --- Also verify DB state ---
    // (We check browser console for save confirmations)
    const saveConfirmations = consoleLogs.filter(l => l.includes('Session saved to database') || l.includes('Saving experiments'));
    console.log(`\n=== SAVE LOGS (${saveConfirmations.length} entries) ===`);
    for (const log of saveConfirmations.slice(-10)) {
      console.log(log);
    }

    // Print any warnings
    const warnings = consoleLogs.filter(l => l.includes('[warning]') || l.includes('Failed to save') || l.includes('not found'));
    if (warnings.length > 0) {
      console.log(`\n=== WARNINGS (${warnings.length}) ===`);
      for (const w of warnings) {
        console.log(w);
      }
    }
  });
});
