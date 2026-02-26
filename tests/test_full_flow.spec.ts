import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'http://localhost:5173';
const MOLECULE = 'CC.O=[N+]([O-])C(COc1nc(OCC([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])nc(C([N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-])n1)([N+](=O)[O-])[N+](=O)[O-].[Na+].[N-]=[N+]=[N-].O';

test.describe.serial('Experiment Save+Load', () => {

  test('Phase 1: Create project, run 3 experiments, observe', async ({ page }) => {
    page.on('console', msg => {
      const t = msg.text();
      if (t.includes('Loading context') || t.includes('Project not found') ||
          t.includes('Experiment not found') || t.includes('[DatabaseDataSource]') ||
          t.includes('Saving experiments') || t.includes('Session') ||
          t.includes('Applying restored') || t.includes('computation') ||
          t.includes('save_session') || t.includes('checkpoint') ||
          t.includes('loadStateFromCurrentExperiment')) {
        console.log(`[${msg.type().toUpperCase()}] ${t.substring(0, 200)}`);
      }
    });

    await page.goto(BASE_URL);
    await page.waitForTimeout(5000);

    // Clear localStorage
    await page.evaluate(() => localStorage.clear());
    await page.reload();
    await page.waitForTimeout(3000);

    // Take initial screenshots
    await page.screenshot({ path: '/tmp/phase1_initial.png', fullPage: true });

    // Find the "+" button or "New Project" button in the sidebar
    // The sidebar should be visible on the left
    console.log('\n--- Creating project ---');
    
    // Click the + button for new project (it's typically a small + icon)
    const addProjectBtn = page.locator('button').filter({ hasText: /\+/ }).first();
    const addProjectVisible = await addProjectBtn.isVisible({ timeout: 2000 }).catch(() => false);
    
    if (addProjectVisible) {
      await addProjectBtn.click();
      await page.waitForTimeout(1000);
    } else {
      // Try finding it by aria label or other selectors
      const newProjBtn = page.locator('[title*="project" i], [aria-label*="project" i]').first();
      if (await newProjBtn.isVisible({ timeout: 2000 }).catch(() => false)) {
        await newProjBtn.click();
        await page.waitForTimeout(1000);
      }
    }

    await page.screenshot({ path: '/tmp/phase1_after_project.png', fullPage: true });

    // Check what buttons/elements exist in the sidebar
    const buttons = await page.evaluate(() => {
      return Array.from(document.querySelectorAll('button')).map(b => ({
        text: b.textContent?.trim().substring(0, 50),
        visible: b.offsetParent !== null,
        class: b.className.substring(0, 80),
      }));
    });
    console.log('Visible buttons:', JSON.stringify(buttons.filter(b => b.visible).map(b => b.text)));

    // Check if we have a project now
    const projs = await page.evaluate(async () => {
      const resp = await fetch('http://localhost:8001/api/projects/');
      return resp.json();
    });
    console.log('Projects in DB:', projs.length);

    // Create experiment  
    console.log('\n--- Creating experiment 1 ---');
    const addExpBtn = page.locator('button').filter({ hasText: /experiment|exp/i }).first();
    if (await addExpBtn.isVisible({ timeout: 2000 }).catch(() => false)) {
      await addExpBtn.click();
      await page.waitForTimeout(1000);
    }
    
    // Wait for experiment to be created and DB to be updated
    await page.waitForTimeout(3000);

    // Enter SMILES
    const smilesInput = page.locator('input[type="text"], textarea').first();
    if (await smilesInput.isVisible()) {
      await smilesInput.fill(MOLECULE);
      await page.waitForTimeout(500);
    }

    // Click Run
    console.log('\n--- Running experiment 1 ---');
    const runBtn = page.locator('button').filter({ hasText: 'Run' }).first();
    if (await runBtn.isVisible({ timeout: 2000 }).catch(() => false)) {
      await runBtn.click();
      
      // Wait for computation to complete (Run button becomes enabled again)
      await page.waitForFunction(() => {
        const btns = Array.from(document.querySelectorAll('button'));
        const run = btns.find(b => b.textContent?.includes('Run'));
        return run && !run.disabled;
      }, {}, { timeout: 120000 });
      
      console.log('Experiment 1 complete');
      await page.waitForTimeout(2000);
    }

    // Count nodes
    const exp1Nodes = await page.evaluate(() => document.querySelectorAll('foreignObject').length);
    console.log(`Experiment 1 rendered nodes: ${exp1Nodes}`);

    await page.screenshot({ path: '/tmp/phase1_exp1_done.png', fullPage: true });

    // Create experiment 2
    console.log('\n--- Creating and running experiment 2 ---');
    if (await addExpBtn.isVisible({ timeout: 2000 }).catch(() => false)) {
      await addExpBtn.click();
      await page.waitForTimeout(1000);
    }
    
    await page.waitForTimeout(2000);
    
    const smilesInput2 = page.locator('input[type="text"], textarea').first();
    if (await smilesInput2.isVisible()) {
      await smilesInput2.fill(MOLECULE);
      await page.waitForTimeout(500);
    }

    const runBtn2 = page.locator('button').filter({ hasText: 'Run' }).first();
    if (await runBtn2.isVisible({ timeout: 2000 }).catch(() => false)) {
      await runBtn2.click();
      
      await page.waitForFunction(() => {
        const btns = Array.from(document.querySelectorAll('button'));
        const run = btns.find(b => b.textContent?.includes('Run'));
        return run && !run.disabled;
      }, {}, { timeout: 120000 });
      
      console.log('Experiment 2 complete');
      await page.waitForTimeout(2000);
    }

    const exp2Nodes = await page.evaluate(() => document.querySelectorAll('foreignObject').length);
    console.log(`Experiment 2 rendered nodes: ${exp2Nodes}`);

    await page.screenshot({ path: '/tmp/phase1_exp2_done.png', fullPage: true });

    // Create experiment 3
    console.log('\n--- Creating and running experiment 3 ---');
    if (await addExpBtn.isVisible({ timeout: 2000 }).catch(() => false)) {
      await addExpBtn.click();
      await page.waitForTimeout(1000);
    }
    
    await page.waitForTimeout(2000);
    
    const smilesInput3 = page.locator('input[type="text"], textarea').first();
    if (await smilesInput3.isVisible()) {
      await smilesInput3.fill(MOLECULE);
      await page.waitForTimeout(500);
    }

    const runBtn3 = page.locator('button').filter({ hasText: 'Run' }).first();
    if (await runBtn3.isVisible({ timeout: 2000 }).catch(() => false)) {
      await runBtn3.click();
      
      await page.waitForFunction(() => {
        const btns = Array.from(document.querySelectorAll('button'));
        const run = btns.find(b => b.textContent?.includes('Run'));
        return run && !run.disabled;
      }, {}, { timeout: 120000 });
      
      console.log('Experiment 3 complete');
      await page.waitForTimeout(2000);
    }

    const exp3Nodes = await page.evaluate(() => document.querySelectorAll('foreignObject').length);
    console.log(`Experiment 3 rendered nodes: ${exp3Nodes}`);

    await page.screenshot({ path: '/tmp/phase1_exp3_done.png', fullPage: true });

    // Now click experiment 1 and check if it loads
    console.log('\n--- Clicking back to Experiment 1 ---');
    const exp1Link = page.locator('text="Experiment 1"').first();
    if (await exp1Link.isVisible({ timeout: 3000 }).catch(() => false)) {
      await exp1Link.click();
      await page.waitForTimeout(3000);
      
      const nodesAfterClick1 = await page.evaluate(() => document.querySelectorAll('foreignObject').length);
      console.log(`After clicking Exp1: ${nodesAfterClick1} nodes`);
      await page.screenshot({ path: '/tmp/phase1_click_exp1.png', fullPage: true });
    }

    // Click experiment 2
    console.log('\n--- Clicking Experiment 2 ---');
    const exp2Link = page.locator('text="Experiment 2"').first();
    if (await exp2Link.isVisible({ timeout: 3000 }).catch(() => false)) {
      await exp2Link.click();
      await page.waitForTimeout(3000);
      
      const nodesAfterClick2 = await page.evaluate(() => document.querySelectorAll('foreignObject').length);
      console.log(`After clicking Exp2: ${nodesAfterClick2} nodes`);
      await page.screenshot({ path: '/tmp/phase1_click_exp2.png', fullPage: true });
    }

    // Click experiment 3
    console.log('\n--- Clicking Experiment 3 ---');
    const exp3Link = page.locator('text="Experiment 3"').first();
    if (await exp3Link.isVisible({ timeout: 3000 }).catch(() => false)) {
      await exp3Link.click();
      await page.waitForTimeout(3000);
      
      const nodesAfterClick3 = await page.evaluate(() => document.querySelectorAll('foreignObject').length);
      console.log(`After clicking Exp3: ${nodesAfterClick3} nodes`);
      await page.screenshot({ path: '/tmp/phase1_click_exp3.png', fullPage: true });
    }

    // Check DB state at the end
    console.log('\n--- Final DB state ---');
    const finalProjs = await page.evaluate(async () => {
      const resp = await fetch('http://localhost:8001/api/projects/');
      return resp.json();
    });
    for (const proj of finalProjs) {
      for (const exp of proj.experiments) {
        const nodes = exp.treeNodes?.length ?? 0;
        const sidebar = exp.sidebarState?.messages?.length ?? 0;
        console.log(`  ${exp.name}: ${nodes} nodes, ${sidebar} sidebar msgs`);
      }
    }
  });
});
