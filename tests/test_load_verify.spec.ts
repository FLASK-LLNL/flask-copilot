import { test, expect, Page } from '@playwright/test';

const BASE_URL = 'http://localhost:5173';

test.describe('Experiment Load Tests', () => {
  test('Load page and verify experiments render from DB', async ({ page }) => {
    // Enable console logging
    page.on('console', msg => {
      const text = msg.text();
      if (text.includes('Loading context') || 
          text.includes('loadStateFromCurrentExperiment') ||
          text.includes('Experiment not found') ||
          text.includes('Project not found') ||
          text.includes('Session loaded') ||
          text.includes('Applying restored') ||
          text.includes('[DatabaseDataSource]') ||
          text.includes('ERROR') ||
          text.includes('Session restored')) {
        console.log(`[BROWSER ${msg.type()}] ${text}`);
      }
    });

    // Go to the app
    await page.goto(BASE_URL);
    
    // Wait for page to fully load and WS to connect
    await page.waitForTimeout(5000);

    // Check what's in projectsRef by evaluating in page context
    // The data comes from the 2s poll of /api/projects
    const apiData = await page.evaluate(async () => {
      const resp = await fetch('http://localhost:8001/api/projects/', {
        headers: { 'X-Forwarded-User': 'marathe1' }
      });
      return await resp.json();
    });
    
    console.log('\n=== API Response ===');
    for (const proj of apiData) {
      console.log(`Project: ${proj.name}`);
      for (const exp of proj.experiments) {
        const nodes = exp.treeNodes?.length ?? 0;
        const sidebar = exp.sidebarState?.messages?.length ?? 0;
        console.log(`  ${exp.name}: ${nodes} nodes, ${sidebar} sidebar msgs`);
      }
    }

    // Wait for project data to load in the app
    await page.waitForTimeout(3000);

    // Check what the app's projectsRef has
    // This would require inspecting React internals, but let's check what's
    // visually rendered instead
    
    // Take a screenshot of the initial state
    await page.screenshot({ path: '/tmp/test_initial.png' });

    // Check how many tree nodes are rendered in the SVG/canvas
    const initialNodeCount = await page.evaluate(() => {
      // Count rendered tree nodes -- they're foreignObject elements in the SVG
      const foreignObjects = document.querySelectorAll('foreignObject');
      return foreignObjects.length;
    });
    console.log(`\nInitial rendered nodes: ${initialNodeCount}`);

    // Check if sidebar has messages
    const initialSidebarMsgs = await page.evaluate(() => {
      const msgs = document.querySelectorAll('[class*="message"]');
      return msgs.length;
    });
    console.log(`Initial sidebar messages: ${initialSidebarMsgs}`);

    // Now find experiment items in the sidebar and click them
    // First, let's examine the sidebar structure
    const sidebarHTML = await page.evaluate(() => {
      const sidebar = document.querySelector('[class*="projectSidebar"], [class*="sidebar"]');
      return sidebar?.innerHTML?.substring(0, 2000) || 'sidebar not found';
    });
    console.log(`\nSidebar HTML (first 500 chars): ${sidebarHTML.substring(0, 500)}`);

    // Get all experiment list items
    const experiments = await page.evaluate(() => {
      // Look for experiment items -- they usually have experiment name text
      const allElements = document.querySelectorAll('*');
      const expElements: string[] = [];
      for (const el of allElements) {
        const text = el.textContent?.trim() || '';
        if (text.match(/^Experiment \d+$/) && el.children.length === 0) {
          expElements.push(text);
        }
      }
      return [...new Set(expElements)];
    });
    console.log(`\nExperiment labels found: ${JSON.stringify(experiments)}`);

    // Click each experiment and check if nodes render
    for (const expName of experiments) {
      console.log(`\n--- Clicking ${expName} ---`);
      
      // Click the element containing this text
      const expLocator = page.locator(`text="${expName}"`).first();
      if (await expLocator.isVisible({ timeout: 2000 })) {
        await expLocator.click();
        await page.waitForTimeout(3000); // Wait for load + render
        
        const nodeCount = await page.evaluate(() => {
          const foreignObjects = document.querySelectorAll('foreignObject');
          return foreignObjects.length;
        });
        console.log(`  Rendered nodes after click: ${nodeCount}`);
        
        // Take screenshot
        await page.screenshot({ path: `/tmp/test_${expName.replace(' ', '_')}.png` });
        
        if (nodeCount === 0) {
          console.log(`  WARNING: ${expName} has 0 rendered nodes!`);
        }
      } else {
        console.log(`  Could not find clickable element for ${expName}`);
      }
    }

    // Final verification
    console.log('\n=== DONE ===');
  });
});
