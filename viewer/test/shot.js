// Renders the viewer with real docking data in headless Chromium and writes screenshots + a console log.
//   node viewer/test/shot.js <outdir> [dark|light]
const { chromium } = require('playwright');
const fs = require('fs'), path = require('path');
const root = path.resolve(__dirname, '..');
const data = (f) => fs.readFileSync(path.join(__dirname, 'data', f));
const b64 = (buf) => Buffer.from(buf).toString('base64');

(async () => {
  const out = process.argv[2] || '/tmp'; const theme = process.argv[3] || 'dark';
  fs.mkdirSync(out, { recursive: true });
  const browser = await chromium.launch({ args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist', '--allow-file-access-from-files'] });
  const page = await browser.newPage({ viewport: { width: 1100, height: 720 }, deviceScaleFactor: 1 });
  const log = [];
  page.on('console', (m) => log.push(`[${m.type()}] ${m.text()}`));
  page.on('pageerror', (e) => log.push('[pageerror] ' + e.message));
  await page.addInitScript(() => { window.__events = []; window.addEventListener('dock-event', (e) => window.__events.push(e.detail)); });
  await page.goto('file://' + path.join(root, 'index.html'));
  await page.waitForFunction(() => window.__events.some((e) => e.type === 'ready' || e.type === 'error'), null, { timeout: 60000 });
  const events = () => page.evaluate(() => window.__events);
  await page.screenshot({ path: path.join(out, '0-empty.png'), timeout: 120000 }); console.log('shot', '0-empty.png');

  const inter = JSON.parse(data('interactions.json'));
  const send = async (cmd) => { log.push('--- ' + cmd.cmd); await page.evaluate((c) => window.dock.command(c), cmd); };
  await send({ cmd: 'setTheme', theme });
  await send({ cmd: 'loadReceptor', data: b64(data('receptor.pdb')), format: 'pdb', name: '1STP · streptavidin' });
  await page.waitForTimeout(1500);
  await page.screenshot({ path: path.join(out, '1-receptor.png'), timeout: 120000 }); console.log('shot', '1-receptor.png');
  await send({ cmd: 'setBox', center: inter.box.center, size: inter.box.size });
  const poses = [1, 2, 3, 4, 5, 6, 7, 8, 9].map((i) => ({ index: i, score: -4.6 + i * 0.03 }));
  await send({ cmd: 'loadPose', data: b64(data('poses.sdf')), format: 'sdf', index: 1, score: -4.575, poses });
  await send({ cmd: 'setInteractions', interactions: inter['1'].interactions });
  await page.waitForTimeout(2500);
  await page.screenshot({ path: path.join(out, '2-complex.png'), timeout: 120000 }); console.log('shot', '2-complex.png');
  await send({ cmd: 'setOptions', surface: true });
  await page.waitForTimeout(2000);
  await page.screenshot({ path: path.join(out, '3-surface.png'), timeout: 120000 }); console.log('shot', '3-surface.png');
  const evs = await events();
  fs.writeFileSync(path.join(out, 'console.log'), log.join('\n') + '\n' + JSON.stringify(evs));
  await browser.close();
  const bad = log.filter((l) => /^\[(error|pageerror)\]/.test(l)).concat(evs.filter((e) => e.type === 'error').map((e) => 'event: ' + e.message));
  if (bad.length) { console.error('FAILED:\n' + bad.join('\n')); process.exit(2); }
  console.log('ok: no errors; events:', evs.map((e) => e.type).join(','));
})().catch((e) => { console.error(e); process.exit(1); });
