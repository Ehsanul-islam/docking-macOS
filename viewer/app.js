/* Dock 3D viewer: a thin, styled shell around Mol* (MIT). The scene is described declaratively (MolViewSpec) and
 * rebuilt from `state` on every change, so Swift only ever sends small JSON commands (see window.dock.command).
 * Events go back through webkit.messageHandlers.dock when hosted in a WKWebView, else to window `dock-event`. */
(() => {
  'use strict';

  const THEMES = {
    dark: {
      palette: ['#6ea8fe', '#7f9cf7', '#9a8ff0', '#b886ea', '#d884da', '#ee88bd', '#f69aa0'],
      ligandCarbon: '#ffc15e', cofactorCarbon: '#8fa1b8', pocketCarbon: '#aab6c8',
      surface: '#8aa4d6', box: '#5eead4', label: '#e8ecf6', labelBg: '#0d1220',
    },
    light: {
      palette: ['#2f6fe0', '#4a63d8', '#6a55cc', '#8c4fbf', '#b04fa8', '#cf5a8c', '#e0687a'],
      ligandCarbon: '#e0850b', cofactorCarbon: '#6b7a90', pocketCarbon: '#7d8aa0',
      surface: '#6f8fc9', box: '#0f9d8a', label: '#1c2438', labelBg: '#ffffff',
    },
  };
  const INTERACTION_STYLE = [
    // prefix match on the helper's type names
    { match: 'Hydrogen bond', key: 'hbond',   name: 'Hydrogen bond', color: '#4da3ff', dash: 0.18, labelled: false },
    { match: 'Salt bridge',   key: 'salt',    name: 'Salt bridge',   color: '#ff5d7a', dash: 0.18, labelled: false },
    { match: 'π-stacking',    key: 'pi',      name: 'π-stacking',    color: '#b57cff', dash: 0.25 },
    { match: 'π-cation',      key: 'cation',  name: 'π-cation',      color: '#ff9f43', dash: 0.25 },
    { match: 'Halogen',       key: 'halogen', name: 'Halogen bond',  color: '#2ed3b7', dash: 0.18 },
    { match: 'Hydrophobic',   key: 'hydro',   name: 'Hydrophobic',   color: '#9aa7b8', dash: 0.35 },
  ];
  const styleFor = (type) => INTERACTION_STYLE.find((s) => type.startsWith(s.match));

  const state = {
    theme: 'dark',
    receptor: null,        // { text, format, name }
    pose: null,            // { text, format, index, score }
    poses: [],             // [{ index, score }]
    box: null,             // { center:[x,y,z], size:[x,y,z] }
    interactions: [],      // helper's interaction records (3D coordinates)
    options: { cartoon: true, surface: false, pocket: true, interactions: true, labels: true, box: true, spin: false },
    framed: false,
    focusPending: false,
  };

  let viewer = null;
  let renderToken = 0;
  let rendering = false;
  const urls = [];

  // ---------------------------------------------------------------- events out
  function post(type, payload = {}) {
    const msg = { type, ...payload };
    try {
      if (window.webkit && window.webkit.messageHandlers && window.webkit.messageHandlers.dock) {
        window.webkit.messageHandlers.dock.postMessage(msg);
      }
    } catch (_) { /* hosted outside WKWebView */ }
    window.dispatchEvent(new CustomEvent('dock-event', { detail: msg }));
  }

  // ---------------------------------------------------------------- helpers
  const $ = (id) => document.getElementById(id);
  const blobUrl = (text, type = 'text/plain') => {
    const u = URL.createObjectURL(new Blob([text], { type }));
    urls.push(u);
    return u;
  };
  function freeUrls() { while (urls.length) URL.revokeObjectURL(urls.pop()); }
  function busy(on, text) { $('busy').classList.toggle('hidden', !on); if (text) $('busy-text').textContent = text; }
  const decode = (b64) => new TextDecoder().decode(Uint8Array.from(atob(b64), (c) => c.charCodeAt(0)));
  const hex = (s) => parseInt(s.slice(1), 16);

  // PDB parsing: only what the shell needs (residue centres, pocket selection)
  function parsePdb(text) {
    const atoms = [];
    for (const line of text.split('\n')) {
      if (!(line.startsWith('ATOM') || line.startsWith('HETATM'))) continue;
      const el = (line.slice(76, 78).trim() || line.slice(12, 16).replace(/[0-9\s]/g, '')[0] || '').toUpperCase();
      if (el === 'H' || el === 'D') continue;
      atoms.push({ chain: line[21], resi: parseInt(line.slice(22, 26), 10), icode: line[26].trim(), resn: line.slice(17, 20).trim(),
                   name: line.slice(12, 16).trim(), het: line.startsWith('HETATM'),
                   x: +line.slice(30, 38), y: +line.slice(38, 46), z: +line.slice(46, 54) });
    }
    return atoms;
  }
  function parseMolblockAtoms(text) {
    const lines = text.split('\n');
    const n = parseInt(lines[3].slice(0, 3), 10);
    const out = [];
    for (let i = 0; i < n; i++) {
      const l = lines[4 + i];
      out.push({ x: +l.slice(0, 10), y: +l.slice(10, 20), z: +l.slice(20, 30), el: l.slice(31, 34).trim() });
    }
    return out;
  }

  // A pose as mmCIF with explicit bond orders (Mol* reads _chem_comp_bond), so double and aromatic bonds are drawn.
  function molblockToCif(text) {
    const lines = text.split('\n');
    const na = parseInt(lines[3].slice(0, 3), 10), nb = parseInt(lines[3].slice(3, 6), 10);
    const count = {};
    const names = [];
    const rows = [];
    for (let i = 0; i < na; i++) {
      const l = lines[4 + i];
      const el = l.slice(31, 34).trim();
      count[el] = (count[el] || 0) + 1;
      const name = el.toUpperCase() + count[el];
      names.push(name);
      rows.push(`HETATM ${i + 1} ${el} ${name} LIG L 1 1 ${(+l.slice(0, 10)).toFixed(3)} ${(+l.slice(10, 20)).toFixed(3)} ${(+l.slice(20, 30)).toFixed(3)} 1 0 1 L 1`);
    }
    const orders = { 1: 'SING', 2: 'DOUB', 3: 'TRIP', 4: 'AROM' };
    const bonds = [];
    for (let i = 0; i < nb; i++) {
      const l = lines[4 + na + i];
      const a = parseInt(l.slice(0, 3), 10) - 1, b = parseInt(l.slice(3, 6), 10) - 1, t = parseInt(l.slice(6, 9), 10);
      bonds.push(`LIG ${names[a]} ${names[b]} ${orders[t] || 'SING'} ${t === 4 ? 'Y' : 'N'}`);
    }
    return ['data_LIG', '#', 'loop_', '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol', '_atom_site.label_atom_id',
      '_atom_site.label_comp_id', '_atom_site.label_asym_id', '_atom_site.label_entity_id', '_atom_site.label_seq_id',
      '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy', '_atom_site.B_iso_or_equiv',
      '_atom_site.auth_seq_id', '_atom_site.auth_asym_id', '_atom_site.pdbx_PDB_model_num', ...rows.map((r) => r.replace(/ 1 1 (-?\d)/, ' 1 . $1')),
      '#', 'loop_', '_chem_comp_bond.comp_id', '_chem_comp_bond.atom_id_1', '_chem_comp_bond.atom_id_2', '_chem_comp_bond.value_order',
      '_chem_comp_bond.pdbx_aromatic_flag', ...bonds, '#', ''].join('\n');
  }

  function poseAtoms() {
    if (!state.pose) return [];
    try {
      if (state.pose.format === 'sdf' || state.pose.format === 'mol') return parseMolblockAtoms(state.pose.text);
    } catch (_) { /* fall through */ }
    return [];
  }
  function pocketResidues(radius = 4.5) {
    if (!state.receptor) return [];
    const lig = poseAtoms();
    if (!lig.length) return [];
    const seen = new Map();
    const r2 = radius * radius;
    for (const a of parsePdb(state.receptor.text)) {
      if (a.het) continue;
      const key = `${a.chain}:${a.resi}${a.icode}`;
      if (seen.has(key)) continue;
      for (const l of lig) {
        const dx = a.x - l.x, dy = a.y - l.y, dz = a.z - l.z;
        if (dx * dx + dy * dy + dz * dz <= r2) { seen.set(key, { chain: a.chain, resi: a.resi, resn: a.resn }); break; }
      }
    }
    return [...seen.values()];
  }

  // Gradient colouring N -> C per chain, supplied to Mol* as an annotation file.
  function gradientAnnotation(palette) {
    const atoms = parsePdb(state.receptor.text).filter((a) => !a.het && a.name === 'CA');
    const byChain = new Map();
    for (const a of atoms) { if (!byChain.has(a.chain)) byChain.set(a.chain, []); byChain.get(a.chain).push(a); }
    const rows = [];
    const stops = palette.map((c) => [parseInt(c.slice(1, 3), 16), parseInt(c.slice(3, 5), 16), parseInt(c.slice(5, 7), 16)]);
    for (const [chain, list] of byChain) {
      list.forEach((a, i) => {
        const t = list.length > 1 ? i / (list.length - 1) : 0;
        const p = t * (stops.length - 1), k = Math.min(Math.floor(p), stops.length - 2), f = p - k;
        const c = stops[k].map((v, j) => Math.round(v + (stops[k + 1][j] - v) * f));
        rows.push({ auth_asym_id: chain, auth_seq_id: a.resi, color: '#' + c.map((v) => v.toString(16).padStart(2, '0')).join('') });
      });
    }
    return rows;
  }

  // ---------------------------------------------------------------- scene
  function buildMvs(cam) {
    const th = THEMES[state.theme];
    const o = state.options;
    const root = [];
    const sel = (r) => ({ auth_asym_id: r.chain, auth_seq_id: r.resi });

    if (state.receptor) {
      const rec = { kind: 'download', params: { url: blobUrl(state.receptor.text) }, children: [] };
      const struct = { kind: 'structure', params: { type: 'model' }, children: [] };
      rec.children.push({ kind: 'parse', params: { format: state.receptor.format || 'pdb' }, children: [struct] });
      root.push(rec);

      if (o.cartoon) {
        const rows = gradientAnnotation(th.palette);
        const comp = { kind: 'component', params: { selector: 'polymer' }, children: [
          { kind: 'representation', params: { type: 'cartoon' }, children: [
            { kind: 'color_from_uri', params: { uri: blobUrl(JSON.stringify(rows), 'application/json'), format: 'json', schema: 'auth_residue' } },
          ] },
        ] };
        struct.children.push(comp);
      }
      const pocket = pocketResidues();
      if (o.surface && pocket.length) {
        struct.children.push({ kind: 'component', params: { selector: pocket.map(sel) }, children: [
          { kind: 'representation', params: { type: 'surface', surface_type: 'molecular', ignore_hydrogens: true }, children: [
            { kind: 'color', params: { color: th.surface } },
            { kind: 'opacity', params: { opacity: 0.42 } },
          ] },
        ] });
      }
      if (o.pocket && pocket.length) {
        struct.children.push({ kind: 'component', params: { selector: pocket.map(sel) }, children: [
          { kind: 'representation', params: { type: 'ball_and_stick', size_factor: 0.42, ignore_hydrogens: true }, children: [
            { kind: 'color', params: { color: th.pocketCarbon } },
            { kind: 'color', params: { selector: { type_symbol: 'N' }, color: '#5b8def' } },
            { kind: 'color', params: { selector: { type_symbol: 'O' }, color: '#ef5b5b' } },
            { kind: 'color', params: { selector: { type_symbol: 'S' }, color: '#e6c34a' } },
          ] },
        ] });
      }
      // cofactors / ions kept in the receptor
      struct.children.push({ kind: 'component', params: { selector: 'ligand' }, children: [
        { kind: 'representation', params: { type: 'ball_and_stick', size_factor: 0.55, ignore_hydrogens: true }, children: [
          { kind: 'color', params: { color: th.cofactorCarbon } },
        ] },
      ] });
      struct.children.push({ kind: 'component', params: { selector: 'ion' }, children: [
        { kind: 'representation', params: { type: 'ball_and_stick' }, children: [{ kind: 'color', params: { color: '#c0c8d8' } }] },
      ] });
    }

    if (state.pose) {
      const poseIsMol = state.pose.format === 'sdf' || state.pose.format === 'mol';
      const lig = { kind: 'download', params: { url: blobUrl(poseIsMol ? molblockToCif(state.pose.text) : state.pose.text) }, children: [
        { kind: 'parse', params: { format: poseIsMol ? 'mmcif' : state.pose.format }, children: [
          { kind: 'structure', params: { type: 'model' }, children: [
            { kind: 'component', params: { selector: 'all' }, children: [
              ...(cam.focus ? [{ kind: 'focus', params: { radius_factor: 1.6 } }] : []),
              { kind: 'representation', params: { type: 'ball_and_stick', size_factor: 0.62, ignore_hydrogens: true }, children: [
                { kind: 'color', params: { color: th.ligandCarbon } },
                { kind: 'color', params: { selector: { type_symbol: 'N' }, color: '#5b8def' } },
                { kind: 'color', params: { selector: { type_symbol: 'O' }, color: '#ff5d5d' } },
                { kind: 'color', params: { selector: { type_symbol: 'S' }, color: '#f2cf4c' } },
                { kind: 'color', params: { selector: { type_symbol: 'F' }, color: '#7be0a5' } },
                { kind: 'color', params: { selector: { type_symbol: 'Cl' }, color: '#5fd35f' } },
              ] },
            ] },
          ] },
        ] },
      ] };
      root.push(lig);
    }

    const prim = [];
    const boxGroups = [];
    if (o.box && state.box) {
      const base = { kind: 'box', center: state.box.center, extent: state.box.size, tooltip: 'Search box' };
      boxGroups.push({ kind: 'primitives', params: { color: th.box, opacity: 0.07 }, children: [
        { kind: 'primitive', params: { ...base, show_faces: true, show_edges: false } }] });
      boxGroups.push({ kind: 'primitives', params: { color: th.box, opacity: 0.9 }, children: [
        { kind: 'primitive', params: { ...base, show_faces: false, show_edges: true, edge_radius: 0.12 } }] });
    }
    if (o.interactions) {
      for (const i of state.interactions) {
        const st = styleFor(i.type);
        if (!st) continue;
        prim.push({ kind: 'primitive', params: { kind: 'distance_measurement', start: i.lig, end: i.prot, color: st.color,
                    radius: 0.11, dash_length: st.dash, label_template: '{{distance}} Å', label_size: st.labelled ? 0.75 : 0.001, label_color: th.label } });
      }
    }
    if (o.labels) {
      const seen = new Set();
      for (const i of state.interactions) {
        if (seen.has(i.residue)) continue;
        seen.add(i.residue);
        prim.push({ kind: 'primitive', params: { kind: 'label', position: i.prot, text: i.residue, label_size: 1.15,
                    label_color: th.label, label_offset: 2.2 } });
      }
    }
    boxGroups.forEach((g) => root.push(g));
    if (prim.length) root.push({ kind: 'primitives', params: { color: th.box, opacity: 1 }, children: prim });

    if (cam.restore) root.push({ kind: 'camera', params: { target: cam.restore.target, position: cam.restore.position, up: cam.restore.up } });
    return { metadata: { version: '1', timestamp: new Date().toISOString(), title: 'Dock' }, root: { kind: 'root', children: root } };
  }

  async function configureCanvas() {
    const dark = state.theme === 'dark';
    const c3d = viewer.plugin.canvas3d;
    try {
      c3d.setProps({
        transparentBackground: true,
        renderer: { backgroundColor: dark ? 0x0d1220 : 0xeef1f8, ambientIntensity: dark ? 0.55 : 0.7, light: undefined,
                    highlightColor: 0x5eead4, selectColor: 0xffffff, },
        postprocessing: {
          occlusion: { name: 'on', params: { multiScale: { name: 'off', params: {} }, radius: 3.5, bias: 0.7, blurKernelSize: 15, samples: 32, resolutionScale: 1, color: 0x000000, transparentThreshold: 0.4 } },
          outline: { name: 'on', params: { scale: 1, threshold: 0.2, color: dark ? 0x05080f : 0x1c2438, includeTransparent: false } },
          shadow: { name: 'off', params: {} },
          antialiasing: { name: 'smaa', params: { edgeThreshold: 0.1, maxSearchSteps: 16 } },
        },
        cameraFog: { name: 'on', params: { intensity: 5 } },
        camera: { helper: { axes: { name: 'off', params: {} } } },
      });
    } catch (e) { post('warning', { message: 'canvas settings: ' + e.message }); }
  }

  async function render() {
    if (!viewer) return;
    const token = ++renderToken;
    const hasContent = !!(state.receptor || state.pose);
    $('empty').classList.toggle('hidden', hasContent);
    refreshChrome();
    if (!hasContent) { try { await viewer.plugin.clear(); } catch (_) {} return; }
    busy(true, 'Building the scene…');
    freeUrls();
    rendering = true;
    const snapshot = state.framed ? state.cameraSnap : null;
    try {
      await viewer.plugin.clear();
      const cam = { focus: state.focusPending && !!state.pose, restore: snapshot && !state.focusPending ? snapshot : null };
      await viewer.loadMvsData(JSON.stringify(buildMvs(cam)), 'mvsj', { sanityChecks: false });
      if (token !== renderToken) return;
      state.focusPending = false;
      state.framed = true;
      await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
      state.cameraSnap = viewer.plugin.canvas3d.camera.getSnapshot();
      setSpin(state.options.spin);
      post('rendered', { receptor: !!state.receptor, pose: state.pose ? state.pose.index : null });
    } catch (e) {
      post('error', { message: 'The 3D view could not be built: ' + (e && e.message ? e.message : e), stack: e && e.stack });
    } finally {
      if (token === renderToken) { busy(false); rendering = false; }
    }
  }

  async function focusLigand(ms = 400) {
    const c3d = viewer.plugin.canvas3d;
    const lig = poseAtoms();
    if (!lig.length) { c3d.requestCameraReset({ durationMs: ms }); return; }
    const m = lig.reduce((a, p) => [a[0] + p.x / lig.length, a[1] + p.y / lig.length, a[2] + p.z / lig.length], [0, 0, 0]);
    const r = Math.max(...lig.map((p) => Math.hypot(p.x - m[0], p.y - m[1], p.z - m[2]))) + 9;
    c3d.camera.focus(m, r, ms);
  }

  function setSpin(on) {
    const c3d = viewer && viewer.plugin.canvas3d;
    if (!c3d) return;
    try { c3d.setProps({ trackball: { animate: on ? { name: 'spin', params: { speed: 0.25 } } : { name: 'off', params: {} } } }); } catch (_) {}
  }

  // ---------------------------------------------------------------- chrome (toolbar, legend, pose chips)
  function refreshChrome() {
    document.documentElement.dataset.theme = state.theme;
    document.querySelectorAll('[data-opt]').forEach((b) => b.classList.toggle('on', !!state.options[b.dataset.opt]));
    const has = {
      cartoon: !!state.receptor, surface: !!state.pose, pocket: !!state.pose, interactions: state.interactions.length > 0,
      labels: state.interactions.length > 0, box: !!state.box,
    };
    document.querySelectorAll('[data-opt]').forEach((b) => { if (b.dataset.opt in has) b.classList.toggle('off-state', !has[b.dataset.opt]); });

    const chip = $('chips');
    const title = state.receptor ? state.receptor.name : '';
    chip.classList.toggle('hidden', !title);
    $('chip-title').textContent = title;
    $('chip-score').textContent = state.pose && state.pose.score != null ? `pose ${state.pose.index} · ${state.pose.score.toFixed(2)} kcal/mol` : '';

    const legend = $('legend');
    const counts = new Map();
    for (const i of state.interactions) { const s = styleFor(i.type); if (s) counts.set(s.key, (counts.get(s.key) || 0) + 1); }
    legend.innerHTML = '';
    for (const s of INTERACTION_STYLE) {
      if (!counts.has(s.key)) continue;
      const row = document.createElement('div'); row.className = 'row';
      row.innerHTML = `<span class="swatch" style="border-color:${s.color}"></span><span>${s.name}</span><span class="count">${counts.get(s.key)}</span>`;
      legend.appendChild(row);
    }
    legend.classList.toggle('hidden', !counts.size || !state.options.interactions);

    const poses = $('poses');
    poses.innerHTML = '';
    for (const p of state.poses) {
      const b = document.createElement('button');
      b.textContent = p.index; b.title = p.score != null ? `${p.score.toFixed(2)} kcal/mol` : '';
      b.className = state.pose && state.pose.index === p.index ? 'on' : '';
      b.setAttribute('role', 'tab');
      b.onclick = () => post('poseSelected', { index: p.index });
      poses.appendChild(b);
    }
    poses.classList.toggle('hidden', state.poses.length < 2);
  }

  // ---------------------------------------------------------------- commands in
  const commands = {
    setTheme: ({ theme }) => { state.theme = theme === 'light' ? 'light' : 'dark'; configureCanvas(); return render(); },
    loadReceptor: ({ data, format = 'pdb', name = '' }) => {
      state.receptor = { text: decode(data), format, name }; state.framed = false;
      state.pose = null; state.poses = []; state.interactions = []; state.box = state.box;
      return render();
    },
    loadPose: ({ data, format = 'sdf', index = 1, score = null, poses = null }) => {
      let text = decode(data);
      if (format === 'sdf' && text.includes('$$$$')) {       // a multi-pose file: show just the requested pose
        const blocks = text.split('$$$$').map((b, i) => (i === 0 ? b : b.replace(/^\r?\n/, ''))).filter((b) => b.trim());
        text = blocks[Math.min(Math.max(index, 1), blocks.length) - 1] + '\n$$$$\n';
      }
      if (!state.pose) state.focusPending = true;
      state.pose = { text, format, index, score };
      if (poses) state.poses = poses;
      return render();
    },
    setInteractions: ({ interactions }) => { state.interactions = interactions || []; return render(); },
    setBox: ({ center, size }) => { state.box = center && size ? { center, size } : null; return render(); },
    setOptions: (opts) => { Object.assign(state.options, opts); return render(); },
    focus: ({ target = 'ligand' }) => (target === 'ligand' && state.pose ? focusLigand() : viewer.plugin.canvas3d.requestCameraReset({ durationMs: 400 })),
    clear: () => { state.receptor = null; state.pose = null; state.poses = []; state.interactions = []; state.box = null; state.framed = false; return render(); },
    screenshot: async () => {
      const png = await viewer.plugin.helpers.viewportScreenshot?.getImageDataUri?.();
      post('screenshot', { dataUri: png || null });
    },
  };

  window.dock = {
    state,
    command(cmd) {
      const run = typeof cmd === 'string' ? JSON.parse(cmd) : cmd;
      const fn = commands[run.cmd];
      if (!fn) { post('error', { message: 'unknown viewer command ' + run.cmd }); return Promise.resolve(); }
      return Promise.resolve(fn(run)).catch((e) => post('error', { message: String(e && e.message || e) }));
    },
    ready: null,
  };

  // ---------------------------------------------------------------- toolbar wiring
  document.querySelectorAll('[data-opt]').forEach((b) => b.addEventListener('click', () => {
    const k = b.dataset.opt;
    if (k === 'spin') { state.options.spin = !state.options.spin; setSpin(state.options.spin); refreshChrome(); return; }
    state.options[k] = !state.options[k]; render();
  }));
  document.querySelectorAll('[data-act]').forEach((b) => b.addEventListener('click', () => {
    const a = b.dataset.act;
    if (a === 'focus') commands.focus({});
    if (a === 'reset') viewer.plugin.canvas3d.requestCameraReset({ durationMs: 400 });
    if (a === 'shot') commands.screenshot();
    if (a === 'theme') commands.setTheme({ theme: state.theme === 'dark' ? 'light' : 'dark' });
  }));

  // ---------------------------------------------------------------- boot
  async function boot() {
    if (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches) state.theme = 'light';
    refreshChrome();
    try {
      viewer = await molstar.Viewer.create('canvas', {
        layoutIsExpanded: false, layoutShowControls: false, layoutShowRemoteState: false, layoutShowSequence: false,
        layoutShowLog: false, layoutShowLeftPanel: false, viewportShowExpand: false, viewportShowControls: false,
        viewportShowSelectionMode: false, viewportShowAnimation: false, viewportShowSettings: false, viewportShowTrajectoryControls: false,
        pixelScale: Math.min(window.devicePixelRatio || 1, 2), collapseLeftPanel: true, collapseRightPanel: true,
        pdbProvider: 'rcsb', emdbProvider: 'rcsb',
      });
      await configureCanvas();
      viewer.plugin.canvas3d.didDraw.subscribe(() => { if (!rendering) state.cameraSnap = viewer.plugin.canvas3d.camera.getSnapshot(); });
      window.addEventListener('resize', () => viewer.handleResize());
      post('ready', { molstar: molstar.version });
    } catch (e) {
      post('error', { message: 'The 3D engine could not start: ' + (e && e.message || e) });
    }
  }
  boot();
})();
