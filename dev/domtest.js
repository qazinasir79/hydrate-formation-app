// Headless DOM smoke test for hydrate_app/index.html — run: node domtest.js
const fs = require('fs');
const path = require('path');
class El {
  constructor(tag, id) {
    this.tag = tag; this.id = id || ''; this.children = []; this.dataset = {};
    this.style = {}; this.classList = {
      _s: new Set(),
      add: c => this.classList._s.add(c),
      remove: c => this.classList._s.delete(c),
      toggle: (c, f) => { f === undefined ? (this.classList._s.has(c) ? this.classList._s.delete(c) : this.classList._s.add(c)) : (f ? this.classList._s.add(c) : this.classList._s.delete(c)); },
      contains: c => this.classList._s.has(c)
    };
    this._html = ''; this._text = ''; this.value = ''; this.checked = false;
    this.selectedIndex = 0; this.options = []; this.handlers = {};
    this.clientWidth = 860; this.disabled = false;
  }
  get firstElementChild() { return this.children[0] || null; }
  set innerHTML(v) { this._html = v; } get innerHTML() { return this._html; }
  set textContent(v) { this._text = v; } get textContent() { return this._text; }
  appendChild(c) { this.children.push(c); return c; }
  addEventListener(t, fn) { (this.handlers[t] = this.handlers[t] || []).push(fn); }
  setAttribute() {}
  fire(t) { (this.handlers[t] || []).forEach(f => f({ target: this, stopPropagation() {} })); }
}
const els = {};
const getEl = id => (els[id] = els[id] || new El('div', id));
const mk = ds => { const e = new El('button'); Object.assign(e.dataset, ds || {}); return e; };
const presetStubs = ['CH4', 'C2H6', 'C3H8', 'N2', 'CO2', 'NG', 'BIO'].map(p => mk({ p }));
const figStubs = ['fig1', 'ch4', 'fig3', 'fig4', 'fig5', 'fig6', 'fig7'].map(f => mk({ fig: f }));
const modeStubs = ['curve', 'px', 'state', 'overview'].map(m => mk({ m }));
const tabStubs = ['calc', 'paper', 'params', 'valid'].map(t => mk({ t }));
const paneStubs = ['tab-calc', 'tab-paper', 'tab-params', 'tab-valid'].map(id => new El('div', id));
const document = {
  getElementById: getEl,
  createElement: t => new El(t),
  createElementNS: (ns, t) => new El(t),
  querySelectorAll(sel) {
    if (sel.includes('.presets')) return presetStubs;
    if (sel.includes('data-fig')) return figStubs;
    if (sel.includes('#modeRow')) return modeStubs;
    if (sel.includes('#tabs')) return tabStubs;
    if (sel.includes('.tabpane')) return paneStubs;
    return [];
  },
  querySelector(sel) {
    let m = sel.match(/data-t="(\w+)"/); if (m) return tabStubs.find(s => s.dataset.t === m[1]);
    m = sel.match(/data-m="(\w+)"/); if (m) return modeStubs.find(s => s.dataset.m === m[1]);
    m = sel.match(/data-fig="(\w+)"/); if (m) return figStubs.find(s => s.dataset.fig === m[1]);
    if (sel === '#cmpRow .btn') return getEl('cmpBtn');
    if (sel.includes('.presets')) return presetStubs[0];
    return new El('div');
  }
};
const window = { addEventListener() {} };
global.document = document; global.window = window;
global.setTimeout = (fn) => { fn(); return 0; };

getEl('c_CH4').value = '100'; ['c_C2H6', 'c_C3H8', 'c_N2', 'c_CO2'].forEach(i => getEl(i).value = '0');
getEl('tmin').value = '263.15'; getEl('tmax').value = '303.15'; getEl('npts').value = '41';
getEl('showanch').checked = true; getEl('logy').checked = false;
getEl('structmode').value = 'auto';
getEl('pxT').value = '277.6'; getEl('pxn').value = '49'; getEl('stateT').value = '277.6';
const pxA = getEl('pxA'), pxB = getEl('pxB');
['CH4', 'C2H6', 'C3H8', 'N2', 'CO2'].forEach(g => { pxA.options.push({ value: g }); pxB.options.push({ value: g }); });
pxA.selectedIndex = 0; pxB.selectedIndex = 2;
for (const sel of [pxA, pxB]) Object.defineProperty(sel, 'value', { get() { return (this.options[this.selectedIndex] || {}).value || ''; } });
getEl('progress').appendChild(new El('div'));

const html = fs.readFileSync(path.join(__dirname, '..', 'index.html'), 'utf8');
const blocks = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m => m[1]);

const nap = () => new Promise(r => setImmediate(r));
let fails = 0;
function chk(name, cond, extra) { console.log((cond ? 'ok   ' : 'FAIL ') + name + (extra ? '  | ' + extra : '')); if (!cond) fails++; }

(async () => {
  (0, eval)(blocks.join('\n;\n'));
  await nap(); await nap();
  console.log('== init (CH4 curve 263.15-303.15, 41 pts) ==');
  chk('curve done', /Done — 41\/41/.test(getEl('status').textContent), getEl('status').textContent.slice(0, 130));
  chk('ice region flagged', getEl('status').textContent.includes('ice') || getEl('status').textContent.includes('H–I–V'));
  chk('legend has structures', getEl('legend').innerHTML.includes('structure I'));
  chk('resultTable rows', (getEl('resultTable').innerHTML.match(/<tr>/g) || []).length >= 40);
  chk('chart svg built', getEl('chart').children.length > 0);
  chk('chart elements', getEl('chart').children[0] && getEl('chart').children[0].children.length > 20);
  chk('occTable lazy (empty before paper tab)', !getEl('occTable').innerHTML);
  chk('cmpRow shown for pure gas', getEl('cmpRow').style.display === 'block');

  console.log('== tabs ==');
  tabStubs[1].fire('click'); await nap();
  chk('paper tab on', tabStubs[1].classList.contains('on') && paneStubs[1].classList.contains('on'));
  const oc = getEl('occTable').innerHTML;
  chk('occTable live rows', oc.includes('273.2') && oc.includes('276.7'), oc.slice(0, 50));
  chk('occTable CH4 273.2 theta', oc.includes('0.97') , oc.match(/0\.9\d\d\d/g)?.slice(0,4)?.join(','));
  tabStubs[0].fire('click'); await nap();

  console.log('== NG preset (Fig 7 gas) ==');
  presetStubs[5].fire('click'); await nap();
  chk('NG comp set', Math.abs(parseFloat(getEl('c_CH4').value) - 93.2) < 1e-9 && Math.abs(parseFloat(getEl('c_N2').value) - 0.43) < 1e-9, `${getEl('c_CH4').value}/${getEl('c_N2').value}`);
  getEl('tmin').value = '277.65'; getEl('tmax').value = '296.65'; getEl('npts').value = '20';
  getEl('calc').fire('click'); await nap(); await nap();
  chk('NG curve done sII', /Done — 20\/20/.test(getEl('status').textContent) && getEl('status').textContent.includes('sII'), getEl('status').textContent.slice(0, 100));

  console.log('== px mode: Fig5 CH4+C3H8 ==');
  modeStubs[1].fire('click'); await nap();
  getEl('pxn').value = '49'; getEl('pxT').value = '277.6';
  getEl('calc').fire('click'); await nap(); await nap();
  const pxst = getEl('status').textContent;
  chk('px transition found', /Structure transitions/.test(pxst), pxst.slice(0, 110));
  chk('px chart drawn', getEl('chart').children.length > 0 && getEl('chart').children[0].children.length > 20);

  console.log('== state mode: CH4 277.6 K ==');
  presetStubs[0].fire('click'); await nap();
  modeStubs[2].fire('click'); await nap();
  getEl('stateT').value = '277.6';
  getEl('calc').fire('click'); await nap(); await nap();
  chk('state done', getEl('status').textContent.startsWith('Done.'), getEl('status').textContent.slice(0, 60));
  const sp = getEl('statePanel').innerHTML;
  chk('statePanel details', sp.includes('Langmuir constants') && sp.includes('Water fugacities') && sp.includes('hydration number'), sp.slice(0, 60));

  console.log('== fig preset fig4 (CO2) + AADP compare ==');
  figStubs[3].fire('click'); await nap(); await nap();
  chk('fig4 comp CO2', parseFloat(getEl('c_CO2').value) === 100, getEl('c_CO2').value);
  modeStubs[0].fire('click'); await nap();
  presetStubs[0].fire('click'); await nap();
  getEl('tmin').value = '273.15'; getEl('tmax').value = '298.15'; getEl('npts').value = '26';
  getEl('calc').fire('click'); await nap(); await nap();
  getEl('cmpBtn').fire('click'); await nap(); await nap();
  const co = getEl('cmpOut').innerHTML;
  chk('cmpOut AADP', co.includes('AADP'), co.replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').slice(0, 120));
  chk('cmpOut dev+parity plots', getEl('chart').children.length >= 3);

  console.log('== overview mode ==');
  modeStubs[3].fire('click'); await nap();
  getEl('npts').value = '11';
  getEl('calc').fire('click'); await nap(); await nap(); await nap();
  chk('overview done', getEl('status').textContent.includes('Overview done'), getEl('status').textContent.slice(0, 60));
  chk('overview chart', getEl('chart').children.length > 0 && getEl('chart').children[0].children.length > 30);

  console.log('== fig1 solubility ==');
  figStubs[0].fire('click'); await nap(); await nap(); await nap();
  chk('solub done', getEl('status').textContent.includes('Solubility diagram done'), getEl('status').textContent.slice(0, 60));
  chk('solub chart', getEl('chart').children.length > 0);

  console.log(fails === 0 ? '=== ALL SMOKE CHECKS PASS ===' : `=== ${fails} FAILURES ===`);
  process.exit(fails ? 1 : 0);
})().catch(e => { console.error('SMOKE FAIL:', e); process.exit(1); });
