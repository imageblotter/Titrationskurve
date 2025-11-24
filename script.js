// ==============================
// Numerische Konstanten
// ==============================
const KW = 1e-14;

// Utilities
function log10(x) { return Math.log(x) / Math.LN10; }
function clamp(x, lo, hi) { return Math.min(Math.max(x, lo), hi); }

// ==========================================
// WURZELFINDER IM LOG[H+]
// ==========================================
function bisectInLogH(F, yGuess = null) {
    let yLo = -14, yHi = 0;
    let FLo = F(yLo), FHi = F(yHi);

    if (FLo === 0) return yLo;
    if (FHi === 0) return yHi;

    let found = (FLo * FHi < 0);
    if (!found) {
        const scanOrder = [];
        const N = 56; // Schritt ~0,25 pH
        if (yGuess !== null && isFinite(yGuess)) {
            const step = (0 - (-14)) / N;
            const ys = Array.from({ length: N + 1 }, (_, i) => -14 + i * step);
            ys.sort((a, b) => Math.abs(a - yGuess) - Math.abs(b - yGuess));
            for (let k = 0; k < ys.length - 1; k++) {
                const a = Math.min(ys[k], ys[k + 1]);
                const b = Math.max(ys[k], ys[k + 1]);
                scanOrder.push([a, b]);
            }
        } else {
            const step = (0 - (-14)) / N;
            for (let i = 0; i < N; i++) {
                scanOrder.push([-14 + i * step, -14 + (i + 1) * step]);
            }
        }
        for (const [a, b] of scanOrder) {
            const Fa = F(a);
            const Fb = F(b);
            if (Fa === 0) return a;
            if (Fb === 0) return b;
            if (Fa * Fb < 0) {
                yLo = a; FLo = Fa;
                yHi = b; FHi = Fb;
                found = true;
                break;
            }
        }
    }
    if (!found) {
        const yMid = (yLo + yHi) / 2;
        const FyMid = F(yMid);
        const candidates = [
            { y: yLo, v: Math.abs(FLo) },
            { y: yHi, v: Math.abs(FHi) },
            { y: yMid, v: Math.abs(FyMid) },
        ];
        candidates.sort((a, b) => a.v - b.v);
        return candidates[0].y;
    }

    let it = 0;
    while ((yHi - yLo) > 1e-8 && it < 100) {
        const ym = 0.5 * (yLo + yHi);
        const Fm = F(ym);
        if (Fm === 0) return ym;
        if (FLo * Fm < 0) {
            yHi = ym; FHi = Fm;
        } else {
            yLo = ym; FLo = Fm;
        }
        it++;
    }
    return 0.5 * (yLo + yHi);
}

// H-Klammerung in [1e-14, 1]
function evalWithH(fn, y) {
    const H = Math.pow(10, clamp(y, -14, 0));
    return fn(H);
}

// -----------------------------------------------------------
// Schwache SÄURE (HA) gegen starke BASE – Ladungsbilanz
// -----------------------------------------------------------
function pH_weakAcid_root(Ca, VaL, Ct, pKa, VaddL, yWarm = null) {
    const Va = VaL / 1000;
    const Vadd = VaddL / 1000;
    const Vtot = Va + Vadd;
    const Ka = Math.pow(10, -pKa);
    const n_HA0 = Ca * Va;
    const CT = n_HA0 / Vtot;
    const Cna = (Ct * Vadd) / Vtot;

    // f(H) = H + Cna - (KW/H) - (CT*Ka/(Ka + H)) = 0
    const F = (y) => evalWithH((H) => {
        const OH = KW / H;
        const A = CT * Ka / (Ka + H);
        return H + Cna - OH - A;
    }, y);

    const y = bisectInLogH(F, yWarm);
    return -y;
}

// -----------------------------------------------------------
// Schwache BASE (B) gegen starke SÄURE – Ladungsbilanz
// -----------------------------------------------------------
function pH_weakBase_root(Cb, VaL, Ct, pKb, VaddL, yWarm = null) {
    const Va = VaL / 1000;
    const Vadd = VaddL / 1000;
    const Vtot = Va + Vadd;
    const Kb = Math.pow(10, -pKb);
    const KaPrime = KW / Kb;
    const n_B0 = Cb * Va;
    const CT = n_B0 / Vtot;
    const Ccl = (Ct * Vadd) / Vtot;

    // g(H) = H + (CT*H/(H + KaPrime)) - (KW/H) - Ccl = 0
    const G = (y) => evalWithH((H) => {
        const OH = KW / H;
        const BH = CT * H / (H + KaPrime);
        return H + BH - OH - Ccl;
    }, y);

    const y = bisectInLogH(G, yWarm);
    return -y;
}

// ===========================================================
// STARK/STARK-Logik
// ===========================================================
// ===========================================================
// STARK/STARK-Logik (Exakt mit quadratischer Gleichung)
// ===========================================================
function pH_strongAcid_vs_strongBase(Ca, VaL, Ct, VaddL) {
    const Va = VaL / 1000;
    const Vb = VaddL / 1000;
    const Vtot = Va + Vb;

    // Konzentrationen im Gesamtvolumen
    const C_acid = (Ca * Va) / Vtot;
    const C_base = (Ct * Vb) / Vtot;

    // Ladungsbilanz: [H+] + [Na+] = [OH-] + [Cl-]
    // H + C_base = KW/H + C_acid
    // H^2 + (C_base - C_acid)*H - KW = 0

    const delta = C_base - C_acid;

    // Fall 1: Base im Überschuss (delta > 0) -> Wir lösen nach OH-
    // [OH-] + [Cl-] = [H+] + [Na+]
    // OH + C_acid = KW/OH + C_base
    // OH^2 + (C_acid - C_base)*OH - KW = 0 -> OH^2 - delta*OH - KW = 0
    if (delta > 0) {
        const OH = (delta + Math.sqrt(delta * delta + 4 * KW)) / 2;
        return 14 - (-log10(OH));
    }
    // Fall 2: Säure im Überschuss (delta < 0) -> Wir lösen nach H+
    // H^2 + delta*H - KW = 0. Sei C_exc = -delta > 0.
    // H^2 - C_exc*H - KW = 0
    else if (delta < 0) {
        const C_exc = -delta;
        const H = (C_exc + Math.sqrt(C_exc * C_exc + 4 * KW)) / 2;
        return -log10(H);
    }
    // Fall 3: Äquivalenz
    else {
        return 7.0;
    }
}

function pH_strongBase_vs_strongAcid(Cb, VaL, Ct, VaddL) {
    const Va = VaL / 1000;
    const VaH = VaddL / 1000;
    const Vtot = Va + VaH;

    const C_base = (Cb * Va) / Vtot;
    const C_acid = (Ct * VaH) / Vtot;

    const delta = C_base - C_acid;

    if (delta > 0) {
        const OH = (delta + Math.sqrt(delta * delta + 4 * KW)) / 2;
        return 14 - (-log10(OH));
    } else if (delta < 0) {
        const C_exc = -delta;
        const H = (C_exc + Math.sqrt(C_exc * C_exc + 4 * KW)) / 2;
        return -log10(H);
    } else {
        return 7.0;
    }
}

// ======================================
// Plugin: Gestrichelte Vertikallinie bei Äquivalenz
// ======================================
let __eqIndex = null; // wird in simulate() gesetzt
const eqLinePlugin = {
    id: 'eqLine',
    afterDatasetsDraw(chart) {
        if (__eqIndex == null) return;
        const yScale = chart.scales.y;
        const meta = chart.getDatasetMeta(0);
        if (!meta || !meta.data || !meta.data[__eqIndex]) return;

        const x = meta.data[__eqIndex].x; // exakte Pixelposition des Datenpunkts bei Äquivalenz
        const ctx = chart.ctx;
        ctx.save();
        ctx.setLineDash([6, 6]);
        ctx.lineWidth = 1.5;
        ctx.strokeStyle = '#64748b'; // slate-500
        ctx.beginPath();
        ctx.moveTo(x, yScale.getPixelForValue(yScale.max));
        ctx.lineTo(x, yScale.getPixelForValue(yScale.min));
        ctx.stroke();
        ctx.restore();
    }
};
Chart.register(eqLinePlugin);

// ======================================
// Diagramm-Verdrahtung
// ======================================
const ctx = document.getElementById('titrationChart').getContext('2d');
const chart = new Chart(ctx, {
    type: 'line',
    data: {
        labels: [],
        datasets: [{
            label: 'Titrationskurve',
            data: [],
            fill: false,
            borderColor: '#3b82f6', // Primary Blue
            backgroundColor: '#3b82f6',
            borderWidth: 2,
            pointRadius: 0,
            tension: 0.12
        }]
    },
    options: {
        animation: false,
        responsive: true,
        maintainAspectRatio: false,
        scales: {
            x: {
                title: { display: true, text: 'Zugegebener Titrant (Äquivalente oder mL)' },
                grid: { color: '#e2e8f0' }
            },
            y: {
                min: 0,
                max: 14,
                title: { display: true, text: 'pH' },
                grid: { color: '#e2e8f0' }
            }
        },
        plugins: {
            legend: { display: false },
            tooltip: {
                callbacks: {
                    label: (ctx) => `pH: ${ctx.parsed.y.toFixed(4)}`
                }
            }
        }
    }
});

// Inputs
const analyteTypeEl = document.getElementById('analyteType');
const analyteStrengthEl = document.getElementById('analyteStrength');
const CaCbEl = document.getElementById('CaCb');
const VaEl = document.getElementById('Va');
const CtEl = document.getElementById('Ct');
const pKaKbEl = document.getElementById('pKaKb');
const xUnitsEl = document.getElementById('xUnits');
const initialPHEl = document.getElementById('initialPH');
const eqValueEl = document.getElementById('eqValue');
const pointCountEl = document.getElementById('pointCount');

function updateEquivalenceBadge() {
    const CaCb = parseFloat(CaCbEl.value);
    const Va = parseFloat(VaEl.value);
    const Ct = parseFloat(CtEl.value);
    const Veq_mL = (CaCb * Va) / Math.max(Ct, 1e-12);
    eqValueEl.textContent = Veq_mL.toFixed(2);
}

analyteTypeEl.addEventListener('change', updateEquivalenceBadge);
CaCbEl.addEventListener('input', updateEquivalenceBadge);
VaEl.addEventListener('input', updateEquivalenceBadge);
CtEl.addEventListener('input', updateEquivalenceBadge);
updateEquivalenceBadge();

// CSV-Export
document.getElementById('exportBtn').addEventListener('click', () => {
    if (!chart.data.labels.length) return;
    const xLabel = xUnitsEl.checked ? 'V_Zugabe (mL)' : 'Äquivalente';
    let csv = `${xLabel},pH\n`;
    for (let i = 0; i < chart.data.labels.length; i++) {
        csv += `${chart.data.labels[i]},${chart.data.datasets[0].data[i]}\n`;
    }
    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'titrationskurve.csv';
    a.click();
    URL.revokeObjectURL(url);
});

// Zurücksetzen
document.getElementById('resetBtn').addEventListener('click', () => {
    analyteTypeEl.value = 'acid';
    analyteStrengthEl.value = 'weak';
    CaCbEl.value = '0.10';
    VaEl.value = '25.00';
    CtEl.value = '0.10';
    pKaKbEl.value = '4.76';
    xUnitsEl.checked = false;
    updateEquivalenceBadge();
    chart.data.labels = [];
    chart.data.datasets[0].data = [];
    chart.update();
    initialPHEl.textContent = '—';
    __eqIndex = null;
    chart.update();
});

// Fullscreen Toggle
document.getElementById('fullscreenBtn').addEventListener('click', () => {
    if (!document.fullscreenElement) {
        document.documentElement.requestFullscreen().catch((e) => {
            console.error(`Error attempting to enable fullscreen mode: ${e.message} (${e.name})`);
        });
    } else {
        if (document.exitFullscreen) {
            document.exitFullscreen();
        }
    }
});

// ======================================
// Simulation
// ======================================
function simulate() {
    const analyteType = analyteTypeEl.value;          // 'acid' oder 'base'
    const analyteStrength = analyteStrengthEl.value;  // 'weak' oder 'strong'
    const CaCb = Math.max(parseFloat(CaCbEl.value), 0);
    const Va = Math.max(parseFloat(VaEl.value), 0);
    const Ct = Math.max(parseFloat(CtEl.value), 0);
    const pKaKb = parseFloat(pKaKbEl.value);
    const usemL = xUnitsEl.checked;

    // x-Gitter: 0 bis 2 Äquivalente (inkl.) oder mL-Bereich
    const npts = 401;
    const labels = new Array(npts);
    const data = new Array(npts);
    pointCountEl.textContent = npts.toString();

    const Veq_mL = (CaCb * Va) / Math.max(Ct, 1e-12);
    const maxV_mL = 2 * Veq_mL;

    // Anfangs-pH anzeigen
    let initialPH = NaN;
    if (analyteStrength === 'weak') {
        if (analyteType === 'acid') {
            initialPH = pH_weakAcid_root(CaCb, Va, Ct, pKaKb, 0.0, -Math.max(0, 7 - pKaKb));
        } else {
            initialPH = pH_weakBase_root(CaCb, Va, Ct, pKaKb, 0.0, -Math.max(0, pKaKb - 7));
        }
    } else {
        if (analyteType === 'acid') {
            initialPH = pH_strongAcid_vs_strongBase(CaCb, Va, Ct, 0.0);
        } else {
            initialPH = pH_strongBase_vs_strongAcid(CaCb, Va, Ct, 0.0);
        }
    }
    initialPHEl.textContent = isFinite(initialPH) ? initialPH.toFixed(4) : '—';

    // Iteration über das Gitter
    let lastY = null; // Warmstart
    for (let i = 0; i < npts; i++) {
        const frac = i / (npts - 1);
        const Vadd_mL = frac * maxV_mL;
        const eq = Vadd_mL / Math.max(Veq_mL, 1e-12);

        labels[i] = usemL ? Vadd_mL.toFixed(2) : eq.toFixed(4);

        let pH = 7.0;

        if (analyteStrength === 'weak') {
            if (analyteType === 'acid') {
                const pHval = pH_weakAcid_root(CaCb, Va, Ct, pKaKb, Vadd_mL, lastY);
                pH = clamp(pHval, 0, 14);
                lastY = -pH;
            } else {
                const pHval = pH_weakBase_root(CaCb, Va, Ct, pKaKb, Vadd_mL, lastY);
                pH = clamp(pHval, 0, 14);
                lastY = -pH;
            }
        } else {
            if (analyteType === 'acid') {
                pH = pH_strongAcid_vs_strongBase(CaCb, Va, Ct, Vadd_mL);
            } else {
                pH = pH_strongBase_vs_strongAcid(CaCb, Va, Ct, Vadd_mL);
            }
            pH = clamp(pH, 0, 14);
            lastY = -pH;
        }

        data[i] = pH;
    }

    // Diagramm aktualisieren
    chart.data.labels = labels;
    chart.data.datasets[0].data = data;

    if (xUnitsEl.checked) {
        __eqIndex = Math.round((Veq_mL / (2 * Veq_mL)) * (npts - 1));
    } else {
        __eqIndex = Math.round((1 / 2) * (npts - 1));
    }

    chart.update();
}

document.getElementById('simulateBtn').addEventListener('click', simulate);

// Erste Darstellung
simulate();
