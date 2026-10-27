/**
 * charts.js — Interactive Canvas Charts for PCA and Read Budgets
 */

class PcaPlot {
  constructor(canvasId, statsContainerId) {
    this.canvas = document.getElementById(canvasId);
    this.ctx = this.canvas.getContext('2d');
    this.statsContainer = document.getElementById(statsContainerId);
    this.data = null;
    this.currentMode = 'filtered_12'; // or 'all_16'
    this.hoveredPoint = null;
    this.animationProgress = 1;
    this.points = [];

    this.tissueColors = {
      'Nodule': '#f59e0b',
      'Exudophore': '#38bdf8',
      'Fuzzy mycelium': '#c084fc',
      'Exuding mycelium': '#34d399'
    };

    this.initEvents();
  }

  loadData(pcaJson) {
    this.data = pcaJson;
    this.setMode(this.currentMode);
  }

  initEvents() {
    window.addEventListener('resize', () => this.resize());
    this.canvas.addEventListener('mousemove', (e) => this.handleMouseMove(e));
    this.canvas.addEventListener('mouseleave', () => {
      this.hoveredPoint = null;
      this.draw();
    });
  }

  resize() {
    const rect = this.canvas.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = rect.width * dpr;
    this.canvas.height = rect.height * dpr;
    this.ctx.scale(dpr, dpr);
    this.displayWidth = rect.width;
    this.displayHeight = rect.height;
    this.draw();
  }

  setMode(mode) {
    this.currentMode = mode;
    if (!this.data) return;

    const dataset = this.data[mode];
    this.points = dataset.points;
    this.updateStats(dataset);
    this.resize();
  }

  updateStats(dataset) {
    if (!this.statsContainer) return;
    const isFiltered = this.currentMode === 'filtered_12';
    this.statsContainer.innerHTML = `
      <div class="pca-stats-card">
        <div class="control-label">PC1 vs Sequencing Depth Correlation</div>
        <div class="pca-stats-val" style="color: ${isFiltered ? 'var(--accent-emerald)' : 'var(--accent-rose)'};">
          r = ${dataset.cor_depth > 0 ? '+' : ''}${dataset.cor_depth.toFixed(3)}
        </div>
        <div style="font-size: 0.8rem; color: var(--text-secondary); margin-top: 0.25rem;">
          ${isFiltered 
            ? 'Depth artefact successfully eliminated. PC1 captures biological tissue differentiation.' 
            : 'Severe confounder: PC1 predominantly tracks sequencing depth rather than biology.'}
        </div>
      </div>
      <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.75rem; margin-top: 0.75rem;">
        <div class="pca-stats-card">
          <div class="control-label">Libraries Retained</div>
          <div class="pca-stats-val">${this.points.length} / 16</div>
        </div>
        <div class="pca-stats-card">
          <div class="control-label">LRT DEGs (padj &lt; 0.05)</div>
          <div class="pca-stats-val" style="color: var(--accent-cyan);">${isFiltered ? '281' : '76'}</div>
        </div>
      </div>
    `;
  }

  draw() {
    if (!this.points || !this.points.length) return;
    const ctx = this.ctx;
    const w = this.displayWidth;
    const h = this.displayHeight;

    ctx.clearRect(0, 0, w, h);

    // Padding
    const pad = { top: 40, right: 40, bottom: 50, left: 60 };
    const plotW = w - pad.left - pad.right;
    const plotH = h - pad.top - pad.bottom;

    // Bounds calculation
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    this.points.forEach(p => {
      if (p.pc1 < minX) minX = p.pc1;
      if (p.pc1 > maxX) maxX = p.pc1;
      if (p.pc2 < minY) minY = p.pc2;
      if (p.pc2 > maxY) maxY = p.pc2;
    });

    // Add margin
    const marginX = (maxX - minX) * 0.15 || 5;
    const marginY = (maxY - minY) * 0.15 || 5;
    minX -= marginX;
    maxX += marginX;
    minY -= marginY;
    maxY += marginY;

    const toScreenX = (val) => pad.left + ((val - minX) / (maxX - minX)) * plotW;
    const toScreenY = (val) => pad.top + ((maxY - val) / (maxY - minY)) * plotH;

    // Grid lines & Zero axes
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.07)';
    const isLight = document.documentElement.getAttribute('data-theme') === 'light';
    if (isLight) ctx.strokeStyle = 'rgba(0, 0, 0, 0.08)';

    ctx.lineWidth = 1;

    // Draw zero axes if in range
    if (minX <= 0 && maxX >= 0) {
      const zx = toScreenX(0);
      ctx.beginPath();
      ctx.moveTo(zx, pad.top);
      ctx.lineTo(zx, h - pad.bottom);
      ctx.strokeStyle = isLight ? 'rgba(0,0,0,0.18)' : 'rgba(255,255,255,0.18)';
      ctx.stroke();
    }
    if (minY <= 0 && maxY >= 0) {
      const zy = toScreenY(0);
      ctx.beginPath();
      ctx.moveTo(pad.left, zy);
      ctx.lineTo(w - pad.right, zy);
      ctx.strokeStyle = isLight ? 'rgba(0,0,0,0.18)' : 'rgba(255,255,255,0.18)';
      ctx.stroke();
    }

    // Axis Labels
    ctx.fillStyle = isLight ? '#475569' : '#94a3b8';
    ctx.font = '600 12px ' + getComputedStyle(document.body).fontFamily;
    ctx.textAlign = 'center';
    ctx.fillText('PC1 (Expression Variance)', pad.left + plotW / 2, h - 15);

    ctx.save();
    ctx.translate(18, pad.top + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.fillText('PC2 (Expression Variance)', 0, 0);
    ctx.restore();

    // Draw Points
    this.points.forEach(p => {
      const sx = toScreenX(p.pc1);
      const sy = toScreenY(p.pc2);
      p._screenX = sx;
      p._screenY = sy;

      const isHovered = this.hoveredPoint === p;
      const radius = isHovered ? 9 : 7;
      const color = this.tissueColors[p.tissue] || '#94a3b8';

      ctx.beginPath();
      ctx.arc(sx, sy, radius, 0, Math.PI * 2);
      ctx.fillStyle = color;
      ctx.fill();

      if (p.flagged_low) {
        // Mark flagged low samples with warning ring
        ctx.strokeStyle = '#f43f5e';
        ctx.lineWidth = 2;
        ctx.stroke();

        ctx.fillStyle = '#f43f5e';
        ctx.font = '700 10px monospace';
        ctx.textAlign = 'center';
        ctx.fillText('LOW', sx, sy - 12);
      } else {
        ctx.strokeStyle = isLight ? '#ffffff' : '#111827';
        ctx.lineWidth = 2;
        ctx.stroke();
      }

      // Draw Well label
      ctx.fillStyle = isLight ? '#0f172a' : '#f8fafc';
      ctx.font = '600 10px ' + getComputedStyle(document.body).fontFamily;
      ctx.textAlign = 'center';
      ctx.fillText(p.well, sx, sy + (isHovered ? 18 : 16));
    });

    // Tooltip rendering
    if (this.hoveredPoint) {
      const p = this.hoveredPoint;
      const text1 = `${p.sample} (${p.well})`;
      const text2 = `${p.tissue} • ${p.mrna.toLocaleString()} mRNA counts`;
      const text3 = `PC1: ${p.pc1.toFixed(1)}, PC2: ${p.pc2.toFixed(1)}`;
      const text4 = p.flagged_low ? '⚠️ Excluded (<15k threshold)' : '✓ Retained in analysis';

      ctx.font = '600 12px sans-serif';
      const m1 = ctx.measureText(text1).width;
      const m2 = ctx.measureText(text2).width;
      const m3 = ctx.measureText(text3).width;
      const maxW = Math.max(m1, m2, m3) + 24;
      const tooltipH = p.flagged_low ? 82 : 68;

      let tx = p._screenX + 14;
      let ty = p._screenY - tooltipH / 2;
      if (tx + maxW > w - 10) tx = p._screenX - maxW - 14;
      if (ty < 10) ty = 10;
      if (ty + tooltipH > h - 10) ty = h - tooltipH - 10;

      ctx.fillStyle = isLight ? 'rgba(15, 23, 42, 0.95)' : 'rgba(17, 24, 39, 0.95)';
      ctx.beginPath();
      ctx.roundRect(tx, ty, maxW, tooltipH, 6);
      ctx.fill();
      ctx.strokeStyle = 'rgba(56, 189, 248, 0.4)';
      ctx.lineWidth = 1;
      ctx.stroke();

      ctx.textAlign = 'left';
      ctx.fillStyle = '#f8fafc';
      ctx.fillText(text1, tx + 12, ty + 18);
      ctx.fillStyle = '#94a3b8';
      ctx.font = '11px sans-serif';
      ctx.fillText(text2, tx + 12, ty + 36);
      ctx.fillText(text3, tx + 12, ty + 52);
      if (p.flagged_low) {
        ctx.fillStyle = '#f43f5e';
        ctx.fillText(text4, tx + 12, ty + 68);
      }
    }
  }

  handleMouseMove(e) {
    const rect = this.canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;

    let nearest = null;
    let minDist = 18;

    this.points.forEach(p => {
      const dist = Math.hypot(p._screenX - mx, p._screenY - my);
      if (dist < minDist) {
        minDist = dist;
        nearest = p;
      }
    });

    if (this.hoveredPoint !== nearest) {
      this.hoveredPoint = nearest;
      this.canvas.style.cursor = nearest ? 'pointer' : 'default';
      this.draw();
    }
  }
}

window.PcaPlot = PcaPlot;
