/**
 * main.js — Interactive Portal Controller
 * Manages theme, data loading, tables, figure gallery, and lightbox
 */

document.addEventListener('DOMContentLoaded', async () => {
  initTheme();
  initNav();
  await loadDataAndInit();
});

// Theme Management
function initTheme() {
  const themeToggle = document.getElementById('theme-toggle');
  const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
  const savedTheme = localStorage.getItem('theme') || (prefersDark ? 'dark' : 'light');

  document.documentElement.setAttribute('data-theme', savedTheme);
  updateThemeIcon(savedTheme);

  if (themeToggle) {
    themeToggle.addEventListener('click', () => {
      const current = document.documentElement.getAttribute('data-theme');
      const next = current === 'dark' ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', next);
      localStorage.setItem('theme', next);
      updateThemeIcon(next);
      if (window.pcaInstance) window.pcaInstance.draw();
    });
  }
}

function updateThemeIcon(theme) {
  const icon = document.getElementById('theme-icon');
  if (icon) {
    icon.textContent = theme === 'dark' ? '☀️' : '🌙';
  }
}

function initNav() {
  document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function(e) {
      e.preventDefault();
      const target = document.querySelector(this.getAttribute('href'));
      if (target) {
        target.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }
    });
  });
}

// Global App State
let appData = {
  budget: [],
  pca: null,
  markers: [],
  modules: [],
  figures: []
};

async function loadDataAndInit() {
  try {
    const [budgetRes, pcaRes, markersRes, modulesRes, figuresRes] = await Promise.all([
      fetch('data/read_budget.json'),
      fetch('data/pca_data.json'),
      fetch('data/robust_markers.json'),
      fetch('data/modules.json'),
      fetch('data/figures_meta.json')
    ]);

    appData.budget = await budgetRes.json();
    appData.pca = await pcaRes.json();
    appData.markers = await markersRes.json();
    appData.modules = await modulesRes.json();
    appData.figures = await figuresRes.json();

    initReadBudgetTable(appData.budget);
    initPcaExplorer(appData.pca);
    initMarkersExplorer(appData.markers);
    initFigureGallery(appData.figures);
  } catch (err) {
    console.error('Error loading data:', err);
  }
}

// Read Budget Table
function initReadBudgetTable(budgetList) {
  const tbody = document.getElementById('budget-tbody');
  if (!tbody) return;

  tbody.innerHTML = budgetList.map(item => {
    const tissueClass = getTissueBadgeClass(item.tissue);
    const statusBadge = item.flagged_low
      ? `<span class="badge badge-danger">Excluded (&lt;15k counts)</span>`
      : `<span class="badge badge-success">Retained (n=${item.tissue === 'Fuzzy mycelium' ? 4 : (item.tissue === 'Exudophore' ? 2 : 3)})</span>`;

    return `
      <tr>
        <td><strong>${item.well}</strong></td>
        <td><code>${item.sample}</code></td>
        <td><span class="badge ${tissueClass}">${item.tissue}</span></td>
        <td>${item.raw_reads.toLocaleString()}</td>
        <td>${item.rrna_reads.toLocaleString()}</td>
        <td><strong>${item.mrna_reads.toLocaleString()}</strong></td>
        <td>${item.mrna_pct.toFixed(2)}%</td>
        <td>${item.genes_det10.toLocaleString()}</td>
        <td>${statusBadge}</td>
      </tr>
    `;
  }).join('');
}

function getTissueBadgeClass(tissue) {
  if (tissue.includes('Nodule')) return 'badge-nodule';
  if (tissue.includes('Exudophore')) return 'badge-exudophore';
  if (tissue.includes('Fuzzy')) return 'badge-fuzzy';
  return 'badge-exuding';
}

// PCA Explorer
function initPcaExplorer(pcaData) {
  const pca = new window.PcaPlot('pca-canvas', 'pca-stats-container');
  pca.loadData(pcaData);
  window.pcaInstance = pca;

  const btnFiltered = document.getElementById('btn-pca-filtered');
  const btnAll = document.getElementById('btn-pca-all');

  if (btnFiltered && btnAll) {
    btnFiltered.addEventListener('click', () => {
      btnFiltered.classList.add('active');
      btnAll.classList.remove('active');
      pca.setMode('filtered_12');
    });

    btnAll.addEventListener('click', () => {
      btnAll.classList.add('active');
      btnFiltered.classList.remove('active');
      pca.setMode('all_16');
    });
  }
}

// Markers Explorer
function initMarkersExplorer(markersList) {
  const tbody = document.getElementById('markers-tbody');
  const searchInput = document.getElementById('marker-search');
  const tissueFilter = document.getElementById('marker-tissue-filter');
  const countIndicator = document.getElementById('marker-count');

  if (!tbody) return;

  function render(filtered) {
    if (countIndicator) {
      countIndicator.textContent = `Showing ${filtered.length} of ${markersList.length} robust markers`;
    }

    if (!filtered.length) {
      tbody.innerHTML = `<tr><td colspan="4" style="text-align: center; color: var(--text-muted); padding: 2rem;">No matching markers found.</td></tr>`;
      return;
    }

    tbody.innerHTML = filtered.map(m => {
      const tissueClass = getTissueBadgeClass(m.tissue);
      return `
        <tr>
          <td><span class="badge ${tissueClass}">${m.tissue}</span></td>
          <td><code>${m.gene}</code></td>
          <td><strong>${m.protein_name}</strong></td>
          <td><span style="font-family: var(--font-mono); font-weight: 700; color: var(--accent-cyan);">${m.tau.toFixed(3)}</span></td>
        </tr>
      `;
    }).join('');
  }

  function applyFilter() {
    const query = (searchInput ? searchInput.value : '').toLowerCase().trim();
    const tissue = tissueFilter ? tissueFilter.value : 'all';

    const filtered = markersList.filter(m => {
      const matchTissue = (tissue === 'all') || (m.tissue === tissue);
      const matchQuery = !query || 
        m.gene.toLowerCase().includes(query) || 
        m.protein_name.toLowerCase().includes(query) || 
        m.tissue.toLowerCase().includes(query);
      return matchTissue && matchQuery;
    });

    render(filtered);
  }

  if (searchInput) searchInput.addEventListener('input', applyFilter);
  if (tissueFilter) tissueFilter.addEventListener('change', applyFilter);

  render(markersList);
}

// Publication Figure Gallery & Lightbox
let currentFigureIndex = 0;
let filteredFigures = [];

function initFigureGallery(figuresList) {
  const galleryGrid = document.getElementById('gallery-grid');
  const filterPills = document.querySelectorAll('.filter-pill');
  filteredFigures = figuresList;

  function renderGallery(figs) {
    if (!galleryGrid) return;
    galleryGrid.innerHTML = figs.map((fig, idx) => `
      <div class="fig-card" data-index="${idx}">
        <div class="fig-thumb-wrap">
          <img src="figures/${fig.id}.png" alt="${fig.number}: ${fig.title}" loading="lazy">
        </div>
        <div class="fig-info">
          <div class="fig-meta-badge">${fig.number} • ${fig.paper}</div>
          <div class="fig-card-title">${fig.title}</div>
          <div class="fig-card-desc">${fig.caption}</div>
        </div>
      </div>
    `).join('');

    // Attach click listeners
    document.querySelectorAll('.fig-card').forEach(card => {
      card.addEventListener('click', () => {
        const idx = parseInt(card.getAttribute('data-index'), 10);
        openLightbox(idx);
      });
    });
  }

  filterPills.forEach(pill => {
    pill.addEventListener('click', () => {
      filterPills.forEach(p => p.classList.remove('active'));
      pill.classList.add('active');
      const cat = pill.getAttribute('data-category');
      filteredFigures = (cat === 'all') 
        ? figuresList 
        : figuresList.filter(f => f.category === cat);
      renderGallery(filteredFigures);
    });
  });

  renderGallery(figuresList);
  initLightbox();
}

function initLightbox() {
  const backdrop = document.getElementById('lightbox-backdrop');
  const closeBtn = document.getElementById('lightbox-close');
  const prevBtn = document.getElementById('lightbox-prev');
  const nextBtn = document.getElementById('lightbox-next');

  if (!backdrop) return;

  closeBtn.addEventListener('click', closeLightbox);
  backdrop.addEventListener('click', (e) => {
    if (e.target === backdrop) closeLightbox();
  });

  prevBtn.addEventListener('click', () => navigateLightbox(-1));
  nextBtn.addEventListener('click', () => navigateLightbox(1));

  document.addEventListener('keydown', (e) => {
    if (!backdrop.classList.contains('open')) return;
    if (e.key === 'Escape') closeLightbox();
    if (e.key === 'ArrowLeft') navigateLightbox(-1);
    if (e.key === 'ArrowRight') navigateLightbox(1);
  });
}

function openLightbox(index) {
  currentFigureIndex = index;
  updateLightboxContent();
  const backdrop = document.getElementById('lightbox-backdrop');
  if (backdrop) backdrop.classList.add('open');
  document.body.style.overflow = 'hidden';
}

function closeLightbox() {
  const backdrop = document.getElementById('lightbox-backdrop');
  if (backdrop) backdrop.classList.remove('open');
  document.body.style.overflow = '';
}

function navigateLightbox(dir) {
  currentFigureIndex = (currentFigureIndex + dir + filteredFigures.length) % filteredFigures.length;
  updateLightboxContent();
}

function updateLightboxContent() {
  const fig = filteredFigures[currentFigureIndex];
  if (!fig) return;

  const titleEl = document.getElementById('lightbox-title');
  const imgEl = document.getElementById('lightbox-img');
  const captionEl = document.getElementById('lightbox-caption');
  const downloadLink = document.getElementById('lightbox-download');

  if (titleEl) titleEl.textContent = `${fig.number}: ${fig.title}`;
  if (imgEl) {
    imgEl.src = `figures/${fig.id}.png`;
    imgEl.alt = `${fig.number}: ${fig.title}`;
  }
  if (captionEl) {
    captionEl.innerHTML = `<strong>${fig.paper}</strong> — ${fig.caption}`;
  }
  if (downloadLink) {
    downloadLink.href = `figures/${fig.id}.png`;
    downloadLink.download = `${fig.id}.png`;
  }
}
