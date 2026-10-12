/* cutseq static docs — lightweight client search + back-to-top.
   Search index lives at assets/js/search.json: [{title,url,text}].
   No external dependency; degrades gracefully if the index cannot be read. */
(function () {
  'use strict';

  var input = document.getElementById('site-search');
  var results = document.getElementById('site-search-results');
  var index = null;
  var INDEX_PATH = 'assets/js/search.json';

  function esc(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }

  function loadIndex() {
    fetch(INDEX_PATH)
      .then(function (r) { return r.ok ? r.json() : Promise.reject(new Error('bad status')); })
      .then(function (data) {
        index = Array.isArray(data) ? data : (data.entries || []);
      })
      .catch(function () {
        // Offline / file:// or missing index: disable search quietly.
        index = [];
      });
  }

  function showMessage(msg) {
    if (!results) return;
    results.hidden = false;
    results.innerHTML = '<span class="search-results-empty">' + esc(msg) + '</span>';
  }

  function run(query) {
    if (!results) return;
    query = (query || '').trim().toLowerCase();
    if (!query) {
      results.hidden = true;
      results.innerHTML = '';
      return;
    }
    if (!index) {
      showMessage('Search index still loading…');
      return;
    }
    var terms = query.split(/\s+/).filter(Boolean);
    var found = index.filter(function (entry) {
      var hay = (entry.title + ' ' + (entry.text || '')).toLowerCase();
      return terms.every(function (t) { return hay.indexOf(t) !== -1; });
    }).slice(0, 12);

    if (!found.length) {
      showMessage('No matches.');
      return;
    }
    results.hidden = false;
    results.innerHTML = found.map(function (entry) {
      var title = esc(entry.title || 'Untitled');
      var snippet = esc((entry.text || '').slice(0, 110));
      return '<a class="search-result-item" href="' + esc(entry.url) + '">' +
        '<strong>' + title + '</strong>' +
        (snippet ? '<span class="sr-snippet">' + snippet + '</span>' : '') +
        '</a>';
    }).join('');
  }

  function onClick(e) {
    if (results && !results.contains(e.target)) {
      results.hidden = true;
    }
  }

  if (input) {
    input.addEventListener('input', function () { run(input.value); });
    input.addEventListener('keydown', function (e) { if (e.key === 'Escape') { input.value = ''; run(''); } });
    document.addEventListener('click', onClick);
  }

  /* Back-to-top button */
  var btt = document.getElementById('back-to-top');
  if (btt) {
    btt.addEventListener('click', function (e) {
      e.preventDefault();
      window.scrollTo({ top: 0, behavior: 'smooth' });
    });
  }

  loadIndex();
})();
