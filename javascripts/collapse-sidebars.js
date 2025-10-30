document.addEventListener('DOMContentLoaded', function() {
  // Create collapse buttons
  const navCollapseBtn = document.createElement('button');
  navCollapseBtn.className = 'md-nav-collapse-btn';
  navCollapseBtn.innerHTML = '◀';
  navCollapseBtn.title = 'Collapse navigation';
  
  const tocCollapseBtn = document.createElement('button');
  tocCollapseBtn.className = 'md-toc-collapse-btn';
  tocCollapseBtn.innerHTML = '▶';
  tocCollapseBtn.title = 'Collapse table of contents';
  
  // Insert buttons into body (not sidebars)
  document.body.appendChild(navCollapseBtn);
  document.body.appendChild(tocCollapseBtn);
  
  // Toggle functionality
  navCollapseBtn.addEventListener('click', function() {
    const nav = document.querySelector('.md-sidebar--primary');
    const content = document.querySelector('.md-content');
    const main = document.querySelector('.md-main__inner');
    
    if (nav) {
      nav.classList.toggle('collapsed');
      navCollapseBtn.classList.toggle('collapsed');
      navCollapseBtn.innerHTML = nav.classList.contains('collapsed') ? '▶' : '◀';
      
      // Expand content when collapsed
      if (content) content.classList.toggle('nav-collapsed');
      if (main) main.classList.toggle('nav-collapsed');
    }
  });
  
  tocCollapseBtn.addEventListener('click', function() {
    const toc = document.querySelector('.md-sidebar--secondary');
    const content = document.querySelector('.md-content');
    const main = document.querySelector('.md-main__inner');
    
    if (toc) {
      toc.classList.toggle('collapsed');
      tocCollapseBtn.classList.toggle('collapsed');
      tocCollapseBtn.innerHTML = toc.classList.contains('collapsed') ? '◀' : '▶';
      
      // Expand content when collapsed
      if (content) content.classList.toggle('toc-collapsed');
      if (main) main.classList.toggle('toc-collapsed');
    }
  });
});
