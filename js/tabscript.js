// For sub-tabs
const subTabButtons = document.querySelectorAll('.sub-tab-btn');
subTabButtons.forEach(button => {
  button.addEventListener('click', () => {
    // 1. Identify the parent container (section) of these sub-tabs
    const parentSection = button.closest('.section');

    // 2. Within that parent, find all sub-tab buttons
    const siblingSubTabButtons = parentSection.querySelectorAll('.sub-tab-btn');
    // Remove 'active' from all sibling sub-tab buttons
    siblingSubTabButtons.forEach(btn => btn.classList.remove('active'));

    // 3. Add 'active' to the clicked button
    button.classList.add('active');

    // 4. Hide all sub-tab content inside this parent section
    const subTabContents = parentSection.querySelectorAll('.sub-tab-content');
    subTabContents.forEach(content => content.style.display = 'none');

    // 5. Show the corresponding sub-tab content
    const subTabId = button.getAttribute('data-subtab');
    const subTabContentToShow = document.getElementById(subTabId);
    if (subTabContentToShow) {
      subTabContentToShow.style.display = 'block';
    }
  });
});
 
