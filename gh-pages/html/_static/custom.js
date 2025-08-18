/* Custom JavaScript for Sanger DNA Damage Analysis Pipeline Documentation */

// Add custom JavaScript here
document.addEventListener('DOMContentLoaded', function() {
    // Example: Add smooth scrolling to anchor links
    const links = document.querySelectorAll('a[href^="#"]');
    links.forEach(link => {
        link.addEventListener('click', function(e) {
            e.preventDefault();
            const target = document.querySelector(this.getAttribute('href'));
            if (target) {
                target.scrollIntoView({
                    behavior: 'smooth'
                });
            }
        });
    });
});
