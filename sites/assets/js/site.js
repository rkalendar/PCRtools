(function(){
  // Year
  var y = document.getElementById("year");
  if (y) y.textContent = new Date().getFullYear();

  // Tool search (index only)
  var search = document.getElementById("toolSearch");
  if (search){
    var cards = Array.prototype.slice.call(document.querySelectorAll(".tool-card"));
    var fn = function(){
      var q = (search.value || "").trim().toLowerCase();
      cards.forEach(function(c){
        var hay = c.getAttribute("data-title") || "";
        c.style.display = (!q || hay.indexOf(q) !== -1) ? "" : "none";
      });
    };
    search.addEventListener("input", fn);
  }
})();
