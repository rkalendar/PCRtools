var ticks,thisTick,nextTick,state,current=0,delay=5000;speed=20;value=1;

function rotateNews(){
  var ticker = document.getElementById("front_news_text");
  ticks = ticker.getElementsByTagName("li");
  for (var i=0; i < ticks.length; i++) {
    if (i==0) {ticks[i].style.opacity = 1;} 
    else {
      ticks[i].style.opacity = 0;
      ticks[i].style.filter = "alpha(opacity=0)";
      ticks[i].style.display = 'none';
    }
  }
  loopNews();
}

function loopNews(){
  if(current>=ticks.length) current = 0;
  thisTick = (nextTick)?nextTick:ticks[current];
  nextTick = (ticks[++current])?ticks[current]:ticks[current = 0];
  thisTick.style.display = "block";
  state = 'on';
  value = 1;
  window.setTimeout(fader,delay);
}

function fader(){
  if (state=='on') {
    value = value - 0.1;
    value = parseInt(value*100)/100;
    if (value >= 0) {
      thisTick.style.opacity = value;
      thisTick.style.filter = "alpha(opacity="+(value*100)+")";
      window.setTimeout(fader,speed);
    } else {
      thisTick.style.display = "none";
      nextTick.style.display = "block";
      state='off';
      window.setTimeout(fader,speed);
    }
  } else if (state=='off') {
    value = value + 0.1;
    value = parseInt(value*100)/100;
    if(value <= 1){
      nextTick.style.opacity = value;
      nextTick.style.filter = "alpha(opacity="+(value*100)+")";
      window.setTimeout(fader,speed);
    } else {loopNews();}
  }
}

// event listener
if(window.addEventListener) window.addEventListener('load', rotateNews, false); // gecko, safari, konqueror and standard
else if(document.addEventListener) document.addEventListener('load', rotateNews, false); // opera 
else if(window.attachEvent) { // win/ie
  window.attachEvent('onload', rotateNews);
} else { // mac/ie5
  if(typeof window.onload == 'function') {
    var existing = onload;
    window.onload = function() {
      existing();
      rotateNews();
    }
  } else {
    window.onload = function() {
      rotateNews();
    }
  }
}