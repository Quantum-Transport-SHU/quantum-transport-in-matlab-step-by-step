// a script for runing code on http://octave-online.net/ unlimitly

var el = document.getElementById('add_time_container');
var mother = document.getElementById('runtime_controls_container')

setInterval(()=>{
  if(el.getAttribute('aria-hidden') === 'false' && mother.getAttribute('aria-hidden') === 'false') {
    console.log('going to click')
    el.click()
  }
},1000)