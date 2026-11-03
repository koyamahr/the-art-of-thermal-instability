/* Motion is optional. Reduced-motion and save-data visitors start with a still. */
(() => {
  const video = document.getElementById('thermal-film');
  const button = document.getElementById('film-toggle');
  if (!video || !button) return;
  const preference = matchMedia('(prefers-reduced-motion: reduce)');
  const sync = () => { button.textContent=video.paused?'Play simulation':'Pause simulation';button.setAttribute('aria-pressed',String(!video.paused)); };
  const play = () => { video.play().catch(sync); };
  button.addEventListener('click',()=>video.paused?play():video.pause());
  video.addEventListener('play',sync);video.addEventListener('pause',sync);
  const unavailable=()=>{button.textContent='Video unavailable';button.disabled=true;};
  video.addEventListener('error',unavailable);video.querySelector('source').addEventListener('error',unavailable);
  preference.addEventListener('change',event=>{if(event.matches)video.pause();});
  if(!preference.matches && !navigator.connection?.saveData) play();
})();
