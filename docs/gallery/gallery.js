/* Keep the archive quiet until requested; play one recording at a time. */
(() => {
  const videos = [...document.querySelectorAll('video')];
  for (const video of videos) video.addEventListener('play', () => {
    for (const other of videos) if (other !== video && !other.paused) other.pause();
  });
})();
