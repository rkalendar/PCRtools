document.addEventListener('DOMContentLoaded', () => {
  const fileInput      = document.getElementById('fileInput');
  const fileContent    = document.getElementById('inputText');
  const seqCounterNode = document.getElementById('number_of_sequences');
  const CHUNK_SIZE     = 1 << 21;  //const CHUNK_SIZE = 2 * 1024 * 1024;   // 2MB
  const decoder        = new TextDecoder();

  fileInput.addEventListener('change', async (event) => {
    const file = event.target.files[0];
    if (!file) return;

    try {
      fileContent.value = '';   
      let offset = 0;

      while (offset < file.size) {
        const slice   = file.slice(offset, offset + CHUNK_SIZE);
        const buffer  = await slice.arrayBuffer();  
        fileContent.value += decoder.decode(buffer, { stream: true });

        offset += CHUNK_SIZE;
 
        await new Promise(r => requestAnimationFrame(r));
      }
    } catch (err) {
      console.error(err);
      fileContent.value = 'Error reading file';
    } finally {
      console.log('File loading process complete.');
      seqCounterNode.innerHTML = Display(fileContent.value);
    }
  });
});
