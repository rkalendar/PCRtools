async function fetchFastaSequence() {
  const id = document.getElementById('ncbiId').value.trim();
  const inputText = document.getElementById('inputText');
  if (!id) { inputText.value = 'Please enter a valid NCBI ID.'; return; }
  // Fetch and display the FASTA sequence
  try {
    // NCBI E-utility endpoint for fetching FASTA
    const url = `https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${id}&rettype=fasta&retmode=text`;
    try {
      inputText.value = '';
      const response = await fetch(url);
      if (!response.ok) {
        throw new Error(`Network response was not ok (Status: ${response.status})`);
      }
      inputText.value = await response.text();
    } catch (error) {
      throw new Error(`Could not fetch FASTA for ID: ${id}\n${error.message}`);
    }
  } catch (error) { inputText.value = `Error: ${error.message}`; }

  const number_of_sequences = document.getElementById('number_of_sequences');
  let result = Display(document.getElementById('inputText').value);
  number_of_sequences.innerHTML = result;
}
