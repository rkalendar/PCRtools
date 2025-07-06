async function fetchMultipleRs() {
  // const speciesName = "Homo sapiens";
  const rsIDsInput = document.getElementById('rsIDs').value.trim();
  const output = document.getElementById('inputText');
  let flankSize = parseInt(document.getElementById('flankInput').value);

  output.value = '';
  if (!rsIDsInput) { output.value = '[ERROR] Please enter at least one rsID, for example: rs1357617\n'; return; }
  if (isNaN(flankSize) || flankSize < 10) { flankSize = 200; }

  // Split multiple rsIDs
  let rsIDList = rsIDsInput.split(/[\s,]+/).map(s => s.trim()).filter(Boolean);
  if (!rsIDList.length) { output.value = '[ERROR] No valid rsIDs found.\n'; return; }
  for (const rsID of rsIDList) {
    //   output.value += `=== RSID: ${rsID} ===\n`;
    try {
      // Fetch dbSNP XML
      const dbSnpUrl = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=${rsID}&retmode=xml`;
      const snpResp = await fetch(dbSnpUrl);
      if (!snpResp.ok) {
        throw new Error(`dbSNP fetch error: ${snpResp.status} (${snpResp.statusText})`);
      }
      const xmlText = await snpResp.text();

      // Parse the XML
      const parser = new DOMParser();
      const xmlDoc = parser.parseFromString(xmlText, 'application/xml');
      if (xmlDoc.querySelector('parsererror')) {
        throw new Error('Error parsing dbSNP XML. Possibly invalid rsID.');
      }

      // [Optional] Parse <Rs taxId="..."> to see actual species tax ID
      const rsElem = xmlDoc.querySelector('Rs');
      if (rsElem && rsElem.hasAttribute('taxId')) {
        const taxId = rsElem.getAttribute('taxId');
        output.value += `taxId from <Rs>: ${taxId}\n`;
      }

      // <SPDI> may have multiple comma-separated notations
      const spdiNode = xmlDoc.querySelector('SPDI');
      if (!spdiNode) {
        output.value += `[INFO] No <SPDI> found for rsID=${rsID}.\n\n`;
        continue;
      }
      const spdiText = spdiNode.textContent.trim();
      const spdiItems = spdiText.split(',').map(s => s.trim()).filter(Boolean);
      if (!spdiItems.length) {
        output.value += `[INFO] <SPDI> empty/invalid for rsID=${rsID}: "${spdiText}"\n\n`;
        continue;
      }

      // <ALLELE> tags (REF/ALT)
      const alleleNodes = xmlDoc.querySelectorAll('ALLELE');
      let alleleList = Array.from(alleleNodes).map(n => n.textContent.trim());


      // For each SPDI, fetch partial FASTA from nuccore
      for (const spdiItem of spdiItems) {
        //   output.value += `--- SPDI: ${spdiItem} ---\n`;
        const parts = spdiItem.split(':'); // e.g. ["NC_000009.12","133256041","C","A"]
        if (parts.length !== 4) {
          output.value += `[ERROR] Invalid SPDI format: "${spdiItem}"\n\n`;
          continue;
        }
        const [acc, posStr, spdiRef, spdiAlt] = parts;
        const position = parseInt(posStr, 10);
        if (!acc || isNaN(position) || !spdiRef || !spdiAlt) {
          output.value += `[ERROR] Could not parse SPDI fields from "${spdiItem}".\n\n`;
          continue;
        }

        // We'll fetch region: position ± flank
        const refLen = spdiRef.length;
        let seqStart = position - flankSize;
        if (seqStart < 1) seqStart = 1;
        let seqStop = position + refLen - 1 + flankSize;

        // E-fetch partial FASTA
        const efetchUrl = new URL('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi');
        efetchUrl.searchParams.set('db', 'nuccore');
        efetchUrl.searchParams.set('id', acc);
        efetchUrl.searchParams.set('rettype', 'fasta');
        efetchUrl.searchParams.set('retmode', 'text');
        efetchUrl.searchParams.set('seq_start', seqStart.toString());
        efetchUrl.searchParams.set('seq_stop', seqStop.toString());
        // optionally: efetchUrl.searchParams.set('api_key','YOUR_KEY');

        try {
          const fastaResp = await fetch(efetchUrl);
          if (!fastaResp.ok) {
            throw new Error(`FASTA fetch error: ${fastaResp.status} (${fastaResp.statusText})`);
          }
          const fastaText = await fastaResp.text();

          // Split out header + sequence
          const lines = fastaText.trim().split('\n');
          const originalHeader = lines[0] || '>No FASTA header from NCBI';
          const seq = lines.slice(1).join('').toUpperCase();

          // zero-based index in this chunk
          const localIndex = position - seqStart;
          const localEnd = localIndex + refLen - 1;
          if (localIndex < 0 || localEnd >= seq.length) {
            output.value += `[ERROR] Position out of range in fetched FASTA.\n\n`;
            continue;
          }

          const seqRef = seq.substring(localIndex, localEnd + 1);

          // For each <ALLELE> (REF/ALT), or fallback spdiRef->spdiAlt
          const fallbackAllele = `${spdiRef}/${spdiAlt}`;
          const finalAlleles = (alleleList.length === 0)
            ? [fallbackAllele]
            : alleleList;

          let alleleIndex = 0;
          for (const alleleText of finalAlleles) {
            alleleIndex++;
            const [alleleRef, alleleAlt] = alleleText.split('/');
            const refUp = alleleRef.toUpperCase();

            // bracket annotation
            const prefix = seq.substring(0, localIndex).toLowerCase() + " ";
            const suffix = " " + seq.substring(localEnd + 1).toLowerCase();
            const bracketed = `${prefix}[${refUp}]${suffix}`;

            // Build FASTA-like output  Keep the original NCBI header, but add annotation
            const annotatedHeader = `${originalHeader} (rsID:${rsID})`;
            output.value += `${annotatedHeader}\n${bracketed}\n\n`;
          }

        } catch (fastaErr) {
          output.value += `[ERROR] Could not fetch FASTA for SPDI="${spdiItem}": ${fastaErr.message}\n\n`;
        }
        break;
      }

    } catch (err) {
      output.value += `[ERROR] For rsID="${rsID}": ${err.message}\n\n`;
    }
  }

  const number_of_sequences = document.getElementById('number_of_sequences');
  let result = Display(document.getElementById('inputText').value);
  number_of_sequences.innerHTML = result;
}