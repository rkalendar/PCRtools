async function fetchVariants() {
  let flankSize = parseInt(document.getElementById('flankInput').value);
  const species = document.getElementById('speciesInput').value.trim();
  const variantIDs = document.getElementById('variantInput').value.trim();
  const output = document.getElementById('inputText');
  output.value = '';

  if (isNaN(flankSize) || flankSize < 10) { flankSize = 200; }
  if (!species) { species.value = 'homo_sapiens'; }
  if (!species) { output.value = '[ERROR] Please enter a species, for example: homo_sapiens \n'; return; }
  if (!variantIDs) { output.value = '[ERROR] Please enter at least one variant srID, for example: rs1357617 \n'; return; }

  // Split multiple IDs by comma or whitespace
  let ids = variantIDs.split(/[\s,]+/).map(s => s.trim()).filter(Boolean);
  if (!ids.length) { output.value = '[ERROR] No valid variant srIDs found.\n'; return; }
  //output.value = `Processing ${ids.length} variants for species="${species}" with flank=${flankSize}...\n\n`;

  for (const id of ids) {
    try {
      //  Fetch Variation info from Ensembl e.g. https://rest.ensembl.org/variation/homo_sapiens/rs56116432?content-type=application/json
      const varUrl = `https://rest.ensembl.org/variation/${species}/${id}?content-type=application/json`;
      const varResp = await fetch(varUrl);
      if (!varResp.ok) {
        throw new Error(`Variation fetch error: ${varResp.status} (${varResp.statusText})`);
      }
      const varData = await varResp.json();
      // Check if there's a "mappings" array
      if (!varData.mappings || !varData.mappings.length) {
        output.value += `[INFO] No mappings found for "${id}".\n\n`;
        continue;
      }
      // For simplicity, pick the first mapping
      const map = varData.mappings[0];
      const chrom = map.seq_region_name; // e.g. "9"
      const start = map.start;
      const end = map.end;
      const strand = map.strand; // 1 or -1
      const allelesString = map.allele_string;
      //Build region with flank  If the variant is from start..end, we assume a reference allele length of (end-start+1). We'll fetch that region ± flank.
      const refLen = (end - start + 1);
      let seqStart = start - flankSize;
      if (seqStart < 1) seqStart = 1;
      let seqEnd = end + flankSize;
      // Fetch the partial FASTA from Ensembl e.g. https://rest.ensembl.org/sequence/region/<species>/<chr>:<start>..<end>:<strand>
      const seqUrl = `https://rest.ensembl.org/sequence/region/${species}/${chrom}:${seqStart}..${seqEnd}:${strand}?content-type=text/plain`;
      const seqResp = await fetch(seqUrl);
      if (!seqResp.ok) { throw new Error(`Sequence fetch error: ${seqResp.status} (${seqResp.statusText})`); }
      let rawSeq = await seqResp.text();
      rawSeq = rawSeq.replace(/\n/g, '').toUpperCase(); // remove line breaks
      // zero-based index of the variant region
      const localIndex = (start - seqStart);
      // if variant length is > 1:
      const localEnd = localIndex + (refLen - 1);
      if (localIndex < 0) { localIndex = 0; }
      if (localEnd >= rawSeq.length) { localEnd = rawSeq.length - 1; }
      // Extract the reference substring
      const seqRef = rawSeq.substring(localIndex, localEnd + 1);
      // Build bracket: prefix + [REF->ALT] + suffix
      const prefix = rawSeq.substring(0, localIndex).toLowerCase() + " ";
      const suffix = " " + rawSeq.substring(localEnd + 1).toLowerCase();
      // Replace the reference substring with [REF->ALT]
      const bracketed = `${prefix}[${allelesString}]${suffix}`;
      // We'll create a simple FASTA-like header line e.g. ">chr9:1000-1010 (rs123) REF->ALT"
      const header = `>${chrom}:${start}-${end}(${strand}) ${id}`;
      output.value += `${header}\n${bracketed}\n\n`;
    } catch (err) { output.value += `[ERROR] For variant="${id}": ${err.message}\n\n`; }
  }
  const number_of_sequences = document.getElementById('number_of_sequences');
  let result = Display(document.getElementById('inputText').value);
  number_of_sequences.innerHTML = result;
}
