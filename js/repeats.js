function analysis() {
    const canvas = document.getElementById('myCanvas');
    const ctx = canvas.getContext('2d');
    let sensitive = document.getElementById('sensitive').checked;
    let lowcomplexity = document.getElementById('lowcomplexity').checked;
    let s = (document.getElementById('inputText').value).toLowerCase();
    let kmer = parseInt(document.getElementById('kmer').value);
    let kmer2 = parseInt(document.getElementById('minlen').value);
    let l = s.length;

    if (l < kmer) { return; }
    if (kmer < 9) { kmer = 9; }
    if (kmer > 21) { kmer = 21; }
    if (kmer2 < kmer) { kmer2 = kmer; }
    const top = (sensitive) ? kmer / 2 - 1 : kmer - 1;

    ReadResult = ReadingSeq(s);
    const name_seq = ReadResult.name_seq;
    const seqs = ReadResult.seqs;
    const n_seq = seqs.length;

    const Kmax = 14; // Kmax=11 for SSR  Kmax=14 for telomers
    const result = ["", ""];
    let resultarea1 = "";
    let resultarea2 = "";
    for (let n = 0; n < n_seq; n++) {
        let startTime = performance.now();

        seqs[n] = DNA(seqs[n]);
        const lseqs = seqs[n].length;
        let msk = RepeatMask(kmer, seqs[n]);
        let f = 0;  // amount of blocks
        const x1 = [];
        const x2 = [];
        for (let i = 0; i < lseqs; i++) {
            if (msk[i] > 1) {
                let z1 = i;
                let z2 = i;
                let v = 0;

                while (i + 1 < i < lseqs && msk[i + 1] > 1) {
                    i++;
                    z2 = i;
                    if (msk[i] > top) {
                        v++;
                    }
                }
                if (v > kmer) {
                    z2 = z2 + kmer - 1;
                    if (z2 >= lseqs) {
                        z2 = lseqs - 1;
                    }
                    x1.push(z1);
                    x2.push(z2);
                    f++;
                }
            }
        }
        // joining close blocks
        for (let i = 0; i < f; i++) {
            if (x2[i] - x1[i] < kmer2) {
                x2.splice(i, 1);
                x1.splice(i, 1);
                i--;
            }
        }

        /*
               // Join blocks x2[i] with x2[i+1] if they are adjusted
               for (let i = 0; i < f - 1; i++) {
                   if (x1[i + 1] - x2[i] <= gap) {
                       x2[i] = x2[i + 1];
                       x2.splice(i + 1, 1);
                       x1.splice(i + 1, 1);
                       i--;
                   }
               }
       */


        const canvasWidth = canvas.width;   // Gets the width of the canvas
        const canvasHeight = canvas.height; // Gets the height of the canvas
        const w = lseqs / (canvasWidth - 20);
        let h = (canvasHeight - 20) / x1.length;
        if (x1.length > 100) { h = (canvasHeight - 20) / 100; }

        // Draw a line based on the dynamic coordinates
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.beginPath();
        ctx.moveTo(10, 10);
        ctx.lineTo(canvasWidth - 10, 10);
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 2; // Set line width
        ctx.stroke();
        ctx.closePath();
        const h1 = (canvasWidth - 20) / 10;
        const z = lseqs / 10;
        for (let i = 0; i < 11; i++) {
            ctx.beginPath();
            const h2 = i * h1;
            ctx.moveTo(h2 + 10, 5);  //x1, y1
            ctx.lineTo(h2 + 10, 15); //x2, y2
            ctx.strokeStyle = 'black';
            ctx.lineWidth = 1; // Set line width
            ctx.stroke();
            ctx.closePath();
            ctx.font = '8px Courier';
            ctx.fillStyle = 'blue';
            ctx.fillText((1 + z * i).toFixed(0), h2 + 5 - i, 5);
        }

        let b = seqs[n].split("");
        let trepeats = 0;
        let ssrrepeats = 0;
        let resultarea3 = "";


        if (lowcomplexity) {
            const lowsequence = new LowComplexitySequence(seqs[n], Kmax, kmer2);
            // console.log(lowsequence.getMapBytes());
            ssrrepeats = lowsequence.getTotalRepeats();
            const bloks = lowsequence.getBlocks();
            const data = bloks[0];
            for (let i = 0; i < data.length; i += 2) {
                resultarea3 += "SSR\t" + (1 + data[i]) + "\t" + (1 + data[i] + data[i + 1]) + "\t" + data[i + 1] + "\n";
                for (let k = data[i]; k < data[i] + data[i + 1]; k++) {
                    b[k] = b[k].toUpperCase();
                }
                const startX = 10 + (data[i] / w);
                let endX = startX + (data[i + 1] / w);
                if (endX - startX < 2) { endX = startX + 2; }
                ctx.beginPath();
                ctx.moveTo(startX, 20);
                ctx.lineTo(endX, 20);
                ctx.strokeStyle = 'green';
                ctx.lineWidth = 3;
                ctx.stroke();
                ctx.closePath();
                ctx.beginPath();
                ctx.moveTo(startX, 10);
                ctx.lineTo(endX, 10);
                ctx.strokeStyle = 'red';
                ctx.lineWidth = 4;
                ctx.stroke();
                ctx.closePath();
            }
        }

        //  if (x1.length > 0) {
        const clustering = new SequencesClustering(seqs[n], x1, x2, 60);
        //  const numberOfClusters = clustering.getNcl();
        //  console.log("Number of clusters:", numberOfClusters);
        //  const clusterAssignments = clustering.result();
        //  console.log("Cluster assignments:", clusterAssignments);
        // Get the processed array of start positions and lengths
        //  const processedArray = clustering.resultArray();
        //  console.log("Processed array:", processedArray);
        // Get the combined and sorted array
        const sortedCombinedArray = clustering.getSortedCombinedArray();
        //  console.log(sortedCombinedArray);
        sortedCombinedArray.forEach((item, index) => {
            const clusterId = item.cluster; // Get the cluster ID
            const data = item.sequenceData; // Get the sequence data (array [startIndex, length])
            trepeats += data[1];
            for (let k = data[0]; k < data[0] + data[1]; k++) {
                b[k] = b[k].toUpperCase();
            }
            resultarea3 += clusterId + "\t" + (1 + data[0]) + "\t" + (1 + data[0] + data[1]) + "\t" + data[1] + "\n";
            // console.log(`Cluster ID: ${clusterId}, Data: Start Index - ${data[0]}, Length - ${data[1]}`);
            // Generate dynamic coordinates for line positions based on data and canvas dimensions
            const startX = 10 + (data[0] / w);           // Scale based on canvas width (assuming data[0] is between 0-100)
            const startY = 20 + clusterId * h;           // Start Y position at 30% of the canvas height
            let endX = startX + (data[1] / w);           // Scale line length based on canvas width and data[1]                
            const endY = 20 + clusterId * h;             // End Y position at 70% of the canvas height
            if (endX - startX < 2) { endX = startX + 2; }

            // Draw a line based on the dynamic coordinates
            if (clusterId < 100) {
                ctx.beginPath();
                ctx.moveTo(startX, startY);
                ctx.lineTo(endX, endY);
                ctx.strokeStyle = 'blue';
                ctx.lineWidth = 2; // Set line width
                ctx.stroke();
                ctx.closePath();
            }
            ctx.beginPath();
            ctx.moveTo(startX, 10);
            ctx.lineTo(endX, 10);
            ctx.strokeStyle = 'red';
            ctx.lineWidth = 4; // Set line width
            ctx.stroke();
            ctx.closePath();
        });

        // End time
        let endTime = performance.now();
        let runtime = endTime - startTime;

        if (trepeats > 0 || ssrrepeats > 0) {
            resultarea1 += "\n>" + name_seq[n] + ": " + lseqs + " bp\n";
            resultarea1 += b.join("");
            resultarea2 += name_seq[n] + ": " + lseqs + " bp\n" + "Sequence coverage by repeats =" + (trepeats * 100 / lseqs).toFixed(2) + "%\n";
            if (lowcomplexity) { resultarea2 += "Sequence coverage by SSR repeats =" + (ssrrepeats * 100 / lseqs).toFixed(2) + "%\n"; }
            if (runtime > 999) {
                runtime = runtime / 1000;
                resultarea2 += "The code took " + runtime.toFixed(0) + " seconds to run.\n\n";
            } else {
                resultarea2 += "The code took " + runtime.toFixed(0) + " milliseconds to run.\n\n";
            }
            resultarea2 += "CltrID\tstart\tend\tlength\n";
            resultarea2 += resultarea3;
        }
        resultarea2 += "\n";
        //  }
    }
    result[0] += resultarea1;
    result[1] += resultarea2;
    return result;
}

function RepeatMask(h, r1) {
    const l = r1.length;
    const r2 = ComplementDNA(r1);
    const m = new Array(l).fill(0);
    const d = new Map();
    for (let y = 0; y < r1.length - h + 1; y++) {
        let s = r1.substring(y, y + h);
        if (!s.includes('n')) {
            if (d.has(s)) {
                let k = d.get(s);
                d.set(s, ++k);
            } else {
                d.set(s, 1);
            }
        }
    }
    for (let y = 0; y < r1.length - h + 1; y++) {
        let s = r2.substring(y, y + h);
        if (!s.includes('n')) {
            if (d.has(s)) {
                let k = d.get(s);
                d.set(s, ++k);
            }
        }
    }
    for (let y = 0; y < r1.length - h + 1; y++) {
        let s = r1.substring(y, y + h);
        if (!s.includes('n')) {
            const k = d.get(s);
            if (k > 1) {
                for (var i = y; i < y + h; i++) {
                    m[i]++;
                }
            }
        }
    }
    for (let y = 0; y < r1.length - h + 1; y++) {
        let s = r2.substring(y, y + h);
        if (!s.includes('n')) {
            const k = d.get(s);
            if (k > 1) {
                for (var i = l - y - h; i < l - y; i++) {
                    m[i]++;
                }
            }
        }
    }
    return m;
}