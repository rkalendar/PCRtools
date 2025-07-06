/*
// Example DNA sequence string
const seq = "aatcgaatgtccggtatccttcgagcctgctagcgtacgatagctagctagtcgtgacgtacgt";

// Example array x with pairs of integers
const x = [0, 10, 20, 30, 40, 50, 60, 70]; // Replace with actual data

// Desired similarity value
const similarity = 70; // Choose a similarity between 60 and 80

// Create an instance of the SequencesClustering class
const clustering = new SequencesClustering(seq, x, similarity);

// Get the number of clusters
const numberOfClusters = clustering.getNcl();
console.log("Number of clusters:", numberOfClusters);

// Get the cluster assignments for each sequence
const clusterAssignments = clustering.result();
console.log("Cluster assignments:", clusterAssignments);

// Get the processed array of start positions and lengths
const processedArray = clustering.resultArray();
console.log("Processed array:", processedArray);

*/

class SequencesClustering {
    constructor(seq, x1, x2, similarity) {
        const nseq = x1.length;

        // Clamp similarity between 60 and 80
        if (similarity < 60) similarity = 60;
        if (similarity > 80) similarity = 80;

        this.d = new Array(nseq).fill(0).map(() => new Array(2));
        for (let j = 0; j < nseq; j++) {
            this.d[j][0] = x1[j];
            this.d[j][1] = 1 + x2[j] - x1[j];
        }
        this.d.sort((a, b) => b[1] - a[1]);
        this.spd = 21;    // =16-32 spd parameter
        this.dif = 1.35;  //  dif parameter (deviation)
        this.ncl = this.clustering(seq, nseq, similarity);
    }

    clustering(seq, nseq, sim) {
        // Initialize the pattern map
        const pt = {
            'aatt': 0, 'acgt': 1, 'agct': 2, 'ttaa': 3, 'tgca': 4, 'tcga': 5, 'ccgg': 6, 'catg': 7, 'ctag': 8, 'ggcc': 9, 'gatc': 10, 'gtac': 11
            //'aatt': 0, 'acgt': 1, 'agct': 2, 'ttaa': 3, 'tgca': 4, 'tcga': 5, 'ccgg': 6, 'catg': 7, 'ctag': 8, 'ggcc': 9, 'gatc': 10, 'gtac': 11, 'atta': 12, 'acca': 13, 'agga': 14, 'taat': 15, 'tggt': 16, 'tcct': 17, 'caac': 18, 'cttc': 19, 'cggc': 20, 'gaag': 21, 'gttg': 22, 'gccg': 23
            //'aatc': 0, 'aatg': 1, 'aact': 2, 'aacg': 3, 'aagt': 4, 'aagc': 5,'atac': 6, 'atag': 7, 'attc': 8, 'attg': 9, 'atca': 10, 'atct': 11,'atcc': 12, 'atcg': 13, 'atga': 14, 'atgt': 15, 'atgc': 16, 'atgg': 17,'acat': 18, 'acag': 19, 'acta': 20, 'actt': 21, 'actc': 22, 'actg': 23,'acct': 24, 'accg': 25, 'acga': 26, 'acgc': 27, 'acgg': 28, 'agat': 29,'agac': 30, 'agta': 31, 'agtt': 32, 'agtc': 33, 'agtg': 34, 'agca': 35,'agcc': 36, 'agcg': 37, 'aggt': 38, 'aggc': 39, 'taac': 40, 'taag': 41,'tatc': 42, 'tatg': 43, 'taca': 44, 'tact': 45, 'tacc': 46, 'tacg': 47,'taga': 48, 'tagt': 49, 'tagc': 50, 'tagg': 51, 'ttac': 52, 'ttag': 53,'ttca': 54, 'ttcg': 55, 'ttga': 56, 'ttgc': 57, 'tcaa': 58, 'tcat': 59,'tcac': 60, 'tcag': 61, 'tcta': 62, 'tctg': 63, 'tcca': 64, 'tccg': 65,'tcgt': 66, 'tcgc': 67, 'tcgg': 68, 'tgaa': 69, 'tgat': 70, 'tgac': 71,'tgag': 72, 'tgta': 73, 'tgtc': 74, 'tgct': 75, 'tgcc': 76, 'tgcg': 77,'tgga': 78, 'tggc': 79, 'caat': 80, 'caag': 81, 'cata': 82, 'catt': 83,'catc': 84, 'cact': 85, 'cacg': 86, 'caga': 87, 'cagt': 88, 'cagc': 89,'cagg': 90, 'ctaa': 91, 'ctat': 92, 'ctac': 93, 'ctta': 94, 'cttg': 95,'ctca': 96, 'ctcg': 97, 'ctga': 98, 'ctgt': 99, 'ctgc': 100, 'ctgg': 101,'ccat': 102, 'ccag': 103, 'ccta': 104, 'cctg': 105, 'ccga': 106, 'ccgt': 107,'cgaa': 108, 'cgat': 109, 'cgac': 110, 'cgag': 111, 'cgta': 112, 'cgtt': 113,'cgtc': 114, 'cgtg': 115, 'cgca': 116, 'cgct': 117, 'cgga': 118, 'cggt': 119,'gaat': 120, 'gaac': 121, 'gata': 122, 'gatt': 123, 'gatg': 124, 'gaca': 125,'gact': 126, 'gacc': 127, 'gacg': 128, 'gagt': 129, 'gagc': 130, 'gtaa': 131,'gtat': 132, 'gtag': 133, 'gtta': 134, 'gttc': 135, 'gtca': 136, 'gtct': 137,'gtcc': 138, 'gtcg': 139, 'gtga': 140, 'gtgc': 141, 'gcaa': 142, 'gcat': 143,'gcac': 144, 'gcag': 145, 'gcta': 146, 'gctt': 147, 'gctc': 148, 'gctg': 149,'gcca': 150, 'gcct': 151, 'gcga': 152, 'gcgt': 153, 'ggat': 154, 'ggac': 155,'ggta': 156, 'ggtc': 157, 'ggca': 158, 'ggct': 159
            //'aaatt': 0,'aagtt': 1,'ataat': 2,'atgat': 3,'acagt': 4,'acggt': 5,'agact': 6,'aggct': 7,'ttaaa': 8,'ttgaa': 9,'taata': 10,'tagta': 11,'tgaca': 12,'tggca': 13,'tcaga': 14,'tcgga': 15,'ccagg': 16,'ccggg': 17,'caatg': 18,'cagtg': 19,'ctaag': 20,'ctgag': 21,'cgacg': 22,'cggcg': 23,'ggacc': 24,'gggcc': 25,'gaatc': 26,'gagtc': 27,'gtaac': 28,'gtgac': 29,'gcagc': 30,'gcggc': 31
            //'atttta': 0,'aattaa': 1,'ttaatt': 2,'taaaat': 3,'atccta': 4,'atggta': 5,'agttga': 6,'aaccaa': 7,'aaggaa': 8,'ttcctt': 9,'ttggtt': 10,'tcaact': 11,'tgaagt': 12,'taccat': 13,'taggat': 14,'ctaatc': 15,'cttttc': 16,'caaaac': 17,'cattac': 18,'gtaatg': 19,'gttttg': 20,'gaaaag': 21,'gattag': 22,'acttca': 23,'acccca': 24,'acggca': 25,'agccga': 26,'agggga': 27,'tcccct': 28,'tcggct': 29,'tgccgt': 30,'tggggt': 31,'ctggtc': 32,'ccaacc': 33,'ccttcc': 34,'cgaagc': 35,'cgttgc': 36,'caggac': 37,'gtcctg': 38,'gcaacg': 39,'gcttcg': 40,'ggaagg': 41,'ggttgg': 42,'gaccag': 43,'ccggcc': 44,'cggggc': 45,'gccccg': 46,'ggccgg': 47
            //'aaattt': 0,'aacgtt': 1,'aagctt': 2,'aatatt': 3,'acatgt': 4,'accggt': 5,'acgcgt': 6,'actagt': 7,'agatct': 8,'agcgct': 9,'aggcct': 10,'agtact': 11,'atcgat': 12,'atgcat': 13,'attaat': 14,'caattg': 15,'cacgtg': 16,'cagctg': 17,'catatg': 18,'ccatgg': 19,'cccggg': 20,'ccgcgg': 21,'cctagg': 22,'cgatcg': 23,'cggccg': 24,'cgtacg': 25,'ctatag': 26,'ctcgag': 27,'ctgcag': 28,'cttaag': 29,'gaattc': 30,'gacgtc': 31,'gagctc': 32,'gatatc': 33,'gcatgc': 34,'gccggc': 35,'gcgcgc': 36,'gctagc': 37,'ggatcc': 38,'ggcgcc': 39,'gggccc': 40,'ggtacc': 41,'gtatac': 42,'gtcgac': 43,'gtgcac': 44,'gttaac': 45,'taatta': 46,'tacgta': 47,'tagcta': 48,'tcatga': 49,'tccgga': 50,'tcgcga': 51,'tctaga': 52,'tgatca': 53,'tgcgca': 54,'tggcca': 55,'tgtaca': 56,'ttataa': 57,'ttcgaa': 58,'ttgcaa': 59,'tttaaa': 60
        };
        /*
        const kmers = [
            "aaatt", "aagtt", "acagt", "acggt", "agact", "aggct", "ataat", "atgat",
            "caatg", "cagtg", "ccagg", "ccggg", "cgacg", "cggcg", "ctaag", "ctgag",
            "gaatc", "gagtc", "gcagc", "gcggc", "ggacc", "gggcc", "gtaac", "gtgac",
            "taata", "tagta", "tcaga", "tcgga", "tgaca", "tggca", "ttaaa", "ttgaa",
            "aactt", "aattt", "accgt", "actgt", "agcct", "agtct", "atcat", "attat",
            "cactg", "cattg", "cccgg", "cctgg", "cgccg", "cgtcg", "ctcag", "cttag",
            "gactc", "gattc", "gccgc", "gctgc", "ggccc", "ggtcc", "gtcac", "gttac",
            "tacta", "tatta", "tccga", "tctga", "tgcca", "tgtca", "ttcaa", "tttaa",
            "aaatc", "aagtc", "ccagt", "acgga", "agaca", "atatc", "atgct", "cagac",
            "ccata", "ccgat", "cgaac", "cgaat", "ctaac", "gaccg", "gaccc", "gccat",
            "gcgca", "ggtta", "ggcat", "gtatt", "gtgga", "tacac", "taggc", "tcaaa",
            "tcgcc", "tggat", "tgtta", "ttacc", "atcta", "atgta", "actca", "acgca",
            "agtga", "agcga", "tacat", "tagat", "tgagt", "tgcgt", "tcact", "tcgct",
            "catac", "ctatc", "ctgtc", "cgagc", "cgtgc", "gatag", "gacag", "gtatg",
            "gtctg", "gcacg", "gctcg"
          ];
         
        const pt = kmers.reduce((acc, curr, idx) => {
          acc[curr] = idx;
          return acc;
        }, {});
           */


        //const keys = Object.keys(pt);
        const kmer = 4; //keys[0].length;
        const nkmers = Object.keys(pt).length;
        const m2 = Array.from({ length: nseq }, () => Array(nkmers).fill(0));

        // Filling m2 matrix with counts
        for (let j = 0; j < nseq; j++) {
            for (let i = 0; i < this.d[j][1] - kmer + 1; i++) {
                let s = seq.substring(this.d[j][0] + i, this.d[j][0] + i + kmer);
                if (pt.hasOwnProperty(s)) {
                    m2[j][pt[s]]++;
                }
                //   s = ComplementDNA(s);
                //   if (pt.hasOwnProperty(s)) {
                //       m2[j][pt[s]]++;
                //   }
            }
            for (let i = 0; i < nkmers; i++) {
                if (m2[j][i] > 0) {
                    m2[j][i] = this.d[j][1] / (1 + m2[j][i]);
                }
            }

        }

        let n = 0;                      // Number of clusters
        this.cx = Array(nseq).fill(0);  // Cluster assignment array

        for (let i = 0; i < nseq; i++) {
            if (this.cx[i] === 0) {
                n++;
                this.cx[i] = n;

                for (let j = i + 1; j < nseq; j++) {
                    if (this.cx[j] === 0) {
                        const m = Array(nkmers + 1).fill(0);

                        for (let k = 0; k < nkmers; k++) {
                            if (m2[i][k] > 1 && m2[j][k] > 1) {
                                if (m[0] === 0) {
                                    m[0] = 1;
                                    m[m[0]] = k;
                                } else {
                                    m[0]++;
                                    m[m[0]] = k;
                                }
                            }
                        }

                        let v = 0; // Matches
                        let z = 0; // Theoretical max matches

                        for (let k = 1; k < m[0]; k++) {
                            for (let y = k + 1; y < 1 + m[0]; y++) {
                                let di, dj;
                                z++;
                                if (m2[i][m[k]] < m2[i][m[y]]) {
                                    di = (100 * m2[i][m[k]]) / m2[i][m[y]];
                                    dj = (100 * m2[j][m[k]]) / m2[j][m[y]];
                                } else {
                                    di = (100 * m2[i][m[y]]) / m2[i][m[k]];
                                    dj = (100 * m2[j][m[y]]) / m2[j][m[k]];
                                }
                                if (di === dj) { v++; }
                                else {
                                    if (di > dj) {
                                        if (di <= dj * this.dif) v++;
                                    }
                                    else {
                                        if (dj <= di * this.dif) v++;
                                    }
                                }
                            }
                            if (z > this.spd) break;
                        }

                        if (v > 1 && ((100 * v) / z) > sim) {
                            this.cx[j] = n;
                        }
                    }
                }
            }
        }
        return n;
    }

    getSortedCombinedArray() {
        // Combine cluster assignments with sequence data into a single array of objects
        const combined = this.cx.map((cluster, index) => ({
            cluster: cluster,
            sequenceData: this.d[index]
        }));
        // Sort the combined array by cluster assignments
        combined.sort((a, b) => a.cluster - b.cluster);
        return combined;
    }

    getNcl() {
        return this.ncl;
    }

    result() {
        return this.cx;
    }

    resultArray() {
        return this.d;
    }
}