function analysis(w1) {
        const M5w = new Array(100).fill(0); //new Array(100);
        const M3w = new Array(100).fill(0); //new Array(100);
        const Mxw = new Array(100).fill(0); //new Array(100);
        const M5e = new Array(100).fill(0); //new Array(100);
        const M3e = new Array(100).fill(0); //new Array(100);
        const Mxe = new Array(100).fill(0); //new Array(100);
        const M5v = new Array(100).fill(0); //new Array(100);
        const M3v = new Array(100).fill(0); //new Array(100);
        const Mxv = new Array(100).fill(0); //new Array(100);

        M5w[1] = 247.2;
        M5e[1] = 0;
        M5v[1] = -97;
        M3w[1] = 179.2;
        M3e[1] = 0;
        M3v[1] = -97;
        M5w[2] = 263.3;
        M5e[2] = 0;
        M5v[2] = -97;
        M5w[3] = 346.97;
        M5e[3] = 0;
        M5v[3] = -97;
        M3w[2] = 136.97;
        M3e[2] = 0;
        M3v[2] = -97;
        M5w[4] = 179.2;
        M5e[4] = 0;
        M5v[4] = -97;
        M5w[5] = 458.4;
        M5e[5] = 8700;
        M5v[5] = -97;
        Mxw[1] = 458.4;
        Mxe[1] = 8700;
        Mxv[1] = -97;
        M3w[3] = 192.97;
        M3e[3] = 0;
        M3v[3] = -97;
        M5w[6] = 136.97;
        M5e[6] = 0;
        M5v[6] = -97;
        M3w[4] = 136.97;
        M3e[4] = 0;
        M3v[4] = -97;
        Mxw[2] = 279.2;
        Mxe[2] = 13700;
        Mxv[2] = 4;
        M5w[7] = 695.8;
        M5e[7] = 71000;
        M5v[7] = -97;
        M3w[5] = 695.6;
        M3e[5] = 71000;
        M3v[5] = -97;
        M5w[8] = 787.5;
        M5e[8] = 81000;
        M5v[8] = -97;
        M3w[6] = 787.8;
        M3e[6] = 81000;
        M3v[6] = -97;
        M5w[9] = 1121.9;
        M5e[9] = 104000;
        M5v[9] = -97;
        M3w[7] = 1121.5;
        M3e[7] = 104000;
        M3v[7] = -97;
        M5w[10] = 883.9;
        M5e[10] = 73000;
        M5v[10] = -97;
        M3w[8] = 883.9;
        M3e[8] = 73000;
        M3v[8] = -97;
        M5w[11] = 1020.7;
        M5e[11] = 239000;
        M5v[11] = -97;
        M3w[9] = 1020.2;
        M3e[9] = 239000;
        M3v[9] = -97;
        Mxw[3] = 1020.7;
        Mxe[3] = 239000;
        Mxv[3] = -97;
        M5w[12] = 941.2;
        M5e[12] = 132000;
        M5v[12] = -97;
        M3w[10] = 941.7;
        M3e[10] = 132000;
        M3v[10] = -97;
        M5w[13] = 1047.6;
        M5e[13] = 240000;
        M5v[13] = -97;
        M3w[11] = 1047.2;
        M3e[11] = 240000;
        M3v[11] = -97;
        M5w[14] = 405.4;
        M5e[14] = 0;
        M5v[14] = -97;
        M3w[12] = 437.4;
        M3e[12] = 0;
        M3v[12] = -97;
        M5w[15] = 684.7;
        M5e[15] = 8700;
        M5v[15] = -97;
        Mxw[4] = 684.7;
        Mxe[4] = 8700;
        Mxv[4] = -97;
        M5w[16] = 569.6;
        M5e[16] = 0;
        M5v[16] = -97;
        M3w[13] = 569.6;
        M3e[13] = 0;
        M3v[13] = -97;
        M3w[14] = 556.5;
        M3e[14] = 34000;
        M3v[14] = -97;
        M5w[17] = 554.5;
        M5e[17] = 34000;
        M5v[17] = -97;
        M3w[15] = 554.5;
        M3e[15] = 34000;
        M3v[15] = -97;
        M3w[16] = 556.5;
        M3e[16] = 38000;
        M3v[16] = -97;
        M3w[17] = 556.5;
        M3e[17] = 42700;
        M3v[17] = -97;
        M5w[18] = 871;
        M5e[18] = 0;
        M5v[18] = -97;
        M5w[19] = 597.6;
        M5e[19] = 0;
        M5v[19] = -97;
        M5w[20] = 724.6;
        M5e[20] = 79000;
        M5v[20] = -97;
        M3w[18] = 724.6;
        M3e[18] = 79000;
        M3v[18] = -97;
        Mxw[5] = 724.6;
        Mxe[5] = 79000;
        Mxv[5] = -97;
        M5w[21] = 724.6;
        M5e[21] = 62000;
        M5v[21] = -97;
        M3w[19] = 724.6;
        M3e[19] = 62000;
        M3v[19] = -97;
        Mxw[6] = 724.6;
        Mxe[6] = 62000;
        Mxv[6] = -97;
        M5w[22] = 724.6;
        M5e[22] = 97000;
        M5v[22] = -97;
        M3w[20] = 724.6;
        M3e[20] = 97000;
        M3v[20] = -97;
        Mxw[7] = 724.6;
        Mxe[7] = 97000;
        Mxv[7] = -97;
        M5w[23] = 724.6;
        M5e[23] = 142000;
        M5v[23] = -97;
        M3w[21] = 724.6;
        M3e[21] = 142000;
        M3v[21] = -97;
        Mxw[8] = 724.6;
        Mxe[8] = 142000;
        Mxv[8] = -97;
        M5w[24] = 724.6;
        M5e[24] = 830000;
        M5v[24] = -97;
        M3w[22] = 724.6;
        M3e[22] = 830000;
        M3v[22] = -97;
        Mxw[9] = 724.6;
        Mxe[9] = 830000;
        Mxv[9] = -97;
        M5w[25] = 724.6;
        M5e[25] = 136000;
        M5v[25] = -97;
        M3w[23] = 724.6;
        M3e[23] = 136000;
        M3v[23] = -97;
        Mxw[10] = 724.6;
        Mxe[10] = 136000;
        Mxv[10] = -97;
        M5w[26] = 724.6;
        M5e[26] = 101000;
        M5v[26] = -97;
        M3w[24] = 724.6;
        M3e[24] = 101000;
        M3v[24] = -97;
        M5w[27] = 724.6;
        M5e[27] = 101000;
        M5v[27] = -97;
        M3w[25] = 724.6;
        M3e[25] = 101000;
        M3v[25] = -97;
        Mxw[11] = 724.6;
        Mxe[11] = 101000;
        Mxv[11] = -97;
        M5w[28] = 724.6;
        M5e[28] = 70000;
        M5v[28] = -97;
        M3w[26] = 724.6;
        M3e[26] = 70000;
        M3v[26] = -97;
        Mxw[12] = 724.6;
        Mxe[12] = 70000;
        Mxv[12] = -97;
        M5w[29] = 724.6;
        M5e[29] = 70000;
        M5v[29] = -97;
        M3w[27] = 724.6;
        M3e[27] = 70000;
        M3v[27] = -97;
        Mxw[13] = 724.6;
        Mxe[13] = 70000;
        Mxv[13] = -97;
        M5w[30] = 724.6;
        M5e[30] = 56000;
        M5v[30] = -97;
        M3w[28] = 724.6;
        M3e[28] = 56000;
        M3v[28] = -97;
        Mxw[14] = 724.6;
        Mxe[14] = 56000;
        Mxv[14] = -97;
        M5w[31] = 724.6;
        M5e[31] = 68000;
        M5v[31] = -97;
        M3w[29] = 724.6;
        M3e[29] = 68000;
        M3v[29] = -97;
        Mxw[15] = 724.6;
        Mxe[15] = 68000;
        Mxv[15] = -97;
        M5w[32] = 303.2;
        M5e[32] = 7400;
        M5v[32] = 2;
        Mxw[16] = 303.2;
        Mxe[16] = 7400;
        Mxv[16] = 2;
        M3w[30] = 755.97;
        M3e[30] = 0;
        M3v[30] = -97;
        M5w[33] = 506.6;
        M5e[33] = 150000;
        M5v[33] = -97;
        M3w[31] = 644.6;
        M3e[31] = 150000;
        M3v[31] = -97;
        Mxw[17] = 506.6;
        Mxe[17] = 150000;
        Mxv[17] = -97;
        M5w[34] = 532.6;
        M5e[34] = 250000;
        M5v[34] = -97;
        M3w[32] = 670.7;
        M3e[32] = 250000;
        M3v[32] = -97;
        Mxw[18] = 532.6;
        Mxe[18] = 250000;
        Mxv[18] = -97;
        M5w[35] = 632.7;
        M5e[35] = 250000;
        M5v[35] = -97;
        M3w[33] = 770.8;
        M3e[33] = 250000;
        M3v[33] = -97;
        M5w[36] = 303.2;
        M5e[36] = 7400;
        M5v[36] = 2;
        Mxw[19] = 303.2;
        Mxe[19] = 7400;
        Mxv[19] = 2;
        Mxw[20] = 269.2;
        Mxe[20] = 7300;
        Mxv[20] = 5;
        M5w[37] = 303.2;
        M5e[37] = 7400;
        M5v[37] = 2;
        M3w[34] = 303.2;
        M3e[34] = 7400;
        M3v[34] = 2;
        Mxw[21] = 303.2;
        Mxe[21] = 7400;
        Mxv[21] = 2;
        M3w[35] = 757.8;
        M3e[35] = 32000;
        M3v[35] = -97;
        M5w[38] = 722.9;
        M5e[38] = 13000;
        M5v[38] = -97;
        M3w[36] = 722.9;
        M3e[36] = 13000;
        M3v[36] = -97;
        M5w[39] = 537.6;
        M5e[39] = 83000;
        M5v[39] = -97;
        M3w[37] = 569.5;
        M3e[37] = 83000;
        M3v[37] = -97;
        M5w[40] = 816.7;
        M5e[40] = 73000;
        M5v[40] = -97;
        M3w[38] = 816.7;
        M3e[38] = 73000;
        M3v[38] = -97;
        M5w[41] = 816.7;
        M5e[41] = 83000;
        M5v[41] = -97;
        Mxw[22] = 816.7;
        Mxe[22] = 83000;
        Mxv[22] = -97;
        M5w[42] = 329.2;
        M5e[42] = 11500;
        M5v[42] = 6;
        Mxw[23] = 329.2;
        Mxe[23] = 11500;
        Mxv[23] = 6;
        Mxw[24] = 295.2;
        Mxe[24] = 10800;
        Mxv[24] = 9;
        M5w[43] = 744.1;
        M5e[43] = 73000;
        M5v[43] = -97;
        M5w[44] = 753.9;
        M5e[44] = 170000;
        M5v[44] = -97;
        M5w[45] = 862.1;
        M5e[45] = 200000;
        M5v[45] = -97;
        M5w[46] = 666.4;
        M5e[46] = 73000;
        M5v[46] = -97;
        M3w[39] = 666.4;
        M3e[39] = 73000;
        M3v[39] = -97;
        Mxw[25] = 666.4;
        Mxe[25] = 73000;
        Mxv[25] = -97;
        M5w[47] = 208.2;
        M5e[47] = 0;
        M5v[47] = -97;
        M5w[48] = 313.2;
        M5e[48] = 5000;
        M5v[48] = 0;
        Mxw[26] = 313.2;
        Mxe[26] = 5000;
        Mxv[26] = 0;
        M5w[49] = 79.9;
        M5e[49] = 0;
        M5v[49] = -97;
        M3w[40] = 79.9;
        M3e[40] = 0;
        M3v[40] = -97;
        M5w[50] = 648.7;
        M5e[50] = 74000;
        M5v[50] = -97;
        M3w[41] = 648.7;
        M3e[41] = 74000;
        M3v[41] = -97;
        Mxw[27] = 648.7;
        Mxe[27] = 74000;
        Mxv[27] = -97;
        M5w[51] = 833;
        M5e[51] = 129000;
        M5v[51] = -97;
        M3w[42] = 833;
        M3e[42] = 129000;
        M3v[42] = -97;
        Mxw[28] = 833;
        Mxe[28] = 129000;
        Mxv[28] = -97;
        M5w[52] = 695.8;
        M5e[52] = 82000;
        M5v[52] = -97;
        M3w[43] = 695.8;
        M3e[43] = 82000;
        M3v[43] = -97;
        Mxw[29] = 695.8;
        Mxe[29] = 82000;
        Mxv[29] = -97;
        M5w[53] = 138.1;
        M5e[53] = 0;
        M5v[53] = -97;
        M3w[44] = 138.1;
        M3e[44] = 0;
        M3v[44] = -97;
        Mxw[30] = 138.1;
        Mxe[30] = 0;
        Mxv[30] = -97;
        M5w[54] = 180.1;
        M5e[54] = 0;
        M5v[54] = -97;
        M5w[55] = 344.3;
        M5e[55] = 0;
        M5v[55] = -97;
        Mxw[31] = 344.3;
        Mxe[31] = 0;
        Mxv[31] = -97;
        M5w[56] = 344.3;
        M5e[56] = 0;
        M5v[56] = -97;
        Mxw[32] = 344.3;
        Mxe[32] = 0;
        Mxv[32] = -97;
        M5w[57] = 212.1;
        M5e[57] = 0;
        M5v[57] = -97;
        Mxw[33] = 212.1;
        Mxe[33] = 0;
        Mxv[33] = -97;
        M3w[45] = 304.2;
        M3e[45] = 8700;
        M3v[45] = 19;
        Mxw[34] = 270.2;
        Mxe[34] = 8400;
        Mxv[34] = 14;
        M3w[46] = 1008;
        M3e[46] = 91000;
        M3v[46] = -97;
        M5w[58] = 591.9;
        M5e[58] = 91000;
        M5v[58] = -97;
        M3w[47] = 591.9;
        M3e[47] = 91000;
        M3v[47] = -97;
        Mxw[35] = 870.9;
        Mxe[35] = 91000;
        Mxv[35] = -97;
        M5w[59] = 675.2;
        M5e[59] = 73000;
        M5v[59] = -97;
        M5w[60] = 966.1;
        M5e[60] = 73000;
        M5v[60] = -97;
        M5w[61] = 881.5;
        M5e[61] = 116000;
        M5v[61] = -97;
        M3w[48] = 881.5;
        M3e[48] = 116000;
        M3v[48] = -97;
        Mxw[36] = 881.5;
        Mxe[36] = 116000;
        Mxv[36] = -97;
        M3w[49] = 244.3;
        M3e[49] = 0;
        M3v[49] = -97;
        M5w[62] = 328.4;
        M5e[62] = 0;
        M5v[62] = -97;
        M5w[63] = 369.1;
        M5e[63] = 9900;
        M5v[63] = 20;
        Mxw[37] = 369.1;
        Mxe[37] = 9900;
        Mxv[37] = 20;
        M5w[64] = 209.2;
        M5e[64] = 0;
        M5v[64] = -97;
        Mxw[38] = 209.2;
        Mxe[38] = 0;
        Mxv[38] = -97;
        M5w[65] = 718.22;
        M5e[65] = 84000;
        M5v[65] = -97;

        //DNA        double Mw = (an * 313.21) + (tn * 304.2) + (cn * 289.18) + (gn * 329.21) + (in * 314) + (un * 290.169) - 61.96;
        // a   b   c   d   g   h  i    k   m   n   r   s   t  u   v   w    x   y  z
        //97  98  99  100 103 104 105 107 109 110 114 115 116 117 118 119 120 121 122
        //0   1    2  3   6   7   8   10  12  13  17  18  19  20  21  22  23  24  25

        const allMW = new Array(5).fill(0).map(() => new Array(30).fill(0));
        const q = 0;// 0-DNA; 1,2-RNA;
        allMW[0][0] = 313.209;  //a
        allMW[0][1] = 307.5293; //b
        allMW[0][2] = 289.184;  //c
        allMW[0][3] = 315.5377; //d
        allMW[0][6] = 329.208;  //g
        allMW[0][7] = 302.1963; //h
        allMW[0][8] = 314.194;  //i
        allMW[0][10] = 316.702; //k
        allMW[0][12] = 301.1965;//m
        allMW[0][13] = 308.9492;//n
        allMW[0][17] = 321.2085;//r
        allMW[0][18] = 309.196; //s
        allMW[0][19] = 304.196; //t
        allMW[0][20] = 306.17;  //u
        allMW[0][21] = 310.5337;//v
        allMW[0][22] = 308.7025;//w
        allMW[0][23] = 308.9492;//x
        allMW[0][24] = 296.69;  //y
        // RNA DNA+15.999;
        allMW[1][0] = 330.21;  //a
        allMW[1][1] = 307.5293 + 15.999; //b
        allMW[1][2] = 305.18;  //c
        allMW[1][3] = 315.5377 + 15.999; //d
        allMW[1][6] = 345.21;  //g
        allMW[1][7] = 302.1963 + 15.999; //h
        allMW[1][8] = 314.194 + 15.999;  //i
        allMW[1][10] = 316.702 + 15.999; //k
        allMW[1][12] = 301.1965 + 15.999;//m
        allMW[1][13] = 308.9492 + 15.999;//n
        allMW[1][17] = 321.2085 + 15.999;//r
        allMW[1][18] = 309.196 + 15.999; //s
        allMW[1][19] = 306.17; //t
        allMW[1][20] = 306.17; //u
        allMW[1][21] = 310.5337 + 15.999;//v
        allMW[1][22] = 308.7025 + 15.999;//w
        allMW[1][23] = 308.9492 + 15.999;//x
        allMW[1][24] = 296.69 + 15.999;  //y
        // methyl-RNA: DNA+30.026;
        allMW[2][0] = 313.209 + 30.026;  //a
        allMW[2][1] = 307.5293 + 30.026; //b
        allMW[2][2] = 289.184 + 30.026;  //c
        allMW[2][3] = 315.5377 + 30.026; //d
        allMW[2][6] = 329.208 + 30.026;  //g
        allMW[2][7] = 302.1963 + 30.026; //h
        allMW[2][8] = 314.194 + 30.026;  //i
        allMW[2][10] = 316.702 + 30.026; //k
        allMW[2][12] = 301.1965 + 30.026;//m
        allMW[2][13] = 308.9492 + 30.026;//n
        allMW[2][17] = 321.2085 + 30.026;//r
        allMW[2][18] = 309.196 + 30.026; //s
        allMW[2][19] = 306.17; //t
        allMW[2][20] = 306.17; //u
        allMW[2][21] = 310.5337 + 30.026;//v
        allMW[2][22] = 308.7025 + 30.026;//w
        allMW[2][23] = 308.9492 + 30.026;//x
        allMW[2][24] = 296.69 + 30.026;  //y
        //LNA: DNA+28.011;
        allMW[3][0] = 341.22;  //a
        allMW[3][1] = 307.5293 + 28.011; //b
        allMW[3][2] = 289.184 + 28.011;  //c
        allMW[3][3] = 315.5377 + 28.011; //d
        allMW[3][6] = 329.208 + 28.011;  //g
        allMW[3][7] = 302.1963 + 28.011; //h
        allMW[3][8] = 314.194 + 28.011;  //i
        allMW[3][10] = 316.702 + 28.011; //k
        allMW[3][12] = 301.1965 + 28.011;//m
        allMW[3][13] = 308.9492 + 28.011;//n
        allMW[3][17] = 321.2085 + 28.011;//r
        allMW[3][18] = 309.196 + 28.011; //s
        allMW[3][19] = 304.196 + 28.011; //t
        allMW[3][20] = 290.169 + 28.011; //u
        allMW[3][21] = 310.5337 + 28.011;//v
        allMW[3][22] = 308.7025 + 28.011;//w
        allMW[3][23] = 308.9492 + 28.011;//x
        allMW[3][24] = 296.69 + 28.011;  //y
        //Phos DNA+16.062
        allMW[4][0] = 313.209 + 16.061;  //a
        allMW[4][1] = 307.5293 + 16.061; //b
        allMW[4][2] = 289.184 + 16.061;  //c
        allMW[4][3] = 315.5377 + 16.061; //d
        allMW[4][6] = 329.208 + 16.061;  //g
        allMW[4][7] = 302.1963 + 16.061; //h
        allMW[4][8] = 314.194 + 16.061;  //i
        allMW[4][10] = 316.702 + 16.061; //k
        allMW[4][12] = 301.1965 + 16.061;//m
        allMW[4][13] = 308.9492 + 16.061;//n
        allMW[4][17] = 321.2085 + 16.061;//r
        allMW[4][18] = 309.196 + 16.061; //s
        allMW[4][19] = 304.196 + 16.061; //t
        allMW[4][20] = 290.169 + 16.061; //u
        allMW[4][21] = 310.5337 + 16.061;//v
        allMW[4][22] = 308.7025 + 16.061;//w
        allMW[4][23] = 308.9492 + 16.061;//x
        allMW[4][24] = 296.69 + 16.061;  //y

        let salt_M = 0.055;
        let Mg_M = 0.001;
        let pr_mkM = 0.2;
        let Mg0 = NumToDouble(document.getElementById('mg_concentration').value);
        let salt0 = NumToDouble(document.getElementById('salt_concentration').value);
        let p0 = NumToDouble(document.getElementById('primer_concentration').value);
        let min3 = parseInt(document.getElementById('sensivity').value);

        if (salt0 > 0) {
                salt_M = salt0 / 1000;
        }
        if (Mg0 > 0) {
                Mg_M = Mg0 / 1000;
        }
        if (p0 > 0) {
                pr_mkM = p0;// /1000000;
        }
        const KMg_M = salt_M + (3.795 * Math.sqrt(Mg_M));

        if (w1 === undefined) { w1 = 0; }

        let Mw = - 61.962;
        let M = 0;
        let od = 1;
        let nmol = 0;
        let mg = 0;
        let od0 = NumToDouble(document.getElementById('od').value);
        let mg0 = NumToDouble(document.getElementById('mg').value);
        let nmol0 = NumToDouble(document.getElementById('nmol').value);
        let stock_concentration = NumToDouble(document.getElementById('stock_concentration').value);

        const sq = PrimerSeq(document.getElementById('inputText').value.toLowerCase().trim());

        // sq="cgcggggaatcgcgcgc";

        const lseq = sq.length;

        const result = ["", ""];
        let resultarea = "5'-" + sq + "-3'\n\n";
        let an = 0, tn = 0, cn = 0, gn = 0, inosine = 0, un = 0;

        for (var i = 0; i < lseq; i++) {
                var chr = sq.charCodeAt(i);
                if (chr === 109) { //m
                        an += 0.5;
                        cn += 0.5;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 114) { //r
                        an += 0.5;
                        gn += 0.5;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 119) {//'w'
                        an += 0.5;
                        tn += 0.5;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 115) {//s
                        gn += 0.5;
                        cn += 0.5;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 121) {//y
                        tn += 0.5;
                        cn += 0.5;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 107) {//k
                        gn += 0.5;
                        tn += 0.5;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 118) {//v
                        an += 0.333;
                        cn += 0.333;
                        gn += 0.333;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 104) {//h
                        an += 0.333;
                        cn += 0.333;
                        tn += 0.333;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 98) {//b
                        tn += 0.333;
                        cn += 0.333;
                        gn += 0.333;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 100) {//d
                        an += 0.333;
                        tn += 0.333;
                        gn += 0.333;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 110) {//n
                        an += 0.25;
                        cn += 0.25;
                        gn += 0.25;
                        tn += 0.25;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 97) {//a
                        an++;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 116) {//t
                        tn++;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 99) {//c
                        cn++;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 103) {//g
                        gn++;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 105) {//i
                        inosine++;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 117) {//u
                        un++;
                        Mw += allMW[q][chr - 97];
                }
                if (chr === 101) {//LNA dA e
                        en++;
                        Mw += allMW[3][0];
                }
                if (chr === 102) {//LNA dC f
                        fn++;
                        Mw += allMW[3][2];
                }
                if (chr === 106) {//LNA dG j
                        jn++;
                        Mw += allMW[3][6];
                }
                if (chr === 108) {//LNA dT l
                        ln++;
                        Mw += allMW[3][19];
                }
                if (chr === 120) { //x
                        xn++;
                        Mw += Mxw[vMx];
                        M += Mxe[vMx];
                }
                if (chr === 122) { //z
                        zn++;
                        Mw += Mxw[vMz];
                        M += Mxe[vMz];
                }
        }

        if (lseq > 1) {
                const gc1 = CG(sq);
                M += e260(sq, q);
                if (M > 0) {
                        od = 1;
                        if (w1 === 0) {
                                od = od0;
                        }
                        if (od < 0) {
                                od = 1;
                        }
                        if (w1 === 0) {
                                mg = (od * Mw * 1000) / M;
                                nmol = (od * 1000000) / M;
                                document.getElementById('mg').value = NumberToSeq(mg, 2);
                                document.getElementById('nmol').value = NumberToSeq(nmol, 2);
                        }
                        if (w1 === 1) {
                                mg = mg0;
                        }
                        if (mg < 0) {
                                mg = 0;
                        }
                        if (w1 === 1) {
                                od = (mg * M) / (Mw * 1000);
                                nmol = (mg * 1000) / Mw;
                                document.getElementById('od').value = NumberToSeq(od, 2);
                                document.getElementById('nmol').value = NumberToSeq(nmol, 2);
                        }
                        if (w1 === 2) {
                                nmol = nmol0;
                        }
                        if (nmol < 0) {
                                nmol = 0;
                        }
                        if (w1 === 2) {
                                od = (nmol * M) / 1000000;
                                mg = (nmol * Mw) / 1000;
                                document.getElementById('mg').value = NumberToSeq(mg, 2);
                                document.getElementById('od').value = NumberToSeq(od, 2);
                        }
                }
                resultarea += "Length=" + lseq + " A=" + NumberToSeq(an) + " T=" + NumberToSeq(tn) + " C=" + NumberToSeq(cn) + " G=" + NumberToSeq(gn) + " CG=" + NumberToSeq(gc1) + "%\n"
                resultarea += "Linguistic Complexity = " + LingComplexity(sq) + "%\n\n"
                resultarea += "Tm=" + NumberToSeq(Tm(sq, salt_M, Mg_M, pr_mkM)) + "°C (Allawi's thermodynamics parameters (Biochemistry,1997, 36:10581-10594)\n"
                resultarea += "Tm=" + NumberToSeq(Tm77(sq, KMg_M, 0)) + "°C (Tm = 77.1 + 11.7Log[K+] + (41*(G + C) - 528)/L)\n\n"
                resultarea += "Extinction coefficient = " + NumberToSeq(M) + " L/(mol·cm)\n"
                resultarea += "Molecular weight = " + NumberToSeq(Mw) + " g/mol\n\n"
                resultarea += "OD260 = " + NumberToSeq(od, 3) + "\n";
                resultarea += "mg = " + NumberToSeq(mg, 3) + "\n";
                resultarea += "nmol = " + NumberToSeq(nmol, 3) + "\n\n\n"
                const dilution = (nmol * 1000 * od) / stock_concentration;
                resultarea += stock_concentration + " µM = dissolve in " + NumberToSeq(dilution, 2) + " µl of MQ-water or TE buffer\n\n";
        }
        result[0] += resultarea;


        // Dimer showing
        resultarea = "\n\n";
        let x0 = [];
        let x1 = [];
        let x2 = [];
        let x3 = [];
        k = DimersLook(sq, sq, min3, min3, x0, x1, x2, x3);
        //k = DimerLookX(sq, sq, min3, min3, x0, x1, x2, x3);
        if (k > -1) {
                if (k === 0) { resultarea += "Dimer:\n"; }
                else { resultarea += (k + 1) + " dimers:\n"; }
                for (let w = 0; w <= k; w++) {
                        let d = new Array(4).fill(0);
                        d = Tmelting2(x1[w], x3[w], pr_mkM, salt_M, Mg_M);
                        if (x0[w] === 1) {
                                resultarea += "\n3'dimer: Tm=" + NumberToSeq(d[0], 1) + "°C; dG=" + NumberToSeq(d[1], 1) + " kcal/mol\n";
                        } else {
                                resultarea += "\nTm=" + NumberToSeq(d[0], 1) + "°C; dG=" + NumberToSeq(d[1], 1) + " kcal/mol\n";
                        }
                        resultarea += x1[w] + "\n";
                        resultarea += x2[w] + "\n";
                        resultarea += x3[w] + "\n\n";
                }
        }
        result[1] = resultarea;
        return result;
}

function DilutionDisplay() {
        let ms = NumToDouble(document.getElementById('stock_concentration').value);
        let v = NumToDouble(document.getElementById('target_volume').value);
        let c = NumToDouble(document.getElementById('target_concentration').value);
        let te = 0;
        let vc = 0;
        if (ms > 0) {
                vc = v * c / ms;
        }
        te = v - vc;
        if (te < 0) {
                te = 0;
        }
        document.getElementById('volume_stock').value = NumberToSeq(vc, 1);
        document.getElementById('volume_te').value = NumberToSeq(te, 1);
}

function DilutionDisplay2() {
        let mc = NumToDouble(document.getElementById('stock_concentration').value);
        let c = NumToDouble(document.getElementById('target_concentration').value);
        let vc = NumToDouble(document.getElementById('volume_stock').value);
        if (vc > 0 && c > 0 && mc > 0) {
                let te = 0;
                let v = 0;
                v = (mc * vc) / c;
                document.getElementById('target_volume').value = NumberToSeq(v, 1);
                te = v - vc;
                document.getElementById('volume_te').value = NumberToSeq(te, 1);
        }
}

function DilutionDisplay3() {
        let mc = NumToDouble(document.getElementById('stock_concentration').value);
        let c = NumToDouble(document.getElementById('target_concentration').value);
        let te = NumToDouble(document.getElementById('volume_te').value);
        if (te > 0 && c > 0 && mc > 0) {
                let vc = te / ((mc - c) - 1);
                let v = vc + te;
                let r = te / ((mc / c) - 1)
                document.getElementById('target_volume').value = NumberToSeq((r + te), 1);
                document.getElementById('volume_stock').value = NumberToSeq(r,);
        }
}