function Calculation() {
  const polname = document.getElementById('taqpol').value;
  const buffername = document.getElementById('pcrbuffer').value;
  const forwardprimername = document.getElementById('forward_primer').value;
  const reverseprimername = document.getElementById('reverse_primer').value;
  const probe1name = document.getElementById('probe1').value;
  const probe2name = document.getElementById('probe2').value;
  const probe3name = document.getElementById('probe3').value;
  const probe4name = document.getElementById('probe4').value;
  const probe5name = document.getElementById('probe5').value;

  let numberreactions = document.getElementById('number_reactions');
  numberreactions = NumToDouble(numberreactions.value);
  if (numberreactions < 1) { numberreactions = 1; }
  let reactionvolume = document.getElementById('reaction_volume');
  reactionvolume = NumToDouble(reactionvolume.value);
  if (reactionvolume < 0) { reactionvolume = 0; }
  let stock_buffer = document.getElementById('stock_buffer');
  stock_buffer = NumToDouble(stock_buffer.value);
  if (stock_buffer < 0) { stock_buffer = 0; }
  let final_buffer = document.getElementById('final_buffer');
  final_buffer = NumToDouble(final_buffer.value);
  if (final_buffer < 0) { final_buffer = 0; }
  let stock_mg = document.getElementById('stock_mg');
  stock_mg = NumToDouble(stock_mg.value);
  if (stock_mg < 0) { stock_mg = 0; }
  let final_mg = document.getElementById('final_mg');
  final_mg = NumToDouble(final_mg.value);
  if (final_mg < 0) { final_mg = 0; }
  let stock_dntp = document.getElementById('stock_dntp');
  stock_dntp = NumToDouble(stock_dntp.value);
  if (stock_dntp < 0) { stock_dntp = 0; }
  let final_dntp = document.getElementById('final_dntp');
  final_dntp = NumToDouble(final_dntp.value);
  if (final_dntp < 0) { final_dntp = 0; }
  let stock_f = document.getElementById('stock_f');
  stock_f = NumToDouble(stock_f.value);
  if (stock_f < 0) { stock_f = 0; }
  let final_f = document.getElementById('final_f');
  final_f = NumToDouble(final_f.value);
  if (final_f < 0) { final_f = 0; }
  let stock_r = document.getElementById('stock_r');
  stock_r = NumToDouble(stock_r.value);
  if (stock_r < 0) { stock_r = 0; }
  let final_r = document.getElementById('final_r');
  final_r = NumToDouble(final_r.value);
  if (final_r < 0) { final_r = 0; }
  let stock_probe1 = document.getElementById('stock_probe1');
  stock_probe1 = NumToDouble(stock_probe1.value);
  let final_probe1 = document.getElementById('final_probe1');
  final_probe1 = NumToDouble(final_probe1.value);
  let stock_probe2 = document.getElementById('stock_probe2');
  stock_probe2 = NumToDouble(stock_probe2.value);
  let final_probe2 = document.getElementById('final_probe2');
  final_probe2 = NumToDouble(final_probe2.value);
  let stock_probe3 = document.getElementById('stock_probe3');
  stock_probe3 = NumToDouble(stock_probe3.value);
  let final_probe3 = document.getElementById('final_probe3');
  final_probe3 = NumToDouble(final_probe3.value);
  let stock_probe4 = document.getElementById('stock_probe4');
  stock_probe4 = NumToDouble(stock_probe4.value);
  let final_probe4 = document.getElementById('final_probe4');
  final_probe4 = NumToDouble(final_probe4.value);
  let stock_probe5 = document.getElementById('stock_probe5');
  stock_probe5 = NumToDouble(stock_probe5.value);
  let final_probe5 = document.getElementById('final_probe5');
  final_probe5 = NumToDouble(final_probe5.value);
  let stock_dmso = document.getElementById('stock_dmso');
  stock_dmso = NumToDouble(stock_dmso.value);
  if (stock_dmso < 0) { stock_dmso = 0; }
  let final_dmso = document.getElementById('final_dmso');
  final_dmso = NumToDouble(final_dmso.value);
  if (final_dmso < 0) { final_dmso = 0; }
  let stock_rox = document.getElementById('stock_rox');
  stock_rox = NumToDouble(stock_rox.value);
  if (stock_rox < 0) { stock_rox = 0; }
  let final_rox = document.getElementById('final_rox')
  final_rox = NumToDouble(final_rox.value);
  if (final_rox < 0) { final_rox = 0; }
  let stock_sybr = document.getElementById('stock_sybr');
  stock_sybr = NumToDouble(stock_sybr.value);
  if (stock_sybr < 0) { stock_sybr = 0; }
  let final_sybr = document.getElementById('final_sybr');
  final_sybr = NumToDouble(final_sybr.value);
  if (final_sybr < 0) { final_sybr = 0; }
  let stock_taq = document.getElementById('stock_taq');
  stock_taq = NumToDouble(stock_taq.value);
  if (stock_taq < 0) { stock_taq = 0; }
  let final_taq = document.getElementById('final_taq');
  final_taq = NumToDouble(final_taq.value);
  if (final_taq < 0) { final_taq = 0; }
  let stock_pfu = document.getElementById('stock_pfu');
  stock_pfu = NumToDouble(stock_pfu.value);
  if (stock_pfu < 0) { stock_pfu = 0; }
  let final_pfu = document.getElementById('final_pfu');
  final_pfu = NumToDouble(final_pfu.value);
  if (final_pfu < 0) { final_pfu = 0; }
  let final_dna = document.getElementById('final_dna');
  final_dna = NumToDouble(final_dna.value);
  if (final_dna < 0) { final_dna = 0; }
  if (stock_probe1 < 0) { stock_probe1 = 0; }
  if (final_probe1 < 0) { final_probe1 = 0; }
  if (stock_probe2 < 0) { stock_probe2 = 0; }
  if (final_probe2 < 0) { final_probe2 = 0; }
  if (stock_probe2 < 0) { stock_probe2 = 0; }
  if (final_probe3 < 0) { final_probe3 = 0; }
  if (stock_probe3 < 0) { stock_probe3 = 0; }
  if (final_probe4 < 0) { final_probe4 = 0; }
  if (stock_probe4 < 0) { stock_probe4 = 0; }
  if (final_probe5 < 0) { final_probe5 = 0; }
  if (stock_probe5 < 0) { stock_probe5 = 0; }

  let totalvolume = numberreactions * reactionvolume;
  let result = "";
  let saves = "";
  result = "Number of reactions:\t" + numberreactions + "\n";
  result += "Reaction volume (µl):\t" + (reactionvolume.toFixed(1)) + "\n";
  result += "Total premix volume (µl):\t" + (totalvolume.toFixed(1)) + "\n\n";
  document.getElementById('total_volume').value = (totalvolume.toFixed(1));
  result += "\tPremix solution (µl):\n";

  let buffer = 0;
  let mg = 0;
  let dntp = 0;
  let prf = 0;
  let prr = 0;
  let probe1 = 0;
  let probe2 = 0;
  let probe3 = 0;
  let probe4 = 0;
  let probe5 = 0;
  let dmso = 0;
  let rox = 0;
  let sybr = 0;
  let taq = 0;
  let pfu = 0;
  let dna = 0;

  if (stock_buffer > 0 && final_buffer > 0) {
    buffer = (totalvolume * final_buffer) / stock_buffer
    saves += stock_buffer + "X " + buffername + "\t" + buffer.toFixed(2) + "\n";
  }
  if (stock_mg > 0 && final_mg > 0) {
    mg = (totalvolume * final_mg) / stock_mg
    saves += stock_mg + " mM MgCl2 or MgSO4\t" + mg.toFixed(2) + "\n";
  }
  if (stock_dntp > 0 && final_dntp > 0) {
    dntp = (totalvolume * final_dntp) / stock_dntp
    saves += stock_dntp + " mM dNTP\t" + dntp.toFixed(2) + "\n";
  }
  if (stock_f > 0 && final_f > 0) {
    prf = (totalvolume * final_f) / stock_f
    saves += stock_f + " µM " + forwardprimername + "\t" + prf.toFixed(2) + "\n";
  }
  if (stock_r > 0 && final_r > 0) {
    prr = (totalvolume * final_r) / stock_r
    saves += stock_r + " µM " + reverseprimername + "\t" + prr.toFixed(2) + "\n";
  }
  if (stock_probe1 > 0 && final_probe1 > 0) {
    probe1 = (totalvolume * final_probe1) / stock_probe1
    saves += stock_probe1 + " µM " + probe1name + "\t" + probe1.toFixed(2) + "\n";
  }
  if (stock_probe2 > 0 && final_probe2 > 0) {
    probe2 = (totalvolume * final_probe1) / stock_probe2
    saves += stock_probe2 + " µM " + probe2name + "\t" + probe2.toFixed(2) + "\n";
  }
  if (stock_probe3 > 0 && final_probe3 > 0) {
    probe3 = (totalvolume * final_probe3) / stock_probe3
    saves += stock_probe3 + " µM " + probe3name + "\t" + probe3.toFixed(2) + "\n";
  }
  if (stock_probe4 > 0 && final_probe4 > 0) {
    probe4 = (totalvolume * final_probe4) / stock_probe4
    saves += stock_probe4 + " µM " + probe4name + "\t" + probe4.toFixed(2) + "\n";
  }
  if (stock_probe5 > 0 && final_probe5 > 0) {
    probe5 = (totalvolume * final_probe5) / stock_probe5
    saves += stock_probe5 + " µM " + probe5name + "\t" + probe5.toFixed(2) + "\n";
  }
  if (stock_dmso > 0 && final_dmso > 0) {
    dmso = (totalvolume * final_dmso) / stock_dmso
    saves += stock_dmso + "% DMSO\t" + dmso.toFixed(2) + "\n";
  }
  if (final_rox > 0 && final_rox > 0) {
    rox = (totalvolume * final_rox) / stock_rox
    saves += stock_rox + "% ROX\t" + rox.toFixed(2) + "\n";
  }
  if (stock_sybr > 0 && final_sybr > 0) {
    sybr = (totalvolume * final_sybr) / stock_sybr
    saves += stock_sybr + "% SybrGreen\t" + sybr.toFixed(2) + "\n";
  }
  if (stock_taq > 0 && final_taq > 0) {
    taq = (totalvolume * final_taq) / stock_taq
    saves += stock_taq + " U/µl " + polname + "\t" + taq.toFixed(2) + "\n";
  }
  if (stock_pfu > 0 && final_pfu > 0) {
    pfu = (totalvolume * final_pfu) / stock_pfu
    saves += stock_pfu + " U/µl pfu-polymerase\t" + pfu.toFixed(2) + "\n";
  }
  if (final_dna > 0) { dna = numberreactions * final_dna; }
  let water = totalvolume - (buffer + mg + dntp + prf + prr + probe1 + probe2 + probe3 + probe4 + probe5 + dmso + rox + sybr + pfu + dna + taq)
  result += "Milli-Q water\t" + water.toFixed(2) + "\n";
  result += saves + "\n";
  result += (reactionvolume - final_dna).toFixed(1) + " µl of premix add " + final_dna + " µl DNA"
  return result;
}