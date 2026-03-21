/**
 * MIXER-JS : Robust Concentration and Volume Calculator
 */

class ConcentrationMixer {
  /**
   * Solve for missing variables in a concentration mixing problem:
   * c1*v1 + c2*v2 = cf*vf
   * v1 + v2 = vf
   * 
   * @param {Object} inputs - Map of variable names to numbers (or null/undefined)
   * @returns {Object} - Map of all variables with calculated values
   */
  static solve(inputs) {
    let { c1, c2, v1, v2, cf, vf } = inputs;

    // Helper: Count non-null values
    const count = (arr) => arr.filter(v => v !== null && v !== undefined && Number.isFinite(v)).length;

    // Basic volume relationships: v1 + v2 = vf
    const solveVolumes = () => {
      if (v1 != null && v2 != null && vf == null) vf = v1 + v2;
      else if (v1 != null && vf != null && v2 == null) v2 = vf - v1;
      else if (v2 != null && vf != null && v1 == null) v1 = vf - v2;
    };

    // Run multiple passes to catch dependencies
    for (let pass = 0; pass < 3; pass++) {
      solveVolumes();

      // Case: Standard C1V1 + C2V2 = CfVf
      // 1. Solve for Concentrations
      if (cf == null && count([c1, v1, c2, v2, vf]) === 5 && vf !== 0) {
        cf = (c1 * v1 + c2 * v2) / vf;
      } else if (c1 == null && count([cf, vf, c2, v2, v1]) === 5 && v1 !== 0) {
        c1 = (cf * vf - c2 * v2) / v1;
      } else if (c2 == null && count([cf, vf, c1, v1, v2]) === 5 && v2 !== 0) {
        c2 = (cf * vf - c1 * v1) / v2;
      }

      // 2. Solve for Volumes using Pearson's Square / Allegation method
      // If we have all concentrations and one volume, we can find the others
      if (count([c1, c2, cf]) === 3) {
        const diff1 = Math.abs(cf - c2); // Parts of C1
        const diff2 = Math.abs(cf - c1); // Parts of C2
        const totalParts = diff1 + diff2;

        if (totalParts > 0) {
          if (vf != null && v1 == null && v2 == null) {
            v1 = vf * (diff1 / totalParts);
            v2 = vf * (diff2 / totalParts);
          } else if (v1 != null && v2 == null) {
            v2 = v1 * (diff2 / diff1);
          } else if (v2 != null && v1 == null) {
            v1 = v2 * (diff1 / diff2);
          }
        }
      }
    }

    return { c1, c2, v1, v2, cf, vf };
  }

  /**
   * Formats a number for display, avoiding floating point artifacts
   */
  static format(val) {
    if (val == null || !Number.isFinite(val)) return "";
    // Round to 6 decimal places but remove trailing zeros
    return parseFloat(val.toFixed(6)).toString();
  }
}

/**
 * Main UI Controller
 */
function Calculation() {
  const fields = [
    'ConcentrationStock', 'ConcentrationSolvent',
    'VolumeStock', 'VolumeSolvent',
    'FinalConcentration', 'FinalVolume'
  ];

  const fieldMap = {
    ConcentrationStock: 'c1',
    ConcentrationSolvent: 'c2',
    VolumeStock: 'v1',
    VolumeSolvent: 'v2',
    FinalConcentration: 'cf',
    FinalVolume: 'vf'
  };

  // 1. Extract and normalize inputs
  const inputs = {};
  fields.forEach(id => {
    const el = document.getElementById(id);
    if (!el) return;
    
    const raw = el.value.trim().replace(',', '.');
    const val = raw === "" ? null : Number(raw);
    inputs[fieldMap[id]] = (val !== null && Number.isFinite(val)) ? val : null;
  });

  // 2. Perform calculation
  const results = ConcentrationMixer.solve(inputs);

  // 3. Update empty fields in the UI
  Object.keys(fieldMap).forEach(id => {
    const el = document.getElementById(id);
    if (!el) return;

    const key = fieldMap[id];
    const currentVal = el.value.trim();
    
    // Only update if the field was empty and we have a result
    if (currentVal === "" && results[key] !== null) {
      el.value = ConcentrationMixer.format(results[key]);
      // Optional: Add a class to highlight calculated fields
      el.classList.add('calculated-highlight');
    }
  });
}

// Export for Node.js environments (testing)
if (typeof module !== 'undefined' && module.exports) {
  module.exports = { ConcentrationMixer };
}
