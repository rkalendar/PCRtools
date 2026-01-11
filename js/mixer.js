function Calculation() {
    const getNum = (id) => {
        const raw = document.getElementById(id)?.value?.trim() ?? "";
        if (raw === "") return null;
        const x = Number(raw.replace(",", "."));
        return Number.isFinite(x) ? x : null;
    };

    const setIfEmpty = (id, value) => {
        const el = document.getElementById(id);
        if (!el) return;
        if ((el.value ?? "").trim() === "" && value != null && Number.isFinite(value)) {
            el.value = String(Math.round(value * 1e6) / 1e6);
        }
    };

    let c1 = getNum('ConcentrationStock');
    let c2 = getNum('ConcentrationSolvent');
    let v1 = getNum('VolumeStock');
    let v2 = getNum('VolumeSolvent');
    let cf = getNum('FinalConcentration');
    let vf = getNum('FinalVolume');

    if (v1 !== null && v2 !== null) vf = v1 + v2;
    if (vf !== null && v1 !== null && v2 === null) v2 = vf - v1;
    if (vf !== null && v2 !== null && v1 === null) v1 = vf - v2;
    if (vf === null && v1 !== null && v2 !== null) vf = v1 + v2;

    if (c1 === null && c2 !== null && cf !== null) {
        c1 = (cf * vf - c2 * v2) / v1;
    }
    if (c2 === null && c1 !== null && cf !== null) {
        c2 = (cf * vf - c1 * v1) / v2;
    }
    if (c1 !== null && c2 !== null && cf === null) {
        cf = (c1 * v1 + c2 * v2) / vf;
    }

    if (c1 !== null && c2 !== null && cf !== null && vf !== null && v1 === null && v2 === null) {
        v2 = Math.abs((vf * (cf - c1)) / (c2 - c1));
        v1 = vf - v2;
    }

    if (c1 !== null && c2 !== null && cf !== null && vf === null && v1 !== null && v2 === null) {
        let vt1 = Math.abs(cf - c2);
        let vt2 = Math.abs(cf - c1);
        let k = v1 / vt1;
        v2 = k * vt2;
        vf = v1 + v2;
    }

    if (c1 !== null && c2 !== null && cf !== null && vf === null && v2 !== null && v1 === null) {
        let vt1 = Math.abs(cf - c2);
        let vt2 = Math.abs(cf - c1);
        let k = v2 / vt2;
        v1 = k * vt1;
        vf = v1 + v2;
    }

    if (c1 === null && c2 === null && cf !== null && v2 !== null && v1 !== null) {
        let h = (v1 > v2) ? v1 / v2 : v2 / v1;
        let m = (v1 > v2) ? v2 : v1;
        let vt1 = v1 / m;
        let vt2 = v2 / m;
        c1 = Math.abs(cf - vt2);
        c2 = Math.abs(cf - vt1);
    }


    setIfEmpty("VolumeStock", v1);
    setIfEmpty("VolumeSolvent", v2);
    setIfEmpty("FinalVolume", vf);
    setIfEmpty("FinalConcentration", cf);
    setIfEmpty("ConcentrationStock", c1);
    setIfEmpty("ConcentrationSolvent", c2);
}