function NumToDouble(val) {
        val = val.replace(',', '.');
        val = val.replace(/[^0-9.-]/g, "");
        var out = parseFloat(val);
        if (isNaN(out))
                return 0;
        else
                return out;
}
function NumberToSeq(val, x) {
        if (x === undefined) { x = 1; }

        return Number(val).toFixed(x);
}

var fitToWidth = function (val, len) {
        var length_to = 0;
        if (typeof len === "number")
                length_to = len;
        else if (typeof len === "string")
                length_to = len.length;

        var result;
        if (length_to > val.length)
                result = val.toString() + Array(length_to - val.length + 1).join(' ');
        else
                result = val.toString();

        return result;
}