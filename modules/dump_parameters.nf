
// Custom function to dump pipeline parameters to a TSV file


/*
 * Flatten a nested map into dot-notated key/value pairs
 * Example: [foo:[bar:1]] -> [[ "foo.bar", 1 ]]
 */
def flattenMap(Map m, String prefix = '') {
    def out = []
    m.each { k, v ->
        def key = prefix ? "${prefix}.${k}" : (k as String)
        if( v instanceof Map )
            out.addAll( flattenMap((Map)v, key) )
        else
            out << [ key, v ]
    }
    return out
}

/*
 * Make a value TSV-safe and readable
 * - Converts collections to comma-separated lists
 * - Normalizes tabs/newlines so the file stays valid TSV
 */
def tsvValue(Object v) {
    if( v == null )           return 'null'
    if( v instanceof Map )    return v.collect { kk, vv -> "${kk}:${tsvValue(vv)}" }.join(',')
    if( v instanceof Collection ) return v.collect { tsvValue(it) }.join(',')
    def s = v.toString()
    return s.replace('\t','    ')
            .replace('\r','\\r')
            .replace('\n','\\n')
}

/*
 * Main function (exported in other files)
 * Return a channel emitting one string per line:  "param\tvalue"
 */
def dumpParamsTsv() {
    assert params instanceof Map : 'params must be map-like'
    def pairs = flattenMap((Map)params).sort { a, b -> a[0] <=> b[0] }
    // If you want a header, uncomment the next line
    // def lines = ['param\tvalue'] + pairs.collect { k,v -> "${k}\t${tsvValue(v)}" }
    def lines = pairs.collect { k,v -> "${k}\t${tsvValue(v)}" }
    // Emit each line on the channel (one item per line)
    return channel.fromList(lines)   // value/queue channel factory methods are standard :contentReference[oaicite:0]{index=0}
}
