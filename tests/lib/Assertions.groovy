class Assertions {

    static void assertStep1Core(String outdir) {
        assert path("${outdir}/01_Demux").exists()
        assert path("${outdir}/02_PrimerCheck").exists()
        assert path("${outdir}/03_ITSx").exists()
        assert path("${outdir}/04_Homopolymer").exists()
        assert path("${outdir}/05_Chimera").exists()
        assert path("${outdir}/06_TagJumpFiltration").exists()
        assert path("${outdir}/07_SeqTable/Seqs.parquet").exists()
        assert path("${outdir}/08_RunSummary").exists()
        assert path("${outdir}/09_DB/SeqQualities.parquet").exists()
        assert path("${outdir}/pipeline_info/README_Step1_Methods.txt").exists()
    }

    static void assertStep2Pooled(String outdir) {
        assert path("${outdir}/04.PooledResults/OTUs.fa.gz").exists()
        assert path("${outdir}/04.PooledResults/UC_Pooled.parquet").exists()
        assert path("${outdir}/04.PooledResults/OTU_table_wide.txt.gz").exists()
        assert path("${outdir}/04.PooledResults/OTU_table_long.txt.gz").exists()
        assert path("${outdir}/pipeline_info/README_Step2_Methods.txt").exists()
    }

    // FASTA file assertions
    static int countFastaHeaders(def gzPath) {
        path(gzPath).linesGzip.count { it.startsWith(">") }
    }

    static void assertFastaCount(String gzPath, int expected) {
        def n = countFastaHeaders(gzPath)
        assert n == expected : "Expected ${expected} OTUs in ${gzPath}, found ${n}"
    }

    /*
    // using nft-fasta plugin
    // (but it reads the file into a map, which is not good for large files; and also will return wrong count if there are duplicate sequence IDs?)
    static long countFastaSequences(def fastaPath) {
        def fasta = path(fastaPath).fasta
        return fasta.size() as long
    }
    
    static void assertFastaCount(def fastaPath, long expected) {
        def n = countFastaSequences(fastaPath)
        assert n == expected : "Expected ${expected} sequences in ${fastaPath}, found ${n}"
    }
    */


    // FASTQ file assertions
    /*
    // using nft-fastq plugin - this does not work with gz-compressed files?
    static long countFastqRecords(def fastqPath) {
        def fastqFile = path(fastqPath).fastq
        return fastqFile.getNumberOfRecords() as long
    }
    */

    static long countFastqRecords(def fastqPath) {
        def p = path(fastqPath)
        def name = fastqPath.toString().toLowerCase()
        def lines = name.endsWith(".gz") ? p.linesGzip : p.lines
        long nLines = lines.count { true } as long
        assert nLines % 4 == 0 : "Invalid FASTQ file ${fastqPath}: ${nLines} lines is not divisible by 4"
        return (nLines / 4) as long
    }

    static void assertFastqCount(def fastqPath, long expected) {
        def n = countFastqRecords(fastqPath)
        assert n == expected : "Expected ${expected} sequences in ${fastqPath}, found ${n}"
    }


    // TSV file assertions
    static void assertTableHeader(String gzPath, String expectedHeader) {
        def header = path(gzPath).grepLineGzip(0)
        assert header == expectedHeader : "Unexpected OTU table header: ${header}"
    }

    static long countTxtRows(def txtPath) {
        def p = path(txtPath)
        def fileName = txtPath.toString()
        def lines = fileName.endsWith(".gz") ? p.linesGzip : p.lines
        return lines.count { true } as long
    }

    static void assertTxtRowCount(def txtPath, long expected) {
        def n = countTxtRows(txtPath)
        assert n == expected : "Expected ${expected} rows in ${txtPath}, found ${n}"
    }


}
