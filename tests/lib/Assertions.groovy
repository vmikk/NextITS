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

}
