package org.jetbrains.bio.io

class BisulfiteBamParserTest {
    // The file is a subsample from chr22 of GSM432685.
    // Required true genome or we need to sample bam files in TestOrganismGenerator
//    private val genomeQuery = GenomeQuery("hg19", "chr22")

//    @Test(expected = IllegalStateException::class) fun unindexed() {
//        print(genomeQuery.get())
//        withResource(BisulfiteBamParserTest::class.java, "example-bsseq.bam") { resource ->
//            BisulfiteBamParser.parse(resource, genomeQuery, 0)
//        }
//    }

//    @Test fun example() {
//        withResource(BisulfiteBamParserTest::class.java, "example-bsseq.bam") { resource ->
//            SamReaderFactory.makeDefault()
//                    .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
//                    .open(resource.toFile()).use {
//                BAMIndexer.createIndex(it, resource.withExtension("bai").toFile())
//            }
//
//            BisulfiteBamParser.parse(resource, genomeQuery, 0)
//        }
//    }
}