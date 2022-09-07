library(testthat)

test_that("Check that there is GCbox in MYC promoters", {
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(Biostrings)
  Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  MYC <- GRanges(seqnames="chr8",
                 ranges=IRanges(start=127736231,end=127742951,names="MYC"),
                 strand="+")
  MYC_promoter <- trim(promoters(MYC, upstream=2000, downstream=200))
  MYC_promoter <- GRanges(MYC_promoter, seqinfo=seqinfo(Hsapiens))
  GCbox <- DNAStringSet("CCGCCC")
  result <- 1
  names(result) <- "CCGCCC"
  expect_equal(SeqCounter(Hsapiens,MYC_promoter,GCbox),result)
})

test_that("Test error for one of the arguments missing", {
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(Biostrings)
  Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  MYC <- GRanges(seqnames="chr8",
                 ranges=IRanges(start=127736231,end=127742951,names="MYC"),
                 strand="+")
  MYC_promoter <- trim(promoters(MYC, upstream=2000, downstream=200))
  MYC_promoter <- GRanges(MYC_promoter, seqinfo=seqinfo(Hsapiens))
  GCbox <- DNAStringSet("")
  expect_error(SeqCounter(Hsapiens,MYC_promoter,GCbox))
})

test_that("Test error for class wrong of one argument", {
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(Biostrings)
  Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  MYC <- GRanges(seqnames="chr8",
                 ranges=IRanges(start=127736231,end=127742951,names="MYC"),
                 strand="+")
  MYC_promoter <- trim(promoters(MYC, upstream=2000, downstream=200))
  MYC_promoter <- GRanges(MYC_promoter, seqinfo=seqinfo(Hsapiens))
  GCbox <- "CCGCCC"
  expect_error(SeqCounter(Hsapiens,MYC_promoter,GCbox))
})


