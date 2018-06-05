## ----setupEnv, include=FALSE---------------------------------------------
working.folder = "./GencodeReference"

while( !dir.exists(working.folder) )
   working.folder = readline(prompt="Working folder not found! Please enter a valid path: ")
if( substr(working.folder, nchar(working.folder), nchar(working.folder)) != "/")
   working.folder = paste0(working.folder, "/")

setwd(working.folder)

species = "Mus musculus" # or 'Mus musculus' for mouse
genome = "mm10" # or mm10 for mouse
version = "M14" # or M16 for mouse

## ----sessionInfo, echo=FALSE---------------------------------------------
cat(sessionInfo()$R.version$version.string, fill=TRUE)
cat(paste("Platform", sessionInfo()$platform), fill=TRUE)
cat(paste("Running under", sessionInfo()$running), fill=TRUE)
cat(paste("Last knitted on", date()), fill=TRUE)
cat(paste("Working directory set to", getwd()), fill=TRUE)
cat(paste("Species: ", species), fill=TRUE)
cat(paste("Genome assembly:", genome), fill=TRUE)
cat(paste("Gencode version:", version), fill=TRUE)

## ----loadPackages, message=FALSE, results="hide"-------------------------
require("AnnotationHub")
require("biomaRt")
require("GenomicFeatures")
require("rtracklayer")
require("data.table")
require("knitr")
require("ggplot2")
require("ggpubr")

## ----downloadGencode, results="hide", message=FALSE, warning=FALSE-------
if( species %in% "Homo sapiens" ){
   if( !file.exists(paste0("gencode.v", version, ".annotation.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_", version, "/gencode.v", version, ".annotation.gff3.gz"),
                    paste0("gencode.v", version, ".annotation.gff3.gz"), quiet = TRUE)

   if( !file.exists(paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_", version, "/gencode.v", version, ".long_noncoding_RNAs.gff3.gz"),
                    paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz"), quiet = TRUE)

   if( !file.exists(paste0("gencode.v", version, ".tRNAs.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_", version, "/gencode.v", version, ".tRNAs.gff3.gz"),
                    paste0("gencode.v", version, ".tRNAs.gff3.gz"), quiet = TRUE)

   if( !file.exists(paste0(genome, ".chrom.sizes")) )
      download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/bigZips/", genome, ".chrom.sizes"),
                    paste0(genome, ".chrom.sizes"), quiet = TRUE)
}else{
   if( !file.exists(paste0("gencode.v", version, ".annotation.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_", version, "/gencode.v", version, ".annotation.gff3.gz"),
                    paste0("gencode.v", version, ".annotation.gff3.gz"), quiet = TRUE)

   if( !file.exists(paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_", version, "/gencode.v", version, ".long_noncoding_RNAs.gff3.gz"),
                    paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz"), quiet = TRUE)

   if( !file.exists(paste0("gencode.v", version, ".tRNAs.gff3.gz")) )
      download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_", version, "/gencode.v", version, ".tRNAs.gff3.gz"),
                    paste0("gencode.v", version, ".tRNAs.gff3.gz"), quiet = TRUE)

   if( !file.exists(paste0(genome, ".chrom.sizes")) )
      download.file(paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/bigZips/", genome, ".chrom.sizes"),
                    paste0(genome, ".chrom.sizes"), quiet = TRUE)
}

## ----readAnnotation, message=FALSE, warning=FALSE------------------------
chromosomes = genomeStyles()[[gsub(" ", "_", species)]]$UCSC
chrominfo = sortSeqlevels(with(fread(paste0(genome, ".chrom.sizes"))[V1%in%chromosomes,],
                               Seqinfo(V1, seqlengths = V2, genome = genome)))

TxDb = makeTxDbFromGFF(file = paste0("gencode.v", version, ".annotation.gff3.gz"),
                       format = "gff3",
                       dataSource = paste("Gencode version", version),
                       organism = species,
                       chrominfo = chrominfo)

gff = import.gff3(paste0("gencode.v", version, ".annotation.gff3.gz"))
seqinfo(gff) = chrominfo

# Remove the duplicated genes on the chromosome Y
x.genes = unique(gff[seqnames(gff)%in%"chrX"]$gene_id)
y.genes = unique(gff[seqnames(gff)%in%"chrY"]$gene_id)
y.genes = y.genes[y.genes%in%x.genes]
gff = gff[!(seqnames(gff)%in%"chrY" & gff$gene_id%in%y.genes)]

# # Remove the ribosomal RNA genes
# ribo.gff = gff[gff$gene_type%in%"rRNA",]
# gff = gff[!gff$gene_type%in%"rRNA",]

lnc.gff = import.gff3(paste0("gencode.v", version, ".long_noncoding_RNAs.gff3.gz"))
seqinfo(lnc.gff) = chrominfo

# Remove the duplicated genes on the chromosome Y
x.genes = unique(lnc.gff[seqnames(lnc.gff)%in%"chrX"]$gene_id)
y.genes = unique(lnc.gff[seqnames(lnc.gff)%in%"chrY"]$gene_id)
y.genes = y.genes[y.genes%in%x.genes]
lnc.gff = lnc.gff[!(seqnames(lnc.gff)%in%"chrY" & lnc.gff$gene_id%in%y.genes)]

## ----extractMetadata-----------------------------------------------------
# Genes, transcripts and exons
genes = genes(TxDb)
genes = genes[!grepl("PAR", names(genes)),]
txs = transcriptsBy(TxDb, by="gene")
txs = txs[!grepl("PAR", names(txs)),]
exons = exonsBy(TxDb, by="tx", use.names=TRUE)
exons = exons[!grepl("PAR", names(exons)),]

# Gene biotypes
bio.names = names(sort(table(gff[gff$type%in%"gene"]$gene_type), decreasing=TRUE))
geneBT = lapply(bio.names,
                function(x)
                   unique(gff$gene_id[gff$gene_type%in%x]) )
names(geneBT) = bio.names

gene.metadata = data.table(as.data.frame(gff))[type%in%"gene",
                                               list(gene_id, gene_type, gene_name, level, seqnames, start, end, width, strand)]

# Transcripts biotypes
bio.names = names(sort(table(gff[gff$type%in%"transcript"]$transcript_type), decreasing=TRUE))
txsBT =  lapply(bio.names,
                function(x) unique(gff$transcript_id[gff$transcript_type%in%x]) )
names(txsBT) = bio.names

txs.metadata = data.table(as.data.frame(gff))[type%in%"transcript",
                                              list(transcript_id, transcript_type, gene_id, gene_type, gene_name, level, seqnames, start, end, width, strand)]

save(genes, txs,
     geneBT, txsBT,
     gene.metadata, txs.metadata,
     file=paste0(genome, "_Gencode", version, "_annotations.RData"))

## ----genomicRegions------------------------------------------------------
# get transcript regions (all genes)
# hierarchy: ncRNA > long_ncRNA > cds > utr3 > utr5 > intron > other > intergenic

# identify protein-coding genes, long_ncRNAs, ncRNAs and other
protein.coding = c("protein_coding", paste(rep(c("IG","TR"),each=4), c("C","D","J","V"), "gene", sep="_"), "IG_LV_gene")
long.ncRNA = unique(names(sort(-table(lnc.gff$gene_type))))
ncRNA = c(setdiff(grep("RNA", names(geneBT), value=TRUE), long.ncRNA), "ribozyme")
other = setdiff(names(geneBT), c(ncRNA, long.ncRNA, protein.coding))

protein.coding = as.vector(unlist(geneBT[protein.coding]))
long.ncRNA = as.vector(unlist(geneBT[long.ncRNA]))
ncRNA = as.vector(unlist(geneBT[ncRNA]))
other = as.vector(unlist(geneBT[other]))

gff.dt = data.table(as.data.frame(sort.GenomicRanges(gff)))
gff.dt = gff.dt[, list(seqnames, start, end, strand, type, gene_id, gene_type, gene_name)]

gff.dt[, splitFactor := factor(as.character(NA), levels=c("ncRNA", "long_ncRNA", "protein_coding", "other"))]
gff.dt[gene_id%in%protein.coding, splitFactor := "protein_coding"]
gff.dt[gene_id%in%long.ncRNA, splitFactor := "long_ncRNA"]
gff.dt[gene_id%in%ncRNA, splitFactor := "ncRNA"]
gff.dt[gene_id%in%other, splitFactor := "other"]

gff.dt = gff.dt[order(splitFactor)]

sel = c("gene", "CDS", "three_prime_UTR", "five_prime_UTR", "exon")
gff.dt = gff.dt[type%in%sel,][, `:=`(type = factor(type, levels=sel),
                                    gene_id = factor(gene_id, levels=unique(gene_id)))]

gff.dt = gff.dt[order(gene_id, type)]
gff.dt[, index := 1:nrow(gff.dt)]

data.table_to_GRanges <- function(x){
   gr = with(x, GRanges(seqnames, IRanges(start, end), strand, type, gene_id, gene_type, gene_name))
   return(gr)
}

genes = data.table_to_GRanges(gff.dt[type%in%"gene"])

exons = data.table_to_GRanges(gff.dt[type%in%"exon"])
subregions = data.table_to_GRanges(gff.dt[type%in%c("CDS", "three_prime_UTR", "five_prime_UTR")])
exons = c(subregions, exons[overlapsAny(exons, setdiff(exons, subregions)),])
temp = disjoin(exons)
overlaps = data.table(as.data.frame(findOverlaps(temp, exons)), key="subjectHits")
temp = temp[overlaps[!duplicated(queryHits), queryHits]]
mcols(temp) = mcols(exons)[overlaps[!duplicated(queryHits), subjectHits],]
exons = temp

exons$name = paste(exons$gene_id, exons$type, sep="_")
exons.ref = unique(data.table(as.data.frame(copy(mcols(exons)))))

exons = reduce(split(exons, exons$name))
exons = unlist(exons, use.names=TRUE)

exons$name = names(exons)
exons = data.table(as.data.frame(exons, row.names=NULL), key="name")

setkey(exons.ref, name)
exons = data.table_to_GRanges(exons[exons.ref,])

introns = setdiff(genes, exons)
overlaps = data.table(as.data.frame(findOverlaps(introns, genes)), key="subjectHits")
introns = introns[overlaps[!duplicated(queryHits), queryHits]]
mcols(introns) = mcols(genes)[overlaps[!duplicated(queryHits), subjectHits],]
introns$type = "intron"

# Sanity check
length(c(exons, introns)) == length(reduce(c(exons, introns), min.gapwidth=0))
length(setdiff(genes, reduce(c(exons, introns), min.gapwidth=0))) == 0

tx.regions = data.table(as.data.frame(c(exons, introns)), key="gene_id")

gff.ref = gff.dt[type%in%"gene", list(gene_id, annotation=splitFactor, index)][order(index), index := seq_along(index)]
setkey(gff.ref, gene_id)

tx.regions = tx.regions[gff.ref, nomatch=0][order(index, annotation, gene_id, start, seqnames, end)]
tx.regions = with(tx.regions, GRanges(seqnames, IRanges(start, end), strand, gene_id, gene_type, gene_name, annotation, type))

# Sanity check
length(tx.regions) == length(reduce(tx.regions, min.gapwidth=0))
length(setdiff(genes, tx.regions)) == 0

saveRDS(tx.regions, file=paste0(genome, "_Gencode", version, "_annotations.all.genes.regions.rds"))
